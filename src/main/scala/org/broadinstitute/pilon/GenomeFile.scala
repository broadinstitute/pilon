package org.broadinstitute.pilon

import java.io.{File,PrintWriter,FileWriter,BufferedWriter}
import scala.collection.JavaConversions._
import scala.collection.mutable.Map
import net.sf.picard.reference._

class GenomeFile(val referenceFile: File, val targets : String = "") {
	val referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile)
	val (contigs, contigMap) = {
	  var contigList = List[ReferenceSequence]()
	  var contigMap = Map.empty[String, ReferenceSequence]
	  var refseq = referenceSequenceFile.nextSequence()
	  while (refseq != null) {
	    contigMap(refseq.getName()) = refseq
	    contigList ::= refseq
	    refseq = referenceSequenceFile.nextSequence()
	  } 
	  (contigList reverse, contigMap)
	}

	val regions = {
	 if (targets != "") 
	   parseTargets(targets)
	 else 
	   (contigs map { c => (c.getName(), contigRegions(c)) } )
	}
	
	def getSequence(contig: String): ReferenceSequence = {
	  referenceSequenceFile.getSequence(contig)
	}
	
	def getSubsequenceAt(contig: String, start: Long, stop: Long): ReferenceSequence = {
	  referenceSequenceFile.getSubsequenceAt(contig, start, stop)
	}
	
	def contigRegions(contig: ReferenceSequence, maxSize: Int = 10000000) = {
	  val cLength = contig.length
	  val nChunks : Int = (cLength + maxSize - 1) / maxSize
	  val chunkSize : Int = (cLength + nChunks - 1) / nChunks
	  //println(contig.getName + ": " + cLength + " " + nChunks + " " + chunkSize + " " + (nChunks*chunkSize))
	  Range(1, cLength+1, chunkSize ) map { base => new GenomeRegion(contig,  base, cLength min base+chunkSize-1) }
	}
	 
	def referenceRegions(maxSize: Int = 10000000) = {
	  contigs flatMap { c => contigRegions(c, maxSize) }
	}
	
	
  def processBam(bam: BamFile) = regions foreach { _._2 foreach { _.processBam(bam) } }
  
  def validateBam(bamFile: BamFile) = {
    val seqs = bamFile.getSeqs
    require(!seqs.intersect(contigMap.keySet).isEmpty, bamFile + " doesn't match " + referenceFile)
  }
  
  def writeFastaElement(writer: PrintWriter, header: String, sequence: String) = {
    writer.println(">" + header)
    sequence sliding(80, 80) foreach writer.println
  }

  def processRegions(bamFiles: List[BamFile]) = {
    bamFiles foreach validateBam
    val fastaFile = Pilon.outputFile(".fasta")
	val fastaWriter = if (Pilon.fixList.length > 0) 
	  new PrintWriter(new BufferedWriter(new FileWriter(fastaFile)))
	else null
	
	val vcf: Vcf = if (Pilon.vcf) 
	  new Vcf(Pilon.outputFile(".vcf"), regions map {r => (r._1, r._2 map {_.size} sum)}) 
	else null
	
	regions foreach { reg =>
	  val name = reg._1
	  reg._2 foreach { r => 
	  	println("Processing " + r)
	    r.initializePileUps
	    bamFiles foreach { r.processBam(_) }
	    r.postProcess
	    r.identifyIssues
   	    if (Pilon.vcf) {
   	      println("Writing " + name + " VCF to " + vcf.file)
   	      r.writeVcf(vcf)
   	    }
	    r.finalizePileUps
	  } 
	  if (Pilon.fixList.length > 0) {
	    println("Fixing " + (Pilon.fixList map {_.name} mkString(", ")))
		val fixedRegions = reg._2 map { _.fixIssues }
        val bases = fixedRegions reduceLeft {_ ++ _} map {_.toChar} mkString ""		
	    println("Writing updated " + name + " to " + fastaFile)
        writeFastaElement(fastaWriter, name + "|pilon", bases)
	  }
	}
    if (Pilon.fixList contains 'novel) {
      val contigs = assembleNovel(bamFiles)
      for (n <- 0 until contigs.length) {
        val header = "pilon_novel_%04d".format(n+1)
        writeFastaElement(fastaWriter, header, contigs(n))
      }
    }
    if (Pilon.fixList.length > 0) fastaWriter.close
    if (Pilon.vcf) vcf.close
  }

  def identifyIssues() = regions foreach { _._2 foreach { _.identifyIssues } }
  
  def assembleNovel(bamFiles: List[BamFile]) = {
    println("Assembling novel sequence")
    val assembler = new Assembler
    bamFiles filter {_.bamType != 'jumps} foreach { bam =>
      if (Pilon.verbose) print("# " + bam + " ")
      val reads = bam.getUnaligned
      if (Pilon.verbose) print(reads.length + " reads")
      assembler.addReads(reads)
    }
    val contigs = assembler.novel
    val contigLengths = contigs map {_.length}
    println("Assembled %d novel contigs containing %d bases".format(contigs.length, contigLengths.sum))
    contigs
  }

  def parseTargets(targetString: String) = {
    val targetHelp = "Target string must be of the form \"element:start-stop\""
    val targets = for (target <- targetString.split(",")) yield {
      val t1 = target.split(":")
      require(t1.size <= 2, targetHelp)
      val contig = contigMap(t1(0))
      val region = 
        if (t1.size == 1) {
            new GenomeRegion(contig, 1, contig.length)
        } else {
        	val t2 = t1(1).split("-")
        	require(t2.size <= 2, targetHelp)
        	new GenomeRegion(contig, t2(0).toInt, t2(1).toInt)
        }
      println("Target: " + region)
      (contig.getName, List(region))
    }
    targets.toList
  }
  
  override def toString() = referenceFile.getCanonicalPath()
}
