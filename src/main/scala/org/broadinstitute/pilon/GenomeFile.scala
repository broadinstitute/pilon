/*
 * Copyright 2012-2014 Broad Institute, Inc.
 *
 * This file is part of Pilon.
 *
 * Pilon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2
 * as published by the Free Software Foundation.
 *
 * Pilon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Pilon.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.broadinstitute.pilon

import java.io.{File,PrintWriter,FileWriter,BufferedWriter}

import scala.collection.JavaConversions._
import scala.collection.mutable.Map
import scala.util.Random
import htsjdk.samtools.reference._
import Utils._

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
	  (contigList.reverse, contigMap)
	}

  val genomeSize: Long = {
    // count everything but Ns
    contigs.map(_.getBases.length).sum
  }

  val regions = {
    if (targets != "")
      parseTargets(targets)
    else
      (contigs map { c => (c.getName(), contigRegions(c)) })
  }

  def getSequence(contig: String): ReferenceSequence = {
    referenceSequenceFile.getSequence(contig)
  }

  def getSubsequenceAt(contig: String, start: Long, stop: Long): ReferenceSequence = {
    referenceSequenceFile.getSubsequenceAt(contig, start, stop)
  }

  def contigRegions(contig: ReferenceSequence) = {
    val maxSize = Pilon.chunkSize
    val cLength = contig.length
    val nChunks: Int = (cLength + maxSize - 1) / maxSize
    val chunkSize: Int = (cLength + nChunks - 1) / nChunks
    //println(contig.getName + ": " + cLength + " " + nChunks + " " + chunkSize + " " + (nChunks*chunkSize))
    Range(1, cLength + 1, chunkSize) map { base => new GenomeRegion(contig, base, cLength min base + chunkSize - 1) }
  }

  def processBam(bam: BamFile) = regions foreach { _._2 foreach { _.processBam(bam) } }
  
  def validateBam(bamFile: BamFile) = {
    val seqs = bamFile.getSeqNames
    require(!seqs.intersect(contigMap.keySet).isEmpty, bamFile + " doesn't match " + referenceFile)
  }
  
  def writeFastaElement(writer: PrintWriter, header: String, sequence: String) = {
    writer.println(">" + header)
    sequence sliding(80, 80) foreach writer.println
  }

  def processRegions(bamFiles: List[BamFile]) = {
    println("Input genome size: " + genomeSize)

    bamFiles foreach validateBam

    if (Pilon.strays) {
      println("Scanning BAMs")
      // Scan BAMs in parallel
      bamFiles.filter({_.bamType != 'unpaired}).par.map(_.scan)

      //if (Pilon.fixList contains 'scaffolds)
      //  for (bam <- bamFiles filter {_.bamType == 'jumps})
      //    Scaffold.analyzeStrays(bam)
    }
    
    // If assemble novel sequence up front, so that we can potentially place the
    // contigs into scaffolds when we process them.
    if (Pilon.fixList contains 'novel)
      Pilon.novelContigs = assembleNovel(bamFiles)    

    var chunks = regions.map(_._2).flatten
    if (Pilon.threads > 1) {
      // Do parallel processing randomly to even out load if all the 
      // big chunks are early in the file
      chunks = Random.shuffle(chunks)
    }
    chunks.par foreach { r =>
      println("Processing " + r)
      r.initializePileUps
      bamFiles foreach { r.processBam(_) }
      r.postProcess
      if (Pilon.vcf || !Pilon.fixList.isEmpty) {
    	r.identifyAndFixIssues
    	// If we don't need pileups for VCF later, free up the memory now!
    	if (!Pilon.vcf) r.finalizePileUps
      }
      println(r + " log:")
      r.printLog
      println("Finished processing " + r)
    }

    val changesFile = Pilon.outputFile(".changes")
    val changesWriter = if (Pilon.changes)
	                      new PrintWriter(new BufferedWriter(new FileWriter(changesFile)))
                      else null
    val fastaFile = Pilon.outputFile(".fasta")
    val fastaWriter = if (!Pilon.fixList.isEmpty)
	                      new PrintWriter(new BufferedWriter(new FileWriter(fastaFile)))
                      else null

    val vcf: Vcf = if (Pilon.vcf)
	                   new Vcf(Pilon.outputFile(".vcf"), regions.map({r => (r._1, r._2.map({_.size}).sum)}))
                   else null

    regions foreach { reg =>
      val name = reg._1
      val sep = if (name.indexOf("|") < 0) "_"
    	  else if (name(name.length-1) == '|') ""
    		  else "|"
      val newName = name + sep + "pilon"
      var offset = 0
      reg._2 foreach { r: GenomeRegion =>
        if (Pilon.vcf) {
          println("Writing " + r + " VCF to " + vcf.file)
          r.writeVcf(vcf)
          // free up memory by getting rid of pileups
          r.finalizePileUps
        }
        if (Pilon.changes) {
          println("Writing " + r + " changes to " + changesFile)
          r.writeChanges(changesWriter, newName, offset)
          offset += r.bases.size - r.size
        }
      }
      // Write the FASTA all at once rather than in chunks for formatting reasons
      if (!Pilon.fixList.isEmpty) {
        println("Writing updated " + newName + " to " + fastaFile)
        val fixedRegions = reg._2 map { _.bases }
        val bases = fixedRegions reduceLeft {_ ++ _} map {_.toChar} mkString ""
        writeFastaElement(fastaWriter, newName, bases)
      }
    }

    if (Pilon.fixList contains 'novel) {
      val novelContigs = Pilon.novelContigs
      for (n <- 0 until novelContigs.length) {
        val header = "pilon_novel_%03d".format(n + 1)
        println("Appending " + header + " length " + novelContigs(n).length)
        writeFastaElement(fastaWriter, header, novelContigs(n))
      }
    }
    if (!Pilon.fixList.isEmpty) fastaWriter.close
    if (Pilon.vcf) vcf.close
    if (Pilon.changes) changesWriter.close
    coverageSummary(bamFiles)
  }

  def coverageSummary(bamFiles: List[BamFile]) = {
    val bamTypes = bamFiles.groupBy(_.bamType)
    var totalBaseCount: Long = 0
    for ((t,files) <- bamTypes) {
      val typeBaseCount = files.map(_.baseCount).sum
      println("Mean " + t.name + " coverage: " + roundDiv(typeBaseCount, genomeSize))
      totalBaseCount += typeBaseCount
    }
    println("Mean total coverage: " + roundDiv(totalBaseCount, genomeSize))
  }

  def identifyAndFixIssues() = regions foreach { _._2 foreach { _.identifyAndFixIssues } }
  
  def assembleNovel(bamFiles: List[BamFile]) = {
    print("Assembling novel sequence:")
    val genomeGraph = new Assembler(minDepth = 1)
    print(" graphing genome")
    for (contig <- contigMap.values) {
      genomeGraph.addGraphSeq(GenomeRegion.baseString(contig.getBases))
      if (Pilon.verbose) print("..." + contig.getName)
    }
    //genomeGraph.buildGraph
    if (Pilon.verbose) println
    val assembler = new Assembler()
    bamFiles filter {_.bamType != 'jumps} foreach { bam =>
      val reads = bam.getUnaligned
      print("..." + reads.length + " unmapped reads")
      assembler.addReads(reads)
    }
    print("...assembling contigs")
    val contigs = assembler.novel(genomeGraph)
    println
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
