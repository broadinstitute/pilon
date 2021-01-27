/*
 * Copyright (c) 2012-2018 Broad Institute, Inc.
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
 *
 */

package org.broadinstitute.pilon

import java.io.{File,PrintWriter,FileWriter,BufferedWriter}

import scala.collection.mutable.Map
import scala.util.Random
import scala.io.Source
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
    contigs.map(_.getBases.length.toLong).sum
  }

  val regions = {
    if (targets != "")
      parseTargets(targets)
    else
      (contigs map { c => (c.getName(), contigRegions(c)) })
  }

  // This are the input fasta elements o in
  val contigsOfInterest = (regions map {_._1}).toSet

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
  

  def writeFastaElement(writer: PrintWriter, header: String, sequence: String) = {
    writer.println(">" + header)
    sequence sliding(80, 80) foreach writer.println
  }

  def processRegions(bamFiles: List[BamFile]) = {
    println("Input genome size: " + genomeSize)

    bamFiles foreach {_.validateSeqs(contigsOfInterest)}

    if (Pilon.strays || Pilon.fixCircles) {
      println("Scanning BAMs")
      bamFiles.map(_.scan(contigsOfInterest))
    }

    val circles = Scaffold.findHgapCircles(bamFiles)

    // If assemble novel sequence up front, so that we can potentially place the
    // contigs into scaffolds when we process them.
    if (Pilon.fixNovel)
      Pilon.novelContigs = assembleNovel(bamFiles)    

    var chunks = regions.map(_._2).flatten
    chunks foreach { r =>
      println("Processing " + r)
      r.initializePileUps
      bamFiles foreach { r.processBam(_) }
      r.postProcess
      if ((Pilon.fixCircles) &&
        /*(r.contig.length < 5000 && r.contig.length >= 1000) ||*/ (circles contains r.contig.getName)) {
        if (Pilon.verbose) println(s"$r might be a circle!")
        r.closeCircle(circles.getOrElse(r.contig.getName, 0))
      }
      if (Pilon.vcf || Pilon.fixSnps || Pilon.fixIndels || Pilon.fixGaps || Pilon.fixLocal) {
        r.identifyAndFixIssues
        // If we don't need pileups for VCF later, free up the memory now!
        if (!Pilon.vcf) r.finalizePileUps
      }
      println(s"$r log:")
      r.printLog()
      println("Finished processing " + r)
    }

    val changesFile = Pilon.outputFile(".changes")
    val changesWriter = if (Pilon.changes)
                          new PrintWriter(new BufferedWriter(new FileWriter(changesFile)))
                        else null
    val fastaFile = Pilon.outputFile(".fasta")
    val fastaWriter = if (Pilon.fixSnps || Pilon.fixIndels || Pilon.fixGaps || Pilon.fixLocal || Pilon.fixNovel)
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
      if (fastaWriter != null) {
        println("Writing updated " + newName + " to " + fastaFile)
        val fixedRegions = reg._2 map { _.bases }
        val bases = fixedRegions reduceLeft {_ ++ _} map {_.toChar} mkString ""
        writeFastaElement(fastaWriter, newName, bases)
      }
    }

    if (Pilon.fixNovel) {
      val novelContigs = Pilon.novelContigs
      for (n <- 0 until novelContigs.length) {
        val header = "pilon_novel_%03d".format(n + 1)
        println("Appending " + header + " length " + novelContigs(n).length)
        writeFastaElement(fastaWriter, header, novelContigs(n))
      }
    }
    if (fastaWriter != null) fastaWriter.close
    if (Pilon.vcf) vcf.close
    if (Pilon.changes) changesWriter.close
    coverageSummary(bamFiles)
  }

  def coverageSummary(bamFiles: List[BamFile]) = {
    val bamTypes = bamFiles.groupBy(_.bamType)
    var totalBaseCount: Long = 0
    for ((t,files) <- bamTypes) {
      val typeBaseCount = files.map(_.baseCount).sum
      println("Mean " + t + " coverage: " + roundDiv(typeBaseCount, genomeSize))
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
    if (Pilon.verbose) println()
    val assembler = new Assembler()
    bamFiles filter {_.bamType != "jumps"} foreach { bam =>
      val reads = bam.getUnaligned
      print("..." + reads.length + " unmapped reads")
      assembler.addReads(reads)
    }
    print("...assembling contigs")
    val contigs = assembler.novel(genomeGraph)
    println()
    val contigLengths = contigs map {_.length}
    println("Assembled %d novel contigs containing %d bases".format(contigs.length, contigLengths.sum))
    contigs
  }

  def parseTargetString(targetString: String) = {
    val targetHelp = "Target string must be of the form \"element:start-stop\""
    val targets = for (target <- targetString.split(",")) yield {

      val (contig_name, start, stop) = try {
        // look for trailing coordinate spec
        val pattern = "^\\s*(.+):([0-9]+)-([0-9]+)\\s*$".r
        val pattern(contig_name, start, stop) = target
        (contig_name, start.toInt, stop.toInt)
      } catch {
        //use entire contig to define a region if start and stop are not specified
        case e: MatchError => (target, 1, contigMap(target).length)
      }

      val contig = contigMap(contig_name)
      val region = new GenomeRegion(contig, start, stop)

      println("Target: " + region)
      (contig.getName, List(region))
    }
    targets.toList
  }

  def parseTargetFile(fileName: String) = {
    try {
      Source.fromFile(fileName).getLines().flatMap({parseTargetString(_)}).toList
    } catch {
      case ioe: Exception => Nil
    }
  }

  def parseTargets(targetStr: String) = {
    val targetsFromFile = parseTargetFile(targetStr)
    if (targetsFromFile.isEmpty) parseTargetString(targetStr)
    else targetsFromFile
  }
  
  override def toString() = referenceFile.getCanonicalPath()
}
