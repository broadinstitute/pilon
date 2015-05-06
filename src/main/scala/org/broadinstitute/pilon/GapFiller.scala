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

import collection.JavaConversions._
import collection.mutable.Map
import java.io.File
import htsjdk.samtools._


object GapFiller {
  val minExtend = 20
  val k = 2 * Assembler.K + 1
  //val k = Assembler.K
}

class GapFiller(val region: GenomeRegion) {
  // TODO: change contract of this class to hold results in class state, not be all functional.
  val noSolution = (0, "", "")
  var tandemRepeat = ""

  def fillGap(gap: Region) = {
    if (Pilon.verbose) println("# Filling gap " + gap)
    assembleAcrossBreak(gap, true)
  }

  def fixBreak(break: Region) = {
    if (Pilon.verbose) println("#fixBreak: " + break)
    assembleAcrossBreak(break, false)
  }

  def assembleAcrossBreak(break: Region, isGap: Boolean) = {
    // TODO: ugh, this is ugly and really wants to be re-written.
    //val reads = if (break.size < 100 && isGap) recruitLocalReads(break) else recruitReads(break)
    val reads = recruitReads(break)
    var (start, pathsFromLeft, pathsFromRight, stop, loop) = assembleIntoBreak(break, reads)
    if (Pilon.verbose) println("L=%d R=%d".format(pathsFromLeft.length, pathsFromRight.length))
    tandemRepeat = loop
    val solutions = breakJoins(start, pathsFromLeft, pathsFromRight, stop)
    if (Pilon.verbose) {
      println("S=%d".format(solutions.length))
      if (solutions.length > 1) solutions.foreach(println)
    }

    val solution =
      if (solutions.length == 1 || (Pilon.multiClosure && solutions.length > 1)) solutions.head
      else noSolution
    var solutionOK = (solution != noSolution) && (loop.length == 0 || !Pilon.trSafe)
    if (solutionOK && isGap) {
      // Sanity check gap margin; must be within gapMargin bases of original gap size
      val closedLength = solution._3.length
      val closedDiff = (closedLength - break.size).abs
      if (closedDiff > Pilon.gapMargin) solutionOK = false
      else if (Pilon.debug) println("Gap closed but bad size: " + break + " " + closedLength)
    }
    if (solutionOK) solution
    else if (isGap || ((Pilon.fixList contains 'breaks) && loop.length == 0)) {
      // build partial solution using consensus from each side, opening gap if necessary
      val fromRight = consensusFromRight(pathsFromRight)
      val fromLeft = consensusFromLeft(pathsFromLeft)
      if (Pilon.debug) {
        println("consensusL[%d]: %s".format(fromLeft.length, fromLeft))
        println("consensusR[%d]: %s".format(fromRight.length, fromRight))
      }
      val newStart = start + fromLeft.length
      val newStop = stop - fromRight.length
      // To be worthy, one side or the other must have extended a non-trivial amount,
      // and we must have generated some sequence different from what was already there.
      if ((newStart >= break.start + GapFiller.minExtend ||
        newStop <= break.stop - GapFiller.minExtend) &&
        !partialMatchesReference(start, fromLeft, fromRight, stop, loop.length)) {
        if (Pilon.debug) println((break.start, break.stop, newStart, newStop))
        val newGapLength = if (isGap) (Pilon.minGap max (newStop - newStart)) else Pilon.minGap
        val newGap = (1 to newGapLength) map {_ => "N"} mkString("")
        val seq = fromLeft + newGap + fromRight
        if (Pilon.debug) println("S:" + seq.size + ": " + seq)
        val solution = trimPatch(start, seq, stop)
        if (solution._3 != "") solution else noSolution
      } else noSolution
    } else noSolution
  }

  def consensusFromLeft(seqs: List[String]): String = {
    val s0 = seqs.head
    if (seqs.length == 1) return s0
    val minLength = (seqs map {_.length}).min
    for (i <- 0 until minLength;
         s <- seqs.tail)
      if (s(i) != s0(i)) return s0.substring(0, i)
    return s0.substring(0, minLength)
  }

  def consensusFromRight(seqs: List[String]): String = {
    val s0 = seqs.head
    if (seqs.length == 1) return s0
    val minLength = (seqs map {_.length}).min
    for (i <- 0 until minLength;
         s <- seqs.tail) {
      val si = s.length - 1 - i
      val s0i = s0.length - 1 - i
      if (s(si) != s0(s0i)) return s0.substring(s0i + 1)
    }
    return s0.substring(s0.length - minLength)
  }

  def assembleIntoBreak(break: Region, reads: List[SAMRecord]) = {
    val assembler = new Assembler()
    assembler.addReads(reads)
    if (Pilon.dumpReads) writeBam(break.toString, reads)
    assembler.buildGraph
    if (Pilon.fixList.contains('novel) && Pilon.novelContigs != Nil) {
      assembler.addGraphSeqs(Pilon.novelContigs)
    } 
    if (Pilon.debug) println("assembleIntoBreak: " + break + " " + assembler)
      
    //if (Pilon.fixList contains 'novelbreaks) assembler.novel

    val startOffset = breakRadius
    var start = (break.start - startOffset) max region.start
    var stop = (break.stop + 1 + startOffset) min region.stop
    val left = region.subString(start, break.start - start)
    val right = region.subString(break.stop + 1, stop - break.stop - 1)
    //var forward = assembler.tryForward(left)
    //var reverse = assembler.tryReverse(right)
    val (forward, reverse, loop) = assembler.multiBridge(left, right)
    /*
    // if we didn't make it back to the beginning, use old boundary
    if (start + forward.length < break.start) {
      //start = break.start
      //     forward = ""
      forward = region.genomeSubString(start, break.start - start)
      if (Pilon.debug) println("reverting forward")
    }
    if (stop - reverse.length > break.stop) {
      //    stop = break.stop
      //reverse = ""
      reverse = region.genomeSubString(break.stop + 1, stop - break.stop - 1)
      if (Pilon.debug) println("reverting reverse")
    }
    */
    (start, forward, reverse, stop, loop)
  }

  def breakJoins(start: Int, forwardPaths: List[String], reversePaths: List[String], stop: Int) = {
    var solutionSet = Set[(Int, String, String)]()
    for (f <- forwardPaths;
         r <- reversePaths) {
      val s = joinBreak(start, f, r, stop)
      if (s != noSolution) {
        solutionSet += s
      }
    }
    val solutions = solutionSet.toList sortWith { (a,b) => a._3.length - a._2.length > b._3.length - b._2.length }
    var solutionLengths = Set[Int]()
    solutions map {s => solutionLengths += s._3.length - s._2.length}

    if (Pilon.debug) {
      println("breakJoins: " + solutions.length + " solutions; lengths " + solutionLengths)
      for (s <- solutions) println("  " + s)
    }
    if (solutionLengths.size == 1) {
      if (Pilon.debug) println("...all solutions of same length, using first above")
      List(solutions.head)
    } else{
      solutions
    }
  }
  
  def joinBreak(startArg: Int, forward: String, reverse: String, stopArg: Int) = {
    var start = startArg
    var stop = stopArg
    var interval = stop - start
    val k = GapFiller.k

    var patch = properOverlap(forward, reverse, k)
    if (Pilon.debug) {
      println("joinBreak start=" + startArg + " stop=" + stopArg)
      //println("forward=" + forward.length + " reverse=" + reverse.length)
      //println("F:" + forward)
      //println("R:" + reverse)
      if (patch != "") println("P:" + patch)
    }

    if (patch != "") {
    	val solution = trimPatch(start, patch, stop)
    	start = solution._1
    	
    	if (solution._2 == solution._3) {
    	  if (Pilon.debug) println("...no change")
    	  solution
    	} else {
    	  if (Pilon.debug)
    	    println("...start=" + start + " ref=(" + solution._2.length + ")" + solution._2 
    	    		+ " new=(" + solution._3.length + ")" + solution._3)
    	  solution
    	}
    } else {
    	if (Pilon.debug) println("...no break join")
    	noSolution
    }
  }
  
  def properOverlap(left: String, right: String, minOverlap: Int): String = {

    def substrMatch(a: String, aOffset: Int, b: String, bOffset: Int, len: Int): Boolean = {
      //println(a.slice(aOffset, aOffset+len))
      //println(b.slice(bOffset, bOffset+len))
      for (i <- 0 until len)
        if (a(aOffset + i) != b(bOffset + i)) return false
      return true
    }
    val ll = left.length
    val rl = right.length
    //println("pO: " + ll + " " + rl + " " + minOverlap)
    for (overlap <- minOverlap to ll + rl - 2 * minOverlap) {
      val leftOffset = (ll - overlap) max 0
      val rightOffset = (overlap - ll) max 0
      val len = (ll - leftOffset) min (rl - rightOffset)
      //println(overlap + " " + leftOffset + " " + rightOffset + " " + len)
      if (substrMatch(left, leftOffset, right, rightOffset, len)) {
        if (Pilon.debug)
          println("overlap: " + leftOffset + "/" + left.length + " " +
            rightOffset + "/" + right.length + " " + len)
        return left.substring(0, leftOffset) + right.substring(rightOffset)
      }
    }
    ""
  }

  def trimPatch(startArg: Int, patchArg: String, stopArg: Int) = {
    var start = startArg
    var stop = stopArg
    var patch = patchArg
    while (start < stop && patch.length > 0 &&
      (region.baseAt(start) == patch(0) || region.originalBaseAt(start) == patch(0))) {
    	start += 1
    	patch = patch.substring(1)
    }
    while (start < stop && patch.length > 0 &&
      (region.baseAt(stop-1) == patch(patch.length-1) || region.originalBaseAt(stop-1) == patch(patch.length-1))) {
    	stop -= 1
    	patch = patch.substring(0, patch.length-1)
    }
    (start, region.refSubString(start, stop-start), patch)
  }

  def partialMatchesReference(start: Int, fromLeft: String, fromRight: String, stop: Int, loopLength: Int) = {
    // If we have a non-closed solution, make sure it doesn't just recapitulate the sequence
    // of the genome (we're not learning anything then).

    // don't fall off either end of the region
    if ((start + fromLeft.length > region.stop) || (stop - fromRight.length < region.start)) {
      false
    } else {
      val regionLeft = region.subString(start, fromLeft.length)
      val regionRight = region.subString(stop - fromRight.length, fromRight.length)
      if (Pilon.debug) {
        if (fromLeft == regionLeft) print("left overlap!") else print("left mismatch!")
        if (fromRight == regionRight) print(" right overlap!") else print(" right mismatch!")
        println(" loop " + loopLength)
      }
      fromLeft == regionLeft && fromRight == regionRight
    }
  }

  
  def breakRadius = {
    val minRadius = 3 * Assembler.K
    val inserts = region.bamsOfType('frags) map {bam => bam.insertSizeMean /*+ bam.insertSizeSigma*/}
    val insertMean = if (inserts.length > 0) (inserts.sum / inserts.length).round.toInt else 0
    minRadius max insertMean
  }
    
  def recruitReadsOfType(reg: Region, bamType: Symbol) = {
    val bams = region.bamsOfType(bamType)
    var reads = List[SAMRecord]()
    for (b <- bams) {
      reads ++= b.recruitFlankReads(reg)
    }
    if (Pilon.debug) 
      println("# Recruiting " + bamType.name + " reads: count=" + reads.length)
    reads
  }

  //def recruitFrags(reg: Region) = mateMapOfType(reg, 'frags).values.toList
  def recruitFrags(reg: Region) = recruitReadsOfType(reg, 'frags)

  def recruitJumps(reg: Region) = {
    var reads = List[SAMRecord]()
    for (b <- region.bamsOfType('jumps)) {
      reads ++= b.recruitBadMates(reg)
    }
    if (Pilon.debug) 
      println("# Recruiting jump bad mates: count=" + reads.length)
    reads
  }

  def recruitUnpaired(reg: Region) = recruitReadsOfType(reg, 'unpaired)

  def recruitLocalReads(reg: Region) = recruitFrags(reg) ++ recruitUnpaired(reg)

  def recruitReads(reg: Region) = recruitLocalReads(reg) ++ recruitJumps(reg)

  def writeBam(fileName: String, reads: List[SAMRecord]) {
    val file = new File(fileName + ".sam")
    val header = new SAMFileHeader()
    header.addSequence(new SAMSequenceRecord(region.contig.getName, region.contig.getBases.length))
    var readGroups: List[SAMReadGroupRecord] = Nil
    for (r <- reads) {
      val rg = r.getReadGroup
      if (!readGroups.contains(rg)) {
        header.addReadGroup(rg)
        readGroups = rg :: readGroups
      }
    }
    val writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, file)
    for (r <- reads) writer.addAlignment(r)
    writer.close()
  }
}
