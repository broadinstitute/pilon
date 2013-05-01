/*
 * Copyright 2012, 2013 Broad Institute, Inc.
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
import net.sf.samtools._


object GapFiller {
  val minExtend = 20
}

class GapFiller(val region: GenomeRegion) {
  val noSolution = (0, "", "")

  def fixBreak(break: Region) = {
    if (Pilon.verbose) println("#fixBreak: " + break)
    //val reads = if (break.size < 100) recruitLocalReads(break) else recruitReads(break)
    val reads = recruitReads(break)
    var (start, right, left, stop) = assembleIntoBreak(break, reads)
    val solution = joinBreak(start, right, left, stop)
    //if (Pilon.debug) println("S:" + seq.size + ": " + seq)
    if (solution != noSolution) solution
    else if (Pilon.fixList contains 'breaks) {
      val newStart = start + right.length
      val newStop = stop - left.length
      // To be worthy, one side or the other must have extended a non-trivial amount
      if (newStart >= break.start + GapFiller.minExtend ||
          newStop <= break.stop - GapFiller.minExtend) {
        if (Pilon.debug) println((break.start, break.stop, newStart, newStop))
        val newGapLength = Pilon.minGap //max (newStop - newStart)
        val newGap = (1 to newGapLength) map {_ => "N"} mkString("")
        val seq = right + newGap + left
        if (Pilon.debug) println("S:" + seq.size + ": " + seq)
        val solution = trimPatch(start, seq, stop)
        if (solution._3 != "") solution 
        else noSolution
      } else noSolution
    } else noSolution
  }  
  
  def fillGap(gap: Region) = {
    if (Pilon.verbose) println("# Filling gap " + gap)
    val reads = if (gap.size < 100) recruitLocalReads(gap) else recruitReads(gap)
    //val reads = recruitReads(gap)
    val (start, right, left, stop) = assembleIntoBreak(gap, reads)
    val solution = joinBreak(start, right, left, stop)
    var gapOk = false
    if (solution != noSolution) {
      // Sanity check gap margin; must be within gapMargin bases
      val closedLength = solution._3.length 
      val closedDiff = (closedLength - gap.size).abs
      if (closedDiff <= Pilon.gapMargin) gapOk = true
      else if (Pilon.debug) println("Gap closed but bad size: " + gap + " " + closedLength)
    }
    if (gapOk) {
      if (Pilon.debug) println("Closed gap " + gap)
      solution
    } else {
      val newStart = start + right.length
      val newStop = stop - left.length
      if (newStart >= gap.start + GapFiller.minExtend || 
          newStop <= gap.stop - GapFiller.minExtend) {
        if (Pilon.debug) println((gap.start, gap.stop, newStart,newStop))
        val newGapLength = Pilon.minGap max (newStop - newStart)
        val newGap = (1 to newGapLength) map {_ => "N"} mkString("")
        val seq = right + newGap + left
        if (Pilon.debug) println("S:" + seq.size + ": " + seq)
        val solution = trimPatch(start, seq, stop)
        if (solution._3 != "") solution else noSolution
      } else noSolution
    }
  }
  
  def fill(region: Region) = {
    if (region.size >= 100) fillGap(region)
    else fixBreak(region)
  }

  def assembleIntoBreak(break: Region, reads: List[SAMRecord]) = {
    val assembler = new Assembler(region.minDepth)
    assembler.addReads(reads)
    if (Pilon.debug) println("assembleIntoBreak: " + break + " " + assembler)
    if (Pilon.fixList contains 'novelbreaks) assembler.novel

    val startOffset = breakRadius
    var start = (break.start - startOffset) max region.start
    var stop = (break.stop + 1 + startOffset) min region.stop
    val left = region.subString(start, break.start - start)
    val right = region.subString(break.stop + 1, stop - break.stop - 1)
    //var forward = assembler.tryForward(left)
    //var reverse = assembler.tryReverse(right)
    val (forward, reverse) = assembler.bridge(left, right)
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
    (start, forward, reverse, stop)
  }


  
  def joinBreak(startArg: Int, forward: String, reverse: String, stopArg: Int) = {
    var start = startArg
    var stop = stopArg
    var interval = stop - start
    val k = 2 * Assembler.K + 1
    //val k = Assembler.K
    
    var patch = properOverlap(forward, reverse, k)
    if (Pilon.debug) {
      println("joinBreak start=" + startArg + " stop=" + stopArg)
      //println("forward=" + forward.length + " reverse=" + reverse.length)
      //println("F:" + forward)
      //println("R:" + reverse)
      println("P:" + patch)
    }

    if (patch != "") {
    	val solution = trimPatch(start, patch, stop)
    	start = solution._1
    	
    	if (solution._2.length == 0 && solution._3.length == 0) {
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
    //for (leftOffset <- 0 until left.length if left.length - leftOffset >= minOverlap) {
    for (rightOffset <- 0 until right.length - minOverlap) {
      val leftOffset = left.lastIndexOfSlice(right.slice(rightOffset, rightOffset + minOverlap))
      if (leftOffset >= 0) {
        val ll = left.length - leftOffset
        val rl = right.length - rightOffset
        val len = ll min rl
        if (Pilon.debug) {
          println("overlap: " + leftOffset + "/" + rightOffset + " " + len) 
          println("L: " + left.slice(leftOffset, leftOffset+len))
          println("R: " + right.slice(rightOffset, rightOffset+len))
        }
        if (len >= minOverlap && 
            left.slice(leftOffset, leftOffset+len) == right.slice(rightOffset, rightOffset+len))
          return left.substring(0, leftOffset) + right.substring(rightOffset)
      }
    }
    ""
  } 
  
  def properOverlapOld(left: String, right: String, minOverlap: Int): String = {
    //for (leftOffset <- 0 until left.length if left.length - leftOffset >= minOverlap) {
    for (leftOffset <- 0 until left.length - minOverlap) {
      val rightOffset = right.indexOfSlice(left.slice(leftOffset, leftOffset + minOverlap))
      if (rightOffset >= 0) {
        val ll = left.length - leftOffset
        val rl = right.length - rightOffset
        val len = ll min rl
        if (Pilon.debug) {
          println("overlap: " + leftOffset + "/" + rightOffset + " " + len) 
          println("L: " + left.slice(leftOffset, leftOffset+len))
          println("R: " + right.slice(rightOffset, rightOffset+len))
        }
        if (len >= minOverlap && 
            left.slice(leftOffset, leftOffset+len) == right.slice(rightOffset, rightOffset+len))
          return left.substring(0, leftOffset) + right.substring(rightOffset)
      }
    }
    ""
  }
  
  def trimPatch(startArg: Int, patchArg: String, stopArg: Int) = {
    var start = startArg
    var stop = stopArg
    var patch = patchArg
    while (start < stop && patch.length > 0 && region.baseAt(start) == patch(0)) {
    	start += 1
    	patch = patch.substring(1)
    }
    while (start < stop && patch.length > 0 && region.baseAt(stop-1) == patch(patch.length-1)) {
    	stop -= 1
    	patch = patch.substring(0, patch.length-1)
    }
    (start, region.subString(start, stop-start), patch)
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
  
}
