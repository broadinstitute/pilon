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
  //val k = 2 * Assembler.K + 1
  val k = Assembler.K
}

class GapFiller(val region: GenomeRegion) {
  val noSolution = (0, "", "")

  def fixBreak(break: Region) = {
    if (Pilon.verbose) println("#fixBreak: " + break)
    //val reads = if (break.size < 100) recruitLocalReads(break) else recruitReads(break)
    val reads = recruitReads(break)
    // TODO: fix totally non-intuitive use of left and right
    var (start, rights, lefts, stop) = assembleIntoBreak(break, reads)
    //val solution = joinBreak(start, right, left, stop)
    val solutions = breakJoins(start, rights, lefts, stop)
    //if (Pilon.debug) println("S:" + seq.size + ": " + seq)
    //val solution = if (solutions.isEmpty) noSolution else solutions.head
    val solution = if (solutions.length == 1) solutions(0) else noSolution
    if (solution != noSolution) solution
    else if (Pilon.fixList contains 'breaks) {
      val left = consensusRight(lefts)
      val right = consensusLeft(rights)
      if (Pilon.debug) {
        println("consensusL: " + left)
        println("consensusR: " + right)
      }
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
    val (start, rights, lefts, stop) = assembleIntoBreak(gap, reads)
    //val solution = joinBreak(start, right, left, stop)
    val solutions = breakJoins(start, rights, lefts, stop)
    //val solution = if (solutions.isEmpty) noSolution else solutions.head
    val solution = if (solutions.length == 1) solutions(0) else noSolution
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
      val left = consensusRight(lefts)
      val right = consensusLeft(rights)
      if (Pilon.debug) {
        println("consensusL: " + left)
        println("consensusR: " + right)
      }
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

  def consensusLeft(seqs: List[String]): String = {
    val s0 = seqs.head
    if (seqs.length == 1) return s0
    val minLength = (seqs map {_.length}).min
    for (i <- 0 until minLength;
         s <- seqs.tail)
      if (s(i) != s0(i)) return s0.substring(0, i)
    return s0.substring(0, minLength)
  }

  def consensusRight(seqs: List[String]): String = {
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
    val assembler = new Assembler(region.minDepth)
    assembler.addReads(reads)
    if (Pilon.debug) println("assembleIntoBreak: " + break + " " + assembler)
    //if (Pilon.fixList contains 'novelbreaks) assembler.novel

    val startOffset = breakRadius
    var start = (break.start - startOffset) max region.start
    var stop = (break.stop + 1 + startOffset) min region.stop
    val left = region.subString(start, break.start - start)
    val right = region.subString(break.stop + 1, stop - break.stop - 1)
    //var forward = assembler.tryForward(left)
    //var reverse = assembler.tryReverse(right)
    val (forward, reverse) = assembler.multiBridge(left, right)
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
    if (Pilon.debug) {
      println("breakJoins: " + solutions.length + " solutions")
      for (s <- solutions) println("  " + s)
    }
    solutions
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
  
  def properOverlapScored(left: String, right: String, minOverlap: Int): String = {
    val matchScore = 1
    val mismatchScore = -10
    val scoreCutoff = -20
    val minScore = Assembler.K

    def substrMatch(a: String, aOffset: Int, b: String, bOffset: Int, len: Int): Int = {
      var score = 0
      //println(a.slice(aOffset, aOffset+len))
      //println(b.slice(bOffset, bOffset+len))
      for (i <- 0 until len)
        if (a(aOffset + i) != b(bOffset + i)) {
          score += mismatchScore
          if (score < scoreCutoff) return score
        } else score += matchScore
      score
    }
    val ll = left.length
    val rl = right.length
    //println("pO: " + ll + " " + rl + " " + minOverlap)
    for (overlap <- minOverlap to ll + rl - 2 * minOverlap) {
      val leftOffset = (ll - overlap) max 0
      val rightOffset = (overlap - ll) max 0
      val len = (ll - leftOffset) min (rl - rightOffset)
      //println(overlap + " " + leftOffset + " " + rightOffset + " " + len)
      val score = substrMatch(left, leftOffset, right, rightOffset, len)
      if (score > minScore) {
        if (Pilon.debug)
          println("overlap: " + leftOffset + "/" + left.length + " " +
                  rightOffset + "/" + right.length + " " + len + " " + score)
        return left.substring(0, leftOffset) + right.substring(rightOffset)
      }
    }
    ""
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

  def properOverlapBroken(left: String, right: String, minOverlap: Int): String = {
    /*
    val stopKmer = rightArg.substring(rightArg.length - Assembler.K)
    val stopKmerIndex = leftArg.lastIndexOf(stopKmer)
    val left = if (stopKmerIndex >= 0) leftArg.substring(0, stopKmerIndex + Assembler.K) else leftArg
    if (Pilon.debug && left != leftArg) println("trimmed left:\n" + leftArg + "\n" + left)
    val startKmer = left.substring(0, Assembler.K)
    val startKmerIndex = rightArg.indexOf(startKmer)
    val right = if (startKmerIndex >= 0) rightArg.substring(startKmerIndex) else rightArg
    if (Pilon.debug && right != rightArg) println("trimmed right:\n" + rightArg + "\n" + right)
    */
    def substrMatch(a: String, aOffset: Int, b: String, bOffset: Int, len: Int): Boolean = {
      for (i <- 0 until len)
        if (a(aOffset + i) != b(bOffset + i)) return false
      return true
    }
    for (rightOffset <- 0 until right.length - minOverlap;
         leftOffset <- left.length - minOverlap - 1 to 0 by -1) {
      //val leftOffset = left.lastIndexOfSlice(right.slice(rightOffset, rightOffset + minOverlap))
      val ll = left.length - leftOffset
      val rl = right.length - rightOffset
      val len = ll min rl
      if (substrMatch(left, leftOffset, right, rightOffset, len)) {
        if (Pilon.debug)
          println("overlap: " + leftOffset + "/" + left.length + " " +
            rightOffset + "/" + right.length + " " + len)
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
