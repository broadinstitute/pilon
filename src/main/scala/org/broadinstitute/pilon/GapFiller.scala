package org.broadinstitute.pilon

import collection.JavaConversions._
import collection.mutable.Map
import net.sf.samtools._


object GapFiller {
  val K = 47
  val minDepth = 5
  val minGap = 10
  val minExtend = 20
}

class GapFiller(val region: GenomeRegion) {
  val noSolution = (0, "", "")

  def fixBreak(break: Region) = {
    if (Pilon.verbose) println("# Fix break " + break)
    val reads = if (break.size < 100) recruitLocalReads(break) else recruitReads(break)
    //val reads = recruitReads(break)
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
        val newGapLength = GapFiller.minGap //max (newStop - newStart)
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
    if (solution != noSolution) {
      if (Pilon.debug) println("Closed gap " + gap)
      solution
    }
    else {
      val newStart = start + right.length
      val newStop = stop - left.length
      if (newStart >= gap.start + GapFiller.minExtend || 
          newStop <= gap.stop - GapFiller.minExtend) {
        if (Pilon.debug) println((gap.start, gap.stop, newStart,newStop))
        val newGapLength = GapFiller.minGap max (newStop - newStart)
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
    val assembler = new Assembler
    assembler.addReads(reads)
    if (Pilon.debug) println("assembleIntoBreak: " + break + " " + assembler)
    if (Pilon.fixList contains 'novelbreaks) assembler.novel

    val k = assembler.K
    val startOffset = 3 * k
    var start = (break.start - startOffset) max region.start
    var stop = (break.stop + 1 + startOffset) min region.stop
    var forward = assembler.tryForward(region.genomeSubString(start, break.start - start))
    var reverse = assembler.tryReverse(region.genomeSubString(break.stop + 1, stop - break.stop - 1))
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
    	  noSolution
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
    for (leftOffset <- 0 until left.length - minOverlap) {
      val rightOffset = right.indexOfSlice(left.slice(leftOffset, leftOffset + minOverlap))
      if (rightOffset >= 0) {
        val ll = left.length - leftOffset
        val rl = right.length - rightOffset
        val len = ll min rl
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
    while (start < stop && patch.length > 0 && region.refBase(start) == patch(0)) {
    	start += 1
    	patch = patch.substring(1)
    }
    while (start < stop && patch.length > 0 && region.refBase(stop-1) == patch(patch.length-1)) {
    	stop -= 1
    	patch = patch.substring(0, patch.length-1)
    }
    (start, region.genomeSubString(start, stop-start), patch)
  }

  
  type MateMap = Map[SAMRecord, SAMRecord]

  
  
  def recruitReadsOfType(reg: Region, bamType: Symbol) = {
    val bams = region.bamsOfType(bamType)
    var reads = List[SAMRecord]()
    for (b <- bams) {
      val flank = b.maxInsertSize
      reads ++= readsInInterval(b, region.name, reg.start - flank, reg.stop + flank)
    }
    if (Pilon.debug) 
      println("# Recruiting " + bamType.name + " reads: count=" + reads.length)
    reads
  }

  def mateMapOfType(reg: Region, bamType: Symbol) = mateMap(recruitReadsOfType(reg, bamType))
  
  def recruitFrags(reg: Region) = mateMapOfType(reg, 'frags).values.toList

  def recruitJumps(reg: Region) = {
    val midpoint = region.midpoint

    val mm = mateMapOfType(reg, 'jumps)
    
    // Filter to find pairs where r1 is anchored and r2 is unmapped (we'll include r2)
    val mm2 = mm filter { pair =>  
      val (r1, r2) = pair
      (!r1.getReadUnmappedFlag) && r2.getReadUnmappedFlag
    }
    if (Pilon.debug) 
      println("# Filtered jumps " + mm2.size + "/" + mm.size)
    mm2.values.toList
  }
  
  def recruitUnpaired(reg: Region) = recruitReadsOfType(reg, 'unpaired)

  def recruitLocalReads(reg: Region) = recruitFrags(reg) ++ recruitUnpaired(reg)

  def recruitReads(reg: Region) = recruitLocalReads(reg) ++ recruitJumps(reg)
  
  def readsInInterval(bam: BamFile, name: String, start: Int, stop: Int) = {
    bam.reader.queryOverlapping(name, start, stop).toList
  }
  
  def mateMap(reads: List[SAMRecord]) = {
    val readMap = Map.empty[String, SAMRecord]
    val mm: MateMap = Map.empty
    for (r <- reads) {
      val readName = r.getReadName
      if (readMap contains readName) {
        val mate = readMap(readName)
        mm += r -> mate
        mm += mate -> r
      }
      readMap += readName -> r
    }
    
    mm
  }

  def mateMapForInterval(bam: BamFile, name: String, start: Int, stop: Int) = {
    val readMap = Map.empty[String, SAMRecord]
    val reads = readsInInterval(bam, name, start, stop)
    mateMap(reads)
  }
  
}
