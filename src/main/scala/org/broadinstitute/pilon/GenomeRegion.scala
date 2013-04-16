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

import collection.mutable.Map
import net.sf.picard.reference._
import Utils._
import java.io.PrintWriter

object GenomeRegion {
  def baseString(b: Array[Byte]) = b map { _.toChar.toUpper } mkString ("")
  type Fix = (Int, String, String)
}

class GenomeRegion(val contig: ReferenceSequence, start: Int, stop: Int)
  extends Region(contig.getName, start, stop) {
  require(stop <= contig.length, "GenomeRegion stop point must be within contig")
  val contigBases = contig.getBases
  var bases = refBases

  var pileUpRegion: PileUpRegion = null
  var minDepth = Pilon.minMinDepth

  var meanReadLength = 0

  // disposition flags
  var confirmed = new Array[Boolean](size)
  var ambiguous = new Array[Boolean](size)
  var changed = new Array[Boolean](size)
  var deleted = new Array[Boolean](size)
  //var lowCoverage = new Array[Boolean](size)
  //var lowConfidence = new Array[Boolean](size)

  // sumary stats
  val badCoverage = new Array[Int](size)
  val clips = new Array[Int](size)
  val copyNumber = new Array[Short](size)
  val coverage = new Array[Int](size)
  val insertSize = new Array[Int](size)
  val physCoverage = new Array[Int](size)
  val fragCoverage = new Array[Int](size)
  val weightedQual = new Array[Byte](size)
  val weightedMq = new Array[Byte](size)

  var changeMap = Map.empty[Int, (Symbol, PileUp)]

  def addChange(loc: Int, kind: Symbol, pu: PileUp) = {
    if (kind == 'amb) ambiguous(loc) = true
    else changed(loc) = true
    changeMap += (loc -> (kind, pu))
  }
  
  def changeList = changeMap.keys.toList.sorted

  lazy val physCoverageDist = new NormalDistribution(physCoverage, 2)
  lazy val coverageDist = new NormalDistribution(coverage, 2)
  lazy val fragCoverageDist = new NormalDistribution(fragCoverage, 2)
  lazy val badCoverageDist = new NormalDistribution(badCoverage, 2)
  lazy val insertSizeDist = new NormalDistribution(insertSize, 2)
  lazy val weightedMqDist = new NormalDistribution(weightedMq, 2)

  val gc = new Array[Byte](size)

  def kmerCopyNumber = {
    // really lame quick & dirty implementation, doesn't do rc, 
    // only for this contig, etc
    println("kcn start")
    val seq = (bases map {_.toChar toUpper}) mkString("")
	val k = Assembler.K
	val halfK = k / 2
    val kcn = new Array[Int](size)
    var kc = Map[String, Int]()
    val kmers = seq.sliding(k).toArray
    for (kmer <- kmers)
      kc(kmer) = kc.getOrElse(kmer, 0) + 1
    for (i <- 0 until kmers.size)
      kcn(i + halfK) = kc(kmers(i))
    println("kcn end")
    kcn
  }

  //val bams = Map.empty[Symbol, BamRegion]
  var bams = List[BamFile]()
  def bamsOfType(bamType: Symbol) = bams filter { _.bamType == bamType }

  def initializePileUps = {
    pileUpRegion = new PileUpRegion(name, start, stop)
  }
  
  def finalizePileUps = {
    pileUpRegion = null
  }
  
  
  def isGC(base: Char) = base == 'G' || base == 'C'

  def computeGc(window: Int = 100) = {
    val gcbuf = new Array[Byte](window)
    var count = 0
    var halfWindow = (window + 1) / 2
    // initialize circular buffer to 50%
    for (i <- 0 until window)
      if ((i & 1) == 0) {
        gcbuf(i) = 1
        count += 1
      } else
        gcbuf(i) = 0
    for (locus <- 0 until contig.length) {
      val center = if (locus >= halfWindow) locus - halfWindow else locus
      val bufIndex = locus % window
      val gcBase: Byte = {
        val base = contigBases(locus).toChar
        if (base == 'G' || base == 'C') 1 
        else if (base == 'A' || base == 'T') 0 
        else gcbuf(bufIndex)	//no-op for Ns, IUPAC, etc	
      }
      count += gcBase - gcbuf(bufIndex)
      gcbuf(bufIndex) = gcBase

      if (inRegion(center))
        gc(index(center)) = count.toByte
      if (inRegion(locus) && locus >= contig.length - halfWindow)
        gc(index(locus)) = count.toByte
    }
  }
  
  def smooth(input: Array[Int], window: Int) = {
    val result = new Array[Int](input.size)
    var accum = 0
    val half = window / 2
    for (i <- 0 until input.size) {
      accum += input(i)
      if (i > window) {
        accum -= input(i-window)
        val smoothed = (accum + half) / window
        result(i-half) = smoothed
      }
    }
    for (i <- 0 until window - half) result(i) = result(window-half)
    for (i <- result.size - half until result.size) result(i) = result(result.size-half-1)
    result
  }

  computeGc()

  def postProcess: Unit = {
    pileUpRegion.postProcess
    val meanCoverage = pileUpRegion.coverage
    val nReads = pileUpRegion.readCount

    minDepth = {
      if (Pilon.minDepth >= 1) Pilon.minDepth.toInt
      else (Pilon.minDepth * meanCoverage).round.toInt max Pilon.minMinDepth
    }

    println("Total Reads: " + nReads + ", Coverage: " + meanCoverage + ", minDepth: " + minDepth)
    //println("Non Jump coverage: " + fragCoverageDist.mean + ", median: " + fragCoverageDist.median)

    if (nReads == 0) {
      return
    }
    //meanReadLength = pileUpRegion.meanReadLength

    // Pass 1: pull out values from pileups & call base changes
    val fixamb = Pilon.fixList contains 'amb

    for (i <- 0 until size) {
      val pu = pileUpRegion(i)
      val n = pu.depth
      val bc = pu.baseCall
      val b = bc.base
      val q = bc.q
      val homo = bc.homo
      val m = pu.weightedMq
      //val (b, p, q, m) = bc
      val r = refBase(i + start)
      //val minDepth = meanCoverage * 25 / 100
      val loc = locus(i)
      coverage(i) = n.toInt
      badCoverage(i) = pu.badPair
      physCoverage(i) = pu.physCov
      insertSize(i) = pu.insertSize
      weightedQual(i) = pu.weightedQual.toByte
      weightedMq(i) = pu.weightedMq.toByte
      clips(i) = pu.clips.toShort

      if (n >= minDepth && r != 'N' && b != 'N' && !deleted(i)) {
        if (homo && b == r && bc.highConfidence && !bc.indel)
          confirmed(i) = true
        else if (bc.insertion) addChange(i, 'ins, pu)
        else if (bc.deletion) {
          addChange(i, 'del, pu)
          for (j <- 1 until pu.deletionCall.length) {
            deleted(i + j) = true
            pileUpRegion(i + j).deletions += pu.deletions
          }
        } else if (b != r) {
          // for ambiguous bases, fix them if --fix fixamb or if original base
          // not one of the top two alternatives
          if (homo) addChange(i, 'snp, pu)
          else if (fixamb || bc.altBase != r) addChange(i, 'amb, pu)
        }
      }
    }

    // Pass 2: computed values
    val baseCov = fragCoverageDist.mean
    println("Frag cov mean=" + baseCov + ", median=" + fragCoverageDist.median)
    val smoothCov = smooth(fragCoverage, 200)
    for (i <- 0 until size) {
      val n = smoothCov(i)
      val cn = if (baseCov > 0) (n / baseCov).round.toShort else 0.toShort
      copyNumber(i) = cn
    }

  }

  def processBam(bam: BamFile) = {
    print(bam.bamType.name + " " + bam.bamFile)
    // This is a real kludge...
    val covBefore = new Array[Int](size)
    if (bam.bamType != 'jumps)
      for (i <- 0 until size) 
        covBefore(i) = pileUpRegion.pileups(i).depth.toInt
    bam.process(this)
    if (bam.bamType != 'jumps)
      for (i <- 0 until size) 
        fragCoverage(i) += pileUpRegion.pileups(i).depth.toInt - covBefore(i)
    bams ::= bam
  }

  type FixList = List[GenomeRegion.Fix]
  var snpFixList: FixList = List()
  var smallFixList: FixList = List()
  var bigFixList: FixList = List()

  def identifyIssues = {
    if (Pilon.verbose) {
      println("# IdentifyIssues: " + this)
      identifyIssueRegions
    }
    var snps = 0
    var ins = 0
    var dels = 0
    var insBases = 0
    var delBases = 0
    var amb = 0
    for (i <- changeList) {
      val (kind, pu) = changeMap(i)
      val rBase = refBase(locus(i))
      val cBase = pu.baseCall.base
      kind match {
        case 'snp =>
          snpFixList ::= (locus(i), rBase.toString, cBase.toString)
          snps += 1
        case 'amb =>
          snpFixList ::= (locus(i), rBase.toString, cBase.toString)
          amb += 1
        case 'ins =>
          val insert = pu.insertCall
          smallFixList ::= (locus(i), "", insert)
          ins += 1
          insBases += insert.length
        case 'del =>
          val deletion = pu.deletionCall
          smallFixList ::= (locus(i), deletion, "")
          dels += 1
          delBases += deletion.length
      }
      if (Pilon.verbose) printChange(i)
    }
    
    // Report small changes
    print("Corrected ")
    if (Pilon.diploid) print((snps + amb) + " snps")
    else {
      print(snps + " snps; ")
      print(amb + " ambiguous bases")
    }
    print("; " + ins + " small insertions totaling " + insBases + " bases")
    println("; " + dels + " small deletions totaling " + delBases + " bases")
    
    // Report large collapsed regions (possible segmental duplication)
    val duplications = duplicationEvents
    if (duplications.size > 0) {
      for (d <- duplications) println("Large collapsed region: " + d + " size " + d.size)
    }
    
    // Apply SNP fixes prior to reassemblies...it helps by giving better anchor sequence!  
    // We can't change coords here, though, so no indels!
    fixIssues(snpFixList)

    // Try to fill gaps
    if ((Pilon.fixList contains 'gaps) && gaps.length > 0) {
      println("# Attempting to fill gaps")
      for (gap <- gaps) {
        val filler = new GapFiller(this)
        val (start, ref, patch) = filler.fillGap(gap)
        if (start > 0) {
          bigFixList ::= (start, ref, patch)
          printFix(gap, start, ref, patch, gap.size)
        }
      }
    }

    // Try to reassemble around possible contiguity breaks
    val breaks = possibleBreaks filter { !_.nearAny(gaps) }
    if ((Pilon.fixList contains 'local) && breaks.length > 0) {
      println("# Attempting to fix local continuity breaks")
      for (break <- breaks) {
        val filler = new GapFiller(this)
        val (start, ref, patch) = filler.fixBreak(break)
        if (start > 0 && (ref.length max patch.length) > 10) {
         	bigFixList ::= (start, ref, patch)
         	printFix(break, start, ref, patch, 0)
        } else if (Pilon.verbose || start == 0) {
          print ("# ")
          printFix(break, start, ref, patch, 0)
        }
      }
    }

    // Apply the bigger fixes
    fixIssues(smallFixList ++ bigFixList)
  }
  
  def printFix(reg: Region, loc: Int, ref: String, patch: String, gapSize: Int) = {
    def countNs(s: String) = s count {_ == 'N'}
    val nRef = countNs(ref)
    val nPatch = countNs(patch)
    val nonGapRef = ref.length - nRef
    val nonGapPatch = patch.length - nPatch
    val regStr = reg.start.toString + "(" + reg.size.toString + ")"
    val regType = if (gapSize > 0) "gap" else "break"
    print("fix " + regType + ": " + name + " " + regStr +
        " " + loc + " -" + nonGapRef + " +" + nonGapPatch)
    if (Pilon.verbose) {
    	print(" " + (if (ref.length > 0) ref else "."))
    	print(" " + (if (patch.length > 0) patch else "."))
    }
    if (loc == 0) print(" NoSolution")
    else if (gapSize == 0 && ref.length == 0 && patch.length == 0) print(" NoChange")
    else if (gapSize == 0 && nPatch == 0) print(" BreakFix")
    else if (nPatch > 0 && nRef == 0) print(" OpenedGap")
    else if (gapSize > 0 && nRef == gapSize && nPatch == 0) print(" ClosedGap")
    else if (gapSize > 0) print(" PartialFill")
    else print(" Unknown!") // Shouldn't happen...cases should be above!
    println()
  }
    

  def identifyIssueRegions = {
    println("# size=" + size + " medianCoverage=" + coverageDist.median + " meanCoverage=" + coverageDist.moments(0))
    //prettyPrintRegions("Unconfirmed", unConfirmedRegions)
    prettyPrintRegions("Gaps", gaps)
    prettyPrintRegions("LowCoverage", lowCoverageRegions)
    prettyPrintRegions("HighCN", highCopyNumberRegions)
    prettyPrintRegions("Break?", possibleBreaks)
    prettyPrintRegions("Clip?", clippingRegions)
    prettyPrintRegions("Insertion?", possibleInsertions)
    prettyPrintRegions("Deletion?", possibleDeletions)
    prettyPrintRegions("CollapsedRepeat?", possibleCollapsedRepeats)
    prettyPrintRegions("Ambiguous", ambiguousRegions)

  }

  def printChange(i: Int, endLine: String = "\n") = {
    val (kind, pu) = changeMap(i)
    val rBase = refBase(locus(i))
    val cBase = pu.baseCall.base
    kind match {
      case 'snp =>
        print(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase)
        if (Pilon.debug) print(" " + pu)
        print(endLine)
      case 'ins =>
        print(name + " " + locus(i) + " " + kind.name + " " + "." + " " + pu.insertCall)
        if (Pilon.debug) print(" " + pu)
        print(endLine)
      case 'del =>
        print(name + " " + locus(i) + " " + kind.name + " " + pu.deletionCall + " " + ".")
        if (Pilon.debug) print(" " + pu)
        print(endLine)
      case 'amb =>
        if (Pilon.verbose && rBase != cBase) {
          print(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase)
          if (Pilon.debug) print(" " + pu)
          print(endLine)
        }
    }
  }

  def prettyPrintRegions(header: String, regions: List[Region]) = {
    if (regions != Nil) {
      val totalSize: Int = (regions map { _.size }).sum
      print("# " + header + " n=" + regions.size + " bases=" + totalSize)
      if (Pilon.debug) {
    	regions foreach { r => print("  " + r.start + "-" + r.stop) }
      }
      println
    }
  }

  def fixFixList(inList: FixList): FixList = {
    """Sort and remove overlaps, keeping larger fix"""
    
    var fixes = inList.sortWith({ (x, y) => x._1 < y._1 })
    var outList: FixList = Nil

    while (fixes != Nil) {
      //println(fixes.head)
      fixes match {
        case fix1 :: fix2 :: tail => {
          val region1 = new Region(name, fix1._1, fix1._1 + ((fix1._2.length - 1) max 0))
          val region2 = new Region(name, fix2._1, fix2._1 + ((fix2._2.length - 1) max 0))
          if (region1.overlaps(region2)) {
            // If we have overlapping changes, keep the one with the most impact
            val fix1len = fix1._2.length + fix1._3.length
            val fix2len = fix2._2.length + fix2._3.length
            if (fix1len > fix2len)
              fixes = fix1 :: tail
            else
              fixes = fix2 :: tail
          } else {
            fixes = fix2 :: tail
            outList ::= fix1
          }
        }
        case fix1 :: Nil =>
          fixes = Nil
          outList ::= fix1
        case Nil =>
          outList
      }
    }
    outList.reverse
  }

  def fixIssues(fixList: FixList) = {
    var newBases = bases.clone
    for (fix <- fixFixList(fixList).reverse) {
      if (Pilon.debug) println("Fix " + fix)
      val (locus, was, patch) = fix
      val start = index(locus)
      if (was.length == patch.length) {
        for (i <- 0 until was.length) {
          val iNew = start + i
          val ref = newBases(iNew).toChar.toUpper
          if (ref != was(i)) 
            println("Fix mismatch: loc=" + (locus+i) + " ref=" + ref + " was=" + was(i))
          newBases(start + i) = patch(i).toByte        
        }
      } else {
        val origLength = newBases.length
        val before = newBases.slice(0, start)
        val after = newBases.slice(start + was.length, newBases.length)
        val ref = newBases.slice(start, start + was.length).map(_.toChar).mkString("").toUpperCase
        if (ref != was) println("Fix mismatch: loc=" + locus + " ref=" + ref + " was=" + was)
        newBases = before ++ (patch map { _.toByte }) ++ after
        assert(newBases.length == origLength + patch.length - was.length, "Fix patch length mismatch: " + fix)
        if (Pilon.debug) println("Fixing=" + was.length + " " + patch.length + " " + newBases.length)
      }
    }
    bases = newBases
  }

  def writeVcf(vcf: Vcf) = {
    var bigFixes = fixFixList(bigFixList)
    for (i <- 0 until size) {
      val loc = locus(i)
      if (bigFixes.length > 0 && bigFixes.head._1 == loc) {
        vcf.writeFixRecord(this, bigFixes.head)
        for (j <- 0 until bigFixes.head._2.length) deleted(i+j) = true
        bigFixes = bigFixes.tail
      }
      vcf.writeRecord(this, i, deleted(i))
    }
  }

  def writeChanges(changes: PrintWriter, newName: String = name) {
    val fixes = fixFixList(snpFixList ++ smallFixList ++ bigFixList)
    var delta = 0
    for (fix <- fixes) {
      val (loc, from, to) = fix
      val newLoc = loc + delta
      val oldRegion = new Region(name, loc, loc + from.length - 1)
      val newRegion = new Region(newName, newLoc, newLoc + to.length - 1)
      changes.println(oldRegion.regionString + " " + newRegion.regionString + " " +
        (if (from.isEmpty) "." else from) + " " + (if (to.isEmpty) "." else to))
      delta += to.length - from.length
    }
  }

  def maxBreak(r: Region) = {
    var maxDeltaLoc = r.start
    var maxDelta = deltaFraction(maxDeltaLoc)
    for (i <- maxDeltaLoc + 1 to r.stop) {
      val delta = deltaFraction(i)
      if (delta > maxDelta) {
        maxDelta = delta
        maxDeltaLoc = i
      }
    }
    maxDeltaLoc
  }

  def delta(i: Int, values: Array[Int], radius: Int = 100) = {
    val left = values(0 max (i - radius))
    val right = values((size - 1) min (i + radius))
    val center = values(i)
    (left - right).abs
    //(left - center).abs + (right - center).abs
  }

  def deltaCoverage(i: Int, radius: Int = 100) = delta(i, fragCoverage, radius)

  def deltaPhysicalCoverage(i: Int, radius: Int = 1000) = delta(i, physCoverage, radius)

  def dip(i: Int, values: Array[Int], radius: Int = 100) = {
    val left = values(0 max (i - radius))
    val right = values((size - 1) min (i + radius))
    val center = values(i)
    (left - center) + (right - center)
  }

  def nearEdge(r: Region, radius: Int = 100) = r.start - start < radius || stop - r.stop < radius
  def deltaFraction(i: Int) = deltaCoverage(i) / fragCoverageDist.mean
  def dipCoverage(i: Int, radius: Int = 100) = dip(i, fragCoverage, radius)
  def dipFraction(i: Int) = dipCoverage(i) / fragCoverageDist.mean

  def ambiguousRegions = summaryRegions({ i: Int => ambiguous(i) }) filter { r => r.start != r.stop }
  def changeRegions = summaryRegions({ i: Int => changed(i) || deleted(i) }, 1)
  def gaps = summaryRegions({ i: Int => refBase(locus(i)) == 'N' }) filter { _.size >= 10}
  def highCopyNumberRegions = summaryRegions({ i: Int => copyNumber(i) > 1 })
  def unConfirmedRegions = summaryRegions({ i: Int => !confirmed(i) })

  def lowCoverage(i: Int) = coverage(i) < Pilon.minMinDepth
  def lowCoverageRegions = summaryRegions(lowCoverage)
  def highClipping(i: Int) = coverage(i) >= Pilon.minMinDepth && pct(clips(i), coverage(i)) >= 33
  def clippingRegions = summaryRegions(highClipping)
  def pctBad(i: Int) = {
    val good = coverage(i)
    val bad = badCoverage(i)
    pct(bad, good + bad)
  }

  def breakp(i: Int) = lowCoverage(i) || highClipping(i) || (pctBad(i) >= 50) || (dipFraction(i) >= 1.5)
  def possibleBreaks = summaryRegions(breakp)

  def insertionp(i: Int) = breakp(i) && insertSizeDist.toSigma(insertSize(i)) <= -3.0
  def possibleInsertions = summaryRegions(insertionp)

  def deletionp(i: Int) = breakp(i) && insertSizeDist.toSigma(insertSize(i)) >= 3.0
  def possibleDeletions = summaryRegions(deletionp)

  def possibleCollapsedRepeats = {
    highCopyNumberRegions filter { r =>
      //println(r + " " + deltaFraction(index(r.start)) + " " + deltaFraction(index(r.stop)))
      deltaFraction(index(r.start)) >= 0.5 && deltaFraction(index(r.stop)) >= 0.5
    }
  }
  
  def duplicationEvents = {
    //val smoothCoverage = smooth(fragCoverage, 1000)
    //val median = fragCoverageDist.median
    //val regions = summaryRegions({ i: Int => smoothCoverage(i)/median > 1.5 }, 1000) 
    val regions = summaryRegions({ i: Int => copyNumber(i) > 1 }, 2000) 
    regions filter {_.size > 10000}
    
  }

  def summaryRegions(positionTest: (Int) => Boolean, slop: Int = 100) = {
    var regions = List[Region]()
    var first = -1
    var last = -1

    for (i <- 0 until size) {
      if (positionTest(i)) {
        last = i
        if (first < 0)
          first = i
      } else {
        if (last >= 0 && i > last + slop) {
          regions ::= new Region(name, locus(first), locus(last))
          first = -1
          last = -1
        }
      }
    }
    if (last >= 0) regions ::= new Region(name, locus(first), locus(last))
    regions.reverse filter { !nearEdge(_) }
  }

  
  def baseAt(locus: Int) = {
    require(inRegion(locus))
    bases(index(locus)).toChar.toUpper
  }
  
  def subString(locus: Int, n: Int) = {
    require(inRegion(locus) && inRegion(locus + n - 1))
    GenomeRegion.baseString(bases.slice(index(locus), index(locus+n)))
  } 

  def refBase(locus: Int) = {
    require(inRegion(locus), "can't fetch base outside region")
    //bases(index(locus)).toChar.toUpper
    contigBases(locus - 1).toChar.toUpper
  }
  
  def refBases = {
    contigBases.slice(start - 1, stop)
  }

  def refSlice(a: Int, b: Int) = contigBases.slice(index(a), index(b))
  def refBases(a: Int, b: Int) = refSlice(a, b).toArray

  def refSubString(locus: Int, n: Int) = {
    require(inRegion(locus) && inRegion(locus + n - 1))
    GenomeRegion.baseString(contigBases.slice(locus - 1, locus + n - 1))
  }

  //def apply(locus: Int) = bases(index(locus))
}