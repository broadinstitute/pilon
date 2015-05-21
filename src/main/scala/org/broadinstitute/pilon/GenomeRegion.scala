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

import collection.mutable.Map
import htsjdk.samtools.reference._
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
  val originalBases = refBases
  var bases = originalBases

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

  var logString = ""

  // Hack alert...if we're running multithreaded, buffer our output for atomic
  // write when we're done. Otherwise, we would rather see output as it happens.
  def log(s: String) = {
    if (Pilon.threads > 1)
      logString += s
    else
      print(s)
  }
  
  def logln(s: String = "") = {
    log(s + "\n")
  }
  
  def printLog() = if (Pilon.threads > 1) print(logString)
  
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
  lazy val pctBadOverall = {
    // Compute the overall percentage of bad coverage to total coverage
    // in this GenomeRegion. Used to determine a threshold for evidence
    // of contiguity breaks.

    // do calculations as longs or we'll overflow
    var totalBad = 0L
    var totalGood = 0L
    for (i <- 0 until size) {
      totalBad += badCoverage(i)
      totalGood += coverage(i)
    }
    if (Pilon.verbose)
      logln("pctBadOverall: " + totalBad + " " + totalGood + " " + pct(totalBad, totalGood + totalBad))
    pct(totalBad, totalGood + totalBad)
  }

  val gc = new Array[Byte](size)

  def kmerCopyNumber = {
    // really lame quick & dirty implementation, doesn't do rc, 
    // only for this contig, etc
    logln("kcn start")
    val seq = (bases.map({_.toChar.toUpper})).mkString("")
	val k = Assembler.K
	val halfK = k / 2
    val kcn = new Array[Int](size)
    var kc = Map[String, Int]()
    val kmers = seq.sliding(k).toArray
    for (kmer <- kmers)
      kc(kmer) = kc.getOrElse(kmer, 0) + 1
    for (i <- 0 until kmers.size)
      kcn(i + halfK) = kc(kmers(i))
    logln("kcn end")
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
    val inputSize = input.size
    val result = new Array[Int](inputSize)
    val half = window / 2
    var accum = 0
    for (i <- 0 until inputSize) {
      accum += input(i)
      if (i > window) {
        accum -= input(i-window)
        val smoothed = (accum + half) / window
        result(i-half) = smoothed
      }
    }
    if (inputSize > window) {
      for (i <- 0 until window - half) result(i) = result(window-half)
      for (i <- inputSize - half until inputSize) result(i) = result(inputSize-half-1)
    } else {
      for (i <- 0 until inputSize) result(i) = accum / inputSize
    }
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

    logln("Total Reads: " + nReads + ", Coverage: " + meanCoverage + ", minDepth: " + minDepth)
    //logln("Non Jump coverage: " + fragCoverageDist.mean + ", median: " + fragCoverageDist.median)

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
      val homo = bc.homo
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

      if (n >= minDepth && r != 'N' && !deleted(i) && bc.called) {
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
    if (Pilon.verbose)
      logln("Frag coverage mean=" + baseCov + ", median=" + fragCoverageDist.median)
    val smoothCov = smooth(fragCoverage, 200)
    for (i <- 0 until size) {
      val n = smoothCov(i)
      val cn = if (baseCov > 0) (n / baseCov).round.toShort else 0.toShort
      copyNumber(i) = cn
    }

  }

  def processBam(bam: BamFile) = {
    log(bam.bamType.name + " " + bam.bamFile + ": ")
    // This is a real kludge...
    val covBefore = new Array[Int](size)
    if (bam.bamType != 'jumps)
      for (i <- 0 until size) 
        covBefore(i) = pileUpRegion.pileups(i).depth.toInt
    val coverage = bam.process(this)
    logln("coverage " + coverage.toString)
    if (bam.bamType != 'jumps)
      for (i <- 0 until size) 
        fragCoverage(i) += pileUpRegion.pileups(i).depth.toInt - covBefore(i)
    bams ::= bam
  }

  type FixList = List[GenomeRegion.Fix]
  var snpFixList: FixList = List()
  var smallFixList: FixList = List()
  var bigFixList: FixList = List()

  def identifyAndFixIssues = {
    if (Pilon.verbose) {
      logln("# IdentifyIssues: " + this)
      identifyIssueRegions
    }
    val fix = Pilon.fixList contains 'bases
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
          if (fix) snpFixList ::= (locus(i), rBase.toString, cBase.toString)
          snps += 1
        case 'amb =>
          if (fix) snpFixList ::= (locus(i), rBase.toString, cBase.toString)
          amb += 1
        case 'ins =>
          val insert = pu.insertCall
          if (fix) smallFixList ::= (locus(i), "", insert)
          ins += 1
          insBases += insert.length
        case 'del =>
          val deletion = pu.deletionCall
          if (fix) smallFixList ::= (locus(i), deletion, "")
          dels += 1
          delBases += deletion.length
      }
      if (Pilon.verbose) logChange(i)
    }
    
    // Report some stats
    val nConfirmed = confirmed count {x => x}
    val nonN = originalBases count {x => x != 'N'}
    logln("Confirmed " + nConfirmed + " of " + nonN + " bases (" +
      (nConfirmed * 100.0 / nonN).formatted("%.2f") + "%)")
    if (Pilon.fixList contains 'bases) log("Corrected ") else log("Found ")
    if (Pilon.diploid) log((snps + amb) + " snps")
    else {
      log(snps + " snps; ")
      log(amb + " ambiguous bases")
    }
    log("; " + ins + " small insertions totaling " + insBases + " bases")
    logln("; " + dels + " small deletions totaling " + delBases + " bases")
    
    // Report large collapsed regions (possible segmental duplication)
    val duplications = duplicationEvents
    if (duplications.size > 0) {
      for (d <- duplications) logln("Large collapsed region: " + d + " size " + d.size)
    }
    
    // Apply SNP fixes prior to reassemblies...it helps by giving better anchor sequence!  
    // We can't change coords here, though, so no indels!
    fixIssues(snpFixList)

    // Try to fill gaps
    if ((Pilon.fixList contains 'gaps) && gaps.length > 0) {
      logln("# Attempting to fill gaps")
      for (gap <- gaps) {
        val filler = new GapFiller(this)
        val (start, ref, patch) = filler.fillGap(gap)
        if (start > 0) {
          bigFixList ::= (start, ref, patch)
          logFix(gap, start, ref, patch, gap.size, filler.tandemRepeat)
        }
      }
    }

    // Try to reassemble around possible contiguity breaks, but stay away from gaps
    val breaks = possibleBreaks filter { !_.nearAny(gaps, 300) }
    if ((Pilon.fixList contains 'local) && breaks.length > 0) {
      logln("# Attempting to fix local continuity breaks")
      for (break <- breaks) {
        val filler = new GapFiller(this)
        val (start, ref, patch) = filler.fixBreak(break)
        if (start > 0 && (ref.length max patch.length) > 10) {
         	bigFixList ::= (start, ref, patch)
         	logFix(break, start, ref, patch, 0, filler.tandemRepeat)
        } else if (Pilon.verbose || start == 0) {
          log ("# ")
          logFix(break, start, ref, patch, 0, filler.tandemRepeat)
        }
      }
    }

    // Apply the bigger fixes
    fixIssues(smallFixList ++ bigFixList)
  }
  
  def logFix(reg: Region, loc: Int, ref: String, patch: String, gapSize: Int, tandemRepeat: String = "") = {
    def countNs(s: String) = s count {_ == 'N'}
    val nRef = countNs(ref)
    val nPatch = countNs(patch)
    val nonGapRef = ref.length - nRef
    val nonGapPatch = patch.length - nPatch
    val regStr = reg.start.toString + "(" + reg.size.toString + ")"
    val regType = if (gapSize > 0) "gap" else "break"
    //log("fix " + regType + ": " + name + " " + regStr +
    //      " " + loc + " -" + nonGapRef + " +" + nonGapPatch)
    log("fix " + regType + ": " + reg +
          " " + loc + " -" + nonGapRef + " +" + nonGapPatch)
    if (Pilon.verbose) {
    	log(" " + (if (ref.length > 0) ref else "."))
    	log(" " + (if (patch.length > 0) patch else "."))
    }
    if (loc == 0) log(" NoSolution")
    else if (gapSize == 0 && ref.length == 0 && patch.length == 0) log(" NoChange")
    else if (gapSize == 0 && nPatch == 0) log(" BreakFix")
    else if (nPatch > 0 && nRef == 0) log(" OpenedGap")
    else if (gapSize > 0 && nRef == gapSize && nPatch == 0) log(" ClosedGap")
    else if (gapSize > 0) log(" PartialFill")
    else log(" Unknown!") // Shouldn't happen...cases should be above!
    if (tandemRepeat != "") {
      log (" TandemRepeat " + tandemRepeat.length)
      if (Pilon.verbose) log(" " + tandemRepeat)
    }
    logln()
  }
    

  def identifyIssueRegions = {
    logln("# size=" + size + " medianCoverage=" + coverageDist.median + " meanCoverage=" + coverageDist.moments(0))
    //prettylogRegions("Unconfirmed", unConfirmedRegions)
    prettylogRegions("Gaps", gaps)
    prettylogRegions("LowCoverage", lowCoverageRegions)
    prettylogRegions("HighCN", highCopyNumberRegions)
    prettylogRegions("Break?", possibleBreaks)
    prettylogRegions("Clip?", clippingRegions)
    prettylogRegions("Insertion?", possibleInsertions)
    prettylogRegions("Deletion?", possibleDeletions)
    prettylogRegions("CollapsedRepeat?", possibleCollapsedRepeats)
    prettylogRegions("Ambiguous", ambiguousRegions)

  }

  def logChange(i: Int, endLine: String = "\n") = {
    val (kind, pu) = changeMap(i)
    val rBase = refBase(locus(i))
    val cBase = pu.baseCall.base
    kind match {
      case 'snp =>
        log(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase)
        if (Pilon.debug) log(" " + pu)
        log(endLine)
      case 'ins =>
        log(name + " " + locus(i) + " " + kind.name + " " + "." + " " + pu.insertCall)
        if (Pilon.debug) log(" " + pu)
        log(endLine)
      case 'del =>
        log(name + " " + locus(i) + " " + kind.name + " " + pu.deletionCall + " " + ".")
        if (Pilon.debug) log(" " + pu)
        log(endLine)
      case 'amb =>
        if (Pilon.verbose && rBase != cBase) {
          log(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase)
          if (Pilon.debug) log(" " + pu)
          log(endLine)
        }
    }
  }

  def prettylogRegions(header: String, regions: List[Region]) = {
    if (regions != Nil) {
      val totalSize: Int = (regions map { _.size }).sum
      log("# " + header + " n=" + regions.size + " bases=" + totalSize)
      if (Pilon.debug) {
    	regions foreach { r => log("  " + r.start + "-" + r.stop) }
      }
      logln()
    }
  }

  def fixFixList(inList: FixList): FixList = {
    // Sort and remove overlaps, keeping larger fix
    
    var fixes = inList.sortWith({ (x, y) => x._1 < y._1 })
    var outList: FixList = Nil

    while (fixes != Nil) {
      //logln(fixes.head)
      fixes match {
        case fix1 :: fix2 :: tail => {
          val region1 = new Region(name, fix1._1, fix1._1 + ((fix1._2.length - 1) max 0))
          val region2 = new Region(name, fix2._1, fix2._1 + ((fix2._2.length - 1) max 0))
          if (region1.overlaps(region2)) {
            // If we have overlapping changes, keep the one with the most impact
            val fix1len = fix1._2.length + fix1._3.length
            val fix2len = fix2._2.length + fix2._3.length
            // Prefer fix1 (first in coord order) if same length, as pileup-based indels
            // get shifted to leftmost coord, and we would rather have those than reassemblies,
            // all else being equal.
            if (fix1len >= fix2len)
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
      if (Pilon.debug) logln("Fix " + fix)
      val (locus, was, patch) = fix
      val start = index(locus)
      if (was.length == patch.length) {
        for (i <- 0 until was.length) {
          val iNew = start + i
          val ref = originalBases(iNew).toChar.toUpper
          if (ref != was(i)) 
            logln("Fix mismatch: loc=" + (locus+i) + " ref=" + ref + " was=" + was(i))
          newBases(start + i) = patch(i).toByte        
        }
      } else {
        val origLength = newBases.length
        val before = newBases.slice(0, start)
        val after = newBases.slice(start + was.length, newBases.length)
        val ref = originalBases.slice(start, start + was.length).map(_.toChar).mkString("").toUpperCase
        if (ref != was) logln("Fix mismatch: loc=" + locus + " ref=" + ref + " was=" + was)
        newBases = before ++ (patch map { _.toByte }) ++ after
        assert(newBases.length == origLength + patch.length - was.length, "Fix patch length mismatch: " + fix)
        if (Pilon.debug) logln("Fixing=" + was.length + " " + patch.length + " " + newBases.length)
      }
    }
    bases = newBases
  }

  def writeVcf(vcf: Vcf) = {
    var fixes = fixFixList(snpFixList ++ smallFixList ++ bigFixList)
    var dupes = duplicationEvents
    for (i <- 0 until size) {
      val loc = locus(i)
      if (dupes != Nil && dupes.head.start == loc) {
        vcf.writeDup(this, dupes.head)
        dupes = dupes.tail
      }
      if (fixes.length > 0 && fixes.head._1 == loc) {
        val fix = fixes.head
        // We write a special record if this was a big fix (local reassembly)
        if ((bigFixList contains fix) && !(smallFixList contains fix)) {
          vcf.writeFixRecord(this, fix)
          for (j <- 0 until fix._2.length) deleted(i+j) = true
        }
        fixes = fixes.tail
      }
      vcf.writeRecord(this, i, deleted(i))
    }
  }

  def writeChanges(changes: PrintWriter, newName: String = name, offset: Int = 0) {
    val fixes = fixFixList(snpFixList ++ smallFixList ++ bigFixList)
    var delta = 0
    for (fix <- fixes) {
      val (loc, from, to) = fix
      val newLoc = loc + delta
      val oldRegion = new Region(name, loc, loc + from.length - 1)
      val newRegion = new Region(newName, newLoc + offset, newLoc + offset + to.length - 1)
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

  def tooBad(i: Int) = {
    // Heuristic to see whether there is an abnormal percentage of "bad" coverage at this locus.
    // Taking the overall badness + 20% for now (previously was just pctBad >= 50).
    val good = coverage(i)
    val bad = badCoverage(i)
    val p = pct(bad, good + bad)
    //p >= 50 // (old way)
    p > pctBadOverall + 20
  }

  def breakp(i: Int) = lowCoverage(i) || highClipping(i) || tooBad(i) || (dipFraction(i) >= 1.5)
  def possibleBreaks = summaryRegions(breakp, 200)

  def insertionp(i: Int) = breakp(i) && insertSizeDist.toSigma(insertSize(i)) <= -3.0
  def possibleInsertions = summaryRegions(insertionp)

  def deletionp(i: Int) = breakp(i) && insertSizeDist.toSigma(insertSize(i)) >= 3.0
  def possibleDeletions = summaryRegions(deletionp)

  def possibleCollapsedRepeats = {
    highCopyNumberRegions filter { r =>
      //logln(r + " " + deltaFraction(index(r.start)) + " " + deltaFraction(index(r.stop)))
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

  
  def baseAt(locus: Int, whichBases: Array[Byte] = bases) = {
    require(inRegion(locus))
    whichBases(index(locus)).toChar.toUpper
  }

  def originalBaseAt(locus: Int) = baseAt(locus, originalBases)
  
  def subString(locus: Int, n: Int) = {
    require(inRegion(locus) && inRegion(locus + n - 1))
    GenomeRegion.baseString(bases.slice(index(locus), index(locus+n)))
  }

  def refSubString(locus: Int, n: Int) = {
    require(inRegion(locus) && inRegion(locus + n - 1))
    GenomeRegion.baseString(contigBases.slice(locus - 1, locus + n - 1))
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


  //def apply(locus: Int) = bases(index(locus))
}
