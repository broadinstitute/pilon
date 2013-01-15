package org.broadinstitute.pilon

import collection.mutable.Map
import net.sf.picard.reference._

object GenomeRegion {
  def baseString(b: Array[Byte]) = b map { _.toChar } mkString ("")
  type Fix = (Int, String, String)
}

class GenomeRegion(val contig: ReferenceSequence, start: Int, stop: Int)
  extends Region(contig.getName, start, stop) {
  require(stop <= contig.length, "GenomeRegion stop point must be within contig")
  val bases = contig.getBases
  var fixedBases = refBases.clone

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
  val badCoverage = new Array[Short](size)
  val clips = new Array[Short](size)
  val copyNumber = new Array[Short](size)
  val coverage = new Array[Short](size)
  val insertSize = new Array[Short](size)
  val physCoverage = new Array[Short](size)
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
  lazy val badCoverageDist = new NormalDistribution(badCoverage, 2)
  lazy val insertSizeDist = new NormalDistribution(insertSize, 2)
  lazy val weightedMqDist = new NormalDistribution(weightedMq, 2)

  val gc = new Array[Byte](size)

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
        val base = bases(locus).toChar
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

    if (nReads == 0) {
      return
    }
    //meanReadLength = pileUpRegion.meanReadLength

    // Pass 1: pull out values from pileups & call base changes

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
      coverage(i) = n.toShort
      badCoverage(i) = pu.badPair.toShort
      physCoverage(i) = pu.physCov.toShort
      insertSize(i) = pu.insertSize.toShort
      weightedQual(i) = pu.weightedQual.toByte
      weightedMq(i) = pu.weightedMq.toByte
      clips(i) = pu.clips.toShort

      if (n >= minDepth /*&& r != 'N'*/ ) {
        if (homo && b == r && bc.highConfidence && !bc.indel)
          confirmed(i) = true
        else if (bc.insertion) addChange(i, 'ins, pu)
        else if (bc.deletion) {
          addChange(i, 'del, pu)
          for (j <- 1 until pu.deletionCall.length) deleted(i + j) = true
        } else if (b != r) {
          if (homo) addChange(i, 'snp, pu)
          else addChange(i, 'amb, pu)
        }
      }
    }

    // Pass 2: computed values
    val baseCov = coverageDist.median
    for (i <- 0 until size) {
      val n = coverage(i)
      val cn = if (baseCov > 0) (n / baseCov).round.toShort else 0.toShort
      copyNumber(i) = cn
    }

  }

  def processBam(bam: BamFile) = {
    print(bam.bamType.name + " " + bam.bamFile)
    //val br = new BamRegion(bam, this)
    bam.process(this)
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
          //newBases(i) = cBase.toByte
          snpFixList ::= (locus(i), rBase.toString, cBase.toString)
          snps += 1
        case 'ins =>
          val insert = pu.insertCall
          smallFixList ::= (locus(i), "", insert)
          ins += 1
          insBases += insert.length
        case 'del =>
          val deletion = pu.deletionCall
          smallFixList ::= (locus(i), "", deletion)
          dels += 1
          delBases += deletion.length
        case 'amb =>
          smallFixList ::= (locus(i), rBase.toString, cBase.toString)
          amb += 1
        case _ =>
      }
      if (Pilon.verbose) printChange(i)
    }
    if (Pilon.diploid) print((snps + amb) + " snps")
    else {
      print(snps + " snps; ")
      print(amb + " ambiguous bases")
    }
    print("; " + ins + " small insertions totaling " + insBases + " bases")
    println("; " + dels + " small deletions totaling " + delBases + " bases")
    
    // fix SNPs prior to reassemblies...it helps!  We can't change coords here, though!
    fixIssues(snpFixList)

    if ((Pilon.fixList contains 'gaps) && gaps.length > 0) {
      println("# Attempting to fill gaps")
      for (gap <- gaps) {
        val filler = new GapFiller(this)
        val (start, ref, patch) = filler.fillGap(gap)
        if (start > 0) {
          val old = if (ref != "") ref else "."
          val now = if (patch != "") patch else "."
          //println(name + " " + start + " fix " + ref.size + " " + old + " " + patch.size + " " + now)
          printFix(start, ref, patch)
          bigFixList ::= (start, ref, patch)
        }
      }
    }

    val breaks = possibleBreaks //clippingRegions
    val filteredBreaks = breaks filter { !_.nearAny(gaps) }
    if ((Pilon.fixList contains 'local) && breaks.length > 0) {
      println("# Attempting to fix local continuity breaks")
      //val breaks = possibleBreaks

      for (break <- filteredBreaks) {
        val filler = new GapFiller(this)
        val (start, ref, patch) = filler.fixBreak(break)
        if (start > 0) {
          val old = if (ref != "") ref else "."
          val now = if (patch != "") patch else "."
          if (ref.length + patch.length >= 20) {
            printFix(start, ref, patch)
        	bigFixList ::= (start, ref, patch)
          }
        }
      }
    }
    fixIssues(smallFixList ++ bigFixList)
  }
  
  def printFix(loc: Int, ref: String, patch: String) = {
    def countNs(s: String) = s count {_ == 'N'}
    val nRef = countNs(ref)
    val nPatch = countNs(patch)
    val nonGapRef = ref.length - nRef
    val nonGapPatch = patch.length - nPatch
    
    print("fix: " + name + " " + loc + " -" + nonGapRef + " +" + nonGapPatch)
    if (Pilon.verbose) {
    	print(" " + (if (ref.length > 0) ref else "."))
    	print(" " + (if (patch.length > 0) patch else "."))
    }
    if (nRef > 0 && nPatch == 0) print(" ClosedGap")
    else if (nPatch > 0 && nRef == 0) print (" OpenedGap")
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
        if (Pilon.debug) print(" " + pu + endLine)
        print(endLine)
      case 'ins =>
        print(name + " " + locus(i) + " " + kind.name + " " + "." + " " + pu.insertCall)
        if (Pilon.debug) print(" " + pu + endLine)
        print(endLine)
      case 'del =>
        print(name + " " + locus(i) + " " + kind.name + " " + pu.deletionCall + " " + ".")
        if (Pilon.debug) print(" " + pu + endLine)
        print(endLine)
      case 'amb =>
        if (Pilon.verbose && rBase != cBase) {
          print(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase)
          if (Pilon.debug) print(" " + pu + endLine)
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


  def fixFixListX(inList: FixList, outList: FixList = Nil): FixList = {
    """Remove overlaps, keeping larger fix"""
    inList match {
      case fix1 :: fix2 :: tail => {
        val region1 = new Region(name, fix1._1, fix1._1 + fix1._2.length - 1)
        val region2 = new Region(name, fix2._1, fix2._1 + fix2._2.length - 1)
        if (region1.overlaps(region2)) {
          //println("fixFixList overlap: " + (region1, region2))
          if (region1.size > region2.size)
            fixFixListX(fix1 :: tail, outList)
          else
            fixFixListX(fix2 :: tail, outList)
        } else
          fixFixListX(fix2 :: tail, fix1 :: outList)
      }
      case fix1 :: Nil =>
        fix1 :: outList
      case Nil =>
        outList
    }
  }

  def fixFixList(inList: FixList): FixList = {
    """Remove overlaps, keeping larger fix"""
    var fixes = inList
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
    outList
  }

  def fixIssues(fixList: FixList) = {
    var newBases = fixedBases.clone
    val sortedFixes = fixList.sortWith({ (x, y) => x._1 < y._1 })
    val fixedFixes = fixFixList(sortedFixes)
    for (fix <- fixedFixes) {
      if (Pilon.debug) println("Fix " + fix)
      val (locus, was, patch) = fix
      val start = index(locus)
      if (was.length == patch.length) {
        for (i <- 0 until was.length)
          newBases(start + i) = patch(i).toByte        
      } else {
    	val before = newBases.slice(0, start)
        val after = newBases.slice(start + was.length, newBases.length)
        val ref = newBases.slice(start, start + was.length).map(_.toChar).mkString("").toUpperCase
        if (ref != was) println("Fix mismatch: loc=" + locus + " ref=" + ref + " was=" + was)
        newBases = before ++ (patch map { _.toByte }) ++ after
        //if (Pilon.debug) println("Fixing=" + was.length + " " + patch.length + " " + newBases.length)
      }
    }
    newBases
  }

  def writeVcf(vcf: Vcf) = {
	var bigFixes = bigFixList.reverse
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

  def delta(i: Int, values: Array[Short], radius: Int = 100) = {
    val left = values(0 max (i - radius))
    val right = values((size - 1) min (i + radius))
    val center = values(i)
    (left - right).abs
    //(left - center).abs + (right - center).abs
  }

  def deltaCoverage(i: Int, radius: Int = 100) = delta(i, coverage, radius)

  def deltaPhysicalCoverage(i: Int, radius: Int = 1000) = delta(i, physCoverage, radius)

  def dip(i: Int, values: Array[Short], radius: Int = 100) = {
    val left = values(0 max (i - radius))
    val right = values((size - 1) min (i + radius))
    val center = values(i)
    (left - center) + (right - center)
  }

  def nearEdge(r: Region, radius: Int = 100) = r.start - start < radius || stop - r.stop < radius
  def deltaFraction(i: Int) = deltaCoverage(i) / coverageDist.mean
  def dipCoverage(i: Int, radius: Int = 100) = dip(i, coverage, radius)
  def dipFraction(i: Int) = dipCoverage(i) / coverageDist.mean

  def ambiguousRegions = summaryRegions({ i: Int => ambiguous(i) }) filter { r => r.start != r.stop }
  def changeRegions = summaryRegions({ i: Int => changed(i) })
  def gaps = summaryRegions({ i: Int => refBase(locus(i)) == 'N' }) filter { _.size >= 10}
  def highCopyNumberRegions = summaryRegions({ i: Int => copyNumber(i) > 1 })
  def unConfirmedRegions = summaryRegions({ i: Int => !confirmed(i) })

  def lowCoverage(i: Int) = coverage(i) < minDepth
  def lowCoverageRegions = summaryRegions(lowCoverage)
  def highClipping(i: Int) = clips(i) >= coverage(i)
  def clippingRegions = summaryRegions(highClipping)

  def breakp(i: Int) = highClipping(i) || lowCoverage(i) // dipFraction(i) >= 1.0
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

  def refBase(locus: Int) = {
    require(inRegion(locus), "can't fetch base outside region")
    bases(locus - 1).toChar.toUpper
  }

  def refBases = {
    bases.slice(start - 1, stop)
  }

  def refSlice(a: Int, b: Int) = bases.slice(index(a), index(b))
  def refBases(a: Int, b: Int) = refSlice(a, b).toArray

  def genomeSubString(locus: Int, n: Int) = {
    require(inRegion(locus) && inRegion(locus + n - 1))
    GenomeRegion.baseString(bases.slice(locus - 1, locus + n - 1))
  }

  def apply(i: Int) = refBase(start + i)
}