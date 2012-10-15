package org.broadinstitute.pilon

import java.io.{File,PrintWriter,FileWriter}
import collection.mutable.Map

class BamRegion (val bam: BamFile, val region: GenomeRegion)
	extends Region(region.name, region.start, region.stop) {
  val contig = region.contig
  require(stop <= contig.length, "GenomeRegion stop point must be within contig")
  var bases = contig.getBases
  
  var meanReadLength = 0
  
  // disposition flags
  var confirmed = new Array[Boolean](size)
  var ambiguous = new Array[Boolean](size)
  var changed = new Array[Boolean](size)
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

  var minDepth = Pilon.minMinDepth
  
  // GC if feature of the genome, not the pileup, but we'll keep
  // a reference to it here to make track generation easier
  val gc = region.gc

  lazy val physCoverageDist = new NormalDistribution(physCoverage, 2)
  lazy val coverageDist = new NormalDistribution(coverage, 2)
  lazy val badCoverageDist = new NormalDistribution(badCoverage, 2)
  lazy val insertSizeDist = new NormalDistribution(insertSize, 2)
  lazy val weightedMqDist = new NormalDistribution(weightedMq, 2)
  
  val changeList = region.changeList
  

  def process : Unit = {

	val pileUpRegion = region.pileUpRegion
	
    // Ugh. Is there a way for samtools to do the right thing here and get
    // the pairs for which only one end overlaps?  Adding 10k slop until
    // that's figured out.
    //pileUpRegion.addReads(bam.reader.queryOverlapping(name, start, stop))
	val radius = bam.maxInsertSize
    val reads = bam.reader.queryOverlapping(name, 
        (start - radius) max 0, (stop + radius) min contig.length)
    val readsBefore = pileUpRegion.readCount
    val covBefore = pileUpRegion.coverage
    
    pileUpRegion.addReads(reads, bases)
    pileUpRegion.computePhysCov
    val meanCoverage = pileUpRegion.coverage - covBefore
    val nReads = pileUpRegion.readCount - readsBefore
    
    println(" Reads: " + nReads + ", Coverage: " + meanCoverage + ", minDepth: " + minDepth)
  }
  
  
  def identifyIssues = {
    println("# IdentifyIssues: " + bam.bamType.name + " " + region)
    for (i <- changeList) {
      printChange(i)
      printChangeVcf(i)
    }

    if (Pilon.verbose) {
      identifyIssueRegions
    }
    

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
          print(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase + " " +
        		  //pu.baseCount.sums(pu.baseIndex(cBase)) + " " + pu.baseCount.sums(pu.baseIndex(rBase)) + " " + 
        		  bam.bamType.name + ":" + pu + endLine)

        case 'ins =>
          print(name + " " + locus(i) + " " + kind.name + " " + "." + " " + pu.insertCall + " " +
        		  //pu.insertions + " " + (pu.count-pu.insertions) + " " + 
        		  bam.bamType.name + ":" + pu + endLine)
        case 'del =>
          print(name + " " + locus(i) + " " + kind.name + " " + pu.deletionCall + " " + "." + " " +
        		  //pu.deletions + " " + pu.baseCount.sums(pu.baseIndex(rBase)) + " " + 
        		  bam.bamType.name + ":" + pu + endLine)
        case 'amb =>
          if (Pilon.verbose && rBase != cBase)
            print(name + " " + locus(i) + " " + kind.name + " " + rBase + " " + cBase + " " +
            		bam.bamType.name + ":" + pu + endLine)
      }
  }

  def printChangeVcf(i: Int) = {
      val (kind, pu) = changeMap(i)
      //Vcf.record(region, i, pu)
  }

  
  def prettyPrintRegions(header: String, regions: List[Region]) = {
    if (regions != Nil) {
      val totalSize: Int = (regions map {_.size}).sum
      print("# " + header + " n=" + regions.size + " bases=" + totalSize)
      regions foreach {r => print("  " + r.start + "-" + r.stop)}
      println
    }
  }
  
  type FixList = List[(Int, String, String)]

  def fixIssues = {
	var fixList : FixList = List()

	for (i <- changeList) {
      val (kind, pu) = changeMap(i)
      val rBase = refBase(locus(i))
      val cBase = pu.baseCall.base
      kind match {
        case 'snp =>
          //newBases(i) = cBase.toByte
          fixList ::= (locus(i), rBase.toString, cBase.toString)
        case 'ins =>
          fixList ::= (locus(i), "", pu.insertCall)
        case 'del =>
          fixList ::= (locus(i), "", pu.deletionCall)
        case _ =>
      }
    }

	if (Pilon.fixList contains 'gaps) {
	  for (gap <- gaps) {
		val filler = new GapFiller(region)
    	val (start, ref, patch) = filler.fillGap(gap)
    	if (start > 0) {
    		val old = if (ref != "") ref else "."
    		val now = if (patch != "") patch else "."
    	    println(name + " " + start + " fix " + ref.size + " " + old + " " + patch.size + " " + now)
    	    fixList ::= (start, ref, patch)
    	}
      }    
	}

	if (Pilon.fixList contains 'breaks) {
      //val breaks = possibleBreaks
      val breaks = clippingRegions
      val filteredBreaks = breaks filter { ! _.nearAny(gaps) }

      for (break <- filteredBreaks) {
        val filler = new GapFiller(region)
    	val (start, ref, patch) = filler.fixBreak(break)
    	if (start > 0) {
    		val old = if (ref != "") ref else "."
    		val now = if (patch != "") patch else "."
    	    println(name + " " + start + " fix " + ref.size + " " + old + " " + patch.size + " " + now)
    	    fixList ::= (start, ref, patch)
    	}
      }
	}
	
    applyFixes(fixList)
  }
  
  def fixFixListX(inList: FixList, outList: FixList = Nil): FixList = {
    """Remove overlaps, keeping larger fix"""
    inList match {
      case fix1 :: fix2 :: tail => {
        val region1 = new Region(region.name, fix1._1, fix1._1 + fix1._2.length - 1)
        val region2 = new Region(region.name, fix2._1, fix2._1 + fix2._2.length - 1)
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
    var outList : FixList = Nil
    
    while (fixes != Nil) {
      //println(fixes.head)
      fixes match {
        case fix1 :: fix2 :: tail => {
          val region1 = new Region(region.name, fix1._1, fix1._1 + fix1._2.length - 1)
          val region2 = new Region(region.name, fix2._1, fix2._1 + fix2._2.length - 1)
          if (region1.overlaps(region2)) {
            //println("fixFixList overlap: " + (region1, region2))
            if (region1.size > region2.size)
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
  
  def applyFixes(fixList: FixList) = {
    var newBases = region.refBases.clone
    val sortedFixes = fixList.sortWith({(x, y) => x._1 < y._1})
    val fixedFixes = fixFixList(sortedFixes)
    for (fix <- fixedFixes ) {
      val (locus, was, patch) = fix
      val start = index(locus)
      val before = newBases.slice(0, start)
      val after = newBases.slice(start + was.length, newBases.length)
      val ref = newBases.slice(start, start + was.length).map(_.toChar).mkString("")
      if (ref != was) println("Fix mismatch: ref=" + ref + " was=" + was)
      newBases = before ++ (patch map {_.toByte}) ++ after
      //if (Pilon.debug) println("Fixing=" + was.length + " " + patch.length + " " + newBases.length)
    }
    newBases    
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
  

  def nearEdge(r: Region, radius: Int = 100) = r.start - region.start < radius || region.stop - r.stop < radius
  def deltaFraction(i: Int) = deltaCoverage(i) / coverageDist.mean
  def dipCoverage(i: Int, radius: Int = 100) = dip(i, coverage, radius)
  def dipFraction(i: Int) = dipCoverage(i) / coverageDist.mean

  def ambiguousRegions = summaryRegions({ i: Int => ambiguous(i)}) filter {r => r.start != r.stop}
  def changeRegions = summaryRegions({ i: Int => changed(i)})
  def gaps = summaryRegions({ i: Int => region.refBase(locus(i)) == 'N' })
  def highCopyNumberRegions = summaryRegions({ i: Int => copyNumber(i) > 1})
  def unConfirmedRegions = summaryRegions({ i: Int => !confirmed(i)})
  
  def lowCoverage(i: Int) = coverage(i) < coverageDist.moments(0) * .20
  def lowCoverageRegions = summaryRegions(lowCoverage)
  def highClipping(i: Int) = clips(i) >= coverage(i)
  def clippingRegions = summaryRegions(highClipping)
  
  def breakp(i: Int) =  dipFraction(i) >= 1.0
  def possibleBreaks = summaryRegions(breakp)

  def insertionp(i: Int) = breakp(i) && insertSizeDist.toSigma(insertSize(i)) <= -3.0
  def possibleInsertions = summaryRegions(insertionp)

  def deletionp(i: Int) = breakp(i) && insertSizeDist.toSigma(insertSize(i)) >= 3.0
  def possibleDeletions = summaryRegions(deletionp)
  
  def possibleCollapsedRepeats = {
    highCopyNumberRegions filter {r => 
      //println(r + " " + deltaFraction(index(r.start)) + " " + deltaFraction(index(r.stop)))
      deltaFraction(index(r.start)) >= 0.5 && deltaFraction(index(r.stop)) >= 0.5}
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


  def refBase(i: Int) = region.refBase(i)
  
  def ++(other: BamRegion) = {
    require(name == other.name)
    require(stop == other.start - 1)
    val gr = new BamRegion(bam, region)
    (coverage ++ other.coverage).copyToArray(gr.coverage)
    (physCoverage ++ other.physCoverage).copyToArray(gr.coverage)
    (badCoverage ++ other.badCoverage).copyToArray(gr.badCoverage)
    (insertSize ++ other.insertSize).copyToArray(gr.insertSize)
    (clips ++ other.clips).copyToArray(gr.clips)
    (copyNumber ++ other.copyNumber).copyToArray(gr.copyNumber)
    (confirmed ++ other.confirmed).copyToArray(gr.confirmed)
    (changed ++ other.changed).copyToArray(gr.changed)
    (ambiguous ++ other.ambiguous).copyToArray(gr.ambiguous)
    (weightedQual ++ other.weightedQual).copyToArray(gr.weightedQual)
    (weightedMq ++ other.weightedMq).copyToArray(gr.weightedMq)
    gr.changeMap = changeMap ++ other.changeMap
    gr
  }	
}
