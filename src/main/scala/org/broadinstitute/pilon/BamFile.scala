/*
 * Copyright (c) 2013. The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc.
 * All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package org.broadinstitute.pilon

import collection.mutable.Map
import java.io.File
import scala.collection.JavaConversions._
import net.sf.samtools._

object BamFile {
  val indexSuffix = ".bai"
  val maxInsertSizes = Map(('frags -> 500), ('jumps -> 10000), ('unpaired -> 5000))
}

class BamFile(val bamFile: File, val bamType: Symbol) {
  val path = bamFile.getCanonicalPath()

  def reader = {
    val r = new SAMFileReader(bamFile, new File(path + BamFile.indexSuffix))
    r.setValidationStringency(SAMFileReader.ValidationStringency.SILENT)
    r
  }

  {
    // validate at object creation time by opening file and reader
	val r = reader
    require(r.hasIndex(), path + " must be indexed BAM")
	r.close
  }
  
  def getSeqs = {
	val r = reader
	val seqs = r.getFileHeader.getSequenceDictionary.getSequences
    require(r.hasIndex(), path + " must be indexed BAM")
	r.close
	seqs map { _.getSequenceName } toSet
  }
  
  lazy val header = {
    val r = reader
    val h = r.getFileHeader()    
    r.close
    h
  }
  lazy val sequenceDictionary = header.getSequenceDictionary()
  lazy val readGroups = header.getReadGroups()
  
  def printIndexStats() = BAMIndexMetaData.printIndexStats(bamFile)


  def getUnaligned = {
    val r = reader
    val reads = r.queryUnmapped.toList
    r.close
    reads
  }
  
  def process(region: GenomeRegion, printInterval: Int = 100000) : Unit = {

	val pileUpRegion = region.pileUpRegion
	
    // Ugh. Is there a way for samtools to do the right thing here and get
    // the pairs for which only one end overlaps?  Adding 10k slop until
    // that's figured out.
    //pileUpRegion.addReads(bam.reader.queryOverlapping(name, start, stop))
	val r = reader
    val reads = r.queryOverlapping(region.name, 
        (region.start-10000) max 0, (region.stop+10000) min region.contig.length)
    val readsBefore = pileUpRegion.readCount
    val covBefore = pileUpRegion.coverage
    var lastLoc = 0
    val huge = 10 * BamFile.maxInsertSizes(bamType)
    
    for (read <- reads) {
    	val loc = read.getAlignmentStart
    	if (Pilon.verbose && printInterval > 0 && loc > lastLoc + printInterval) {
    	  lastLoc = printInterval * (loc / printInterval)
    	  print("..." + lastLoc)
    	}
    	if ((!Pilon.pf) || (!read.getReadFailsVendorQualityCheckFlag)) {
    	  val insertSize = pileUpRegion.addRead(read, region.contigBases)
    	  if (insertSize > huge) {
    		if (Pilon.debug) println("WARNING: huge insert " + insertSize + " " + r)
    	  } 
    	  if (insertSize > 0 && insertSize <= huge) addInsert(insertSize)
    	}
    }
    r.close
    val meanCoverage = pileUpRegion.coverage - covBefore
    val nReads = pileUpRegion.readCount - readsBefore
    println(" Reads: " + nReads + ", Coverage: " + meanCoverage + ", Insert Size: " + insertSizeStats)
  }
  
  var insertSizeSum = 0.0
  var insertSizeCount = 0.0
  var insertSizeSumSq = 0.0

  def addInsert(insertSize: Int) = {
    insertSizeCount += 1
    insertSizeSum += insertSize
    insertSizeSumSq += insertSize * insertSize
  }
  
  def insertSizeMean = if (insertSizeSum > 0) (insertSizeSum / insertSizeCount) else 0.0
  
  def insertSizeSigma = {
    if (insertSizeSum > 0) {
      val mean = insertSizeMean
      scala.math.sqrt(((insertSizeSumSq / insertSizeCount) - (mean * mean)).abs)
    } else 0.0 
  }
  
  def insertSizeStats = "%.0f+/-%.0f".format(insertSizeMean, insertSizeSigma)

  def maxInsertSize = {
    """Returns reasonable max insert size for this bam, either computed or defaulted"""
    if (insertSizeCount >= 1000) (insertSizeMean + 3 * insertSizeSigma).round.toInt
    else BamFile.maxInsertSizes(bamType)    
  }
  
  class MateMap(reads: Seq[SAMRecord] = Nil) {
    val readMap = Map[String, SAMRecord]()
    val mateMap = Map[SAMRecord, SAMRecord]()
    var n = 0
    
    addReads(reads)
    
    def addReads(reads: Seq[SAMRecord]) = for (r <- reads) addRead(r)
    
    def addRead(read: SAMRecord) = {
      val readName = read.getReadName
      if (readMap contains readName) {
        val mate = readMap(readName)
        readMap -= readName
        mateMap += read -> mate
        mateMap += mate -> read
      } else readMap += readName -> read
      n += 1
    }
    
    def addStrays = {
      for (mate <- findStrays) {
        val name = mate.getReadName
        val read = readMap(name)
 	      readMap -= name
   	    mateMap += read -> mate
   	    mateMap += mate -> read
   	  }
    }
    
    def findStrays = {
      var nStrays = 0
      val mates = (readMap map {pair => strayMateMap.lookup(pair._2)}) filter {_ != null} toList;
      if (Pilon.debug) println("findStrays: " + mates.length)
      mates
    }
    
    def lookup(read: SAMRecord): SAMRecord = mateMap.getOrElse(read, null)
    
    def printDebug = println("mm: " + readMap.size + "/" + mateMap.size/2)
  }
  
  val strayMateMap = new MateMap()
  
  def mateMap(reads: Seq[SAMRecord]) = new MateMap(reads).mateMap
  
  def buildStrayMateMap = {
    val r = reader
   	// if it's not a proper pair but both ends are mapped, it's a stray mate
    print("Mapping stray pairs in " + bamFile + "...")
    for (read <- r.iterator)
      if (!(read.getProperPairFlag | read.getReadUnmappedFlag | read.getMateUnmappedFlag))
        strayMateMap.addRead(read)
    strayMateMap.printDebug
    r.close
    strayMateMap
  }
  
  def readsInRegion(region: Region) = {
    val r = reader
    val reads = r.queryOverlapping(region.name, region.start, region.stop).toList
    r.close
    reads
  }

  def recruitFlankReads(region: Region) = {
    val flanks = flankRegion(region)
    var reads = readsInRegion(flanks) 
    if (Pilon.debug) println("readsInRegion flanks: " + flanks + " " + reads.length + " reads")
    if (Pilon.strays) {
      val mm = new MateMap(reads)
      if (Pilon.debug) mm.printDebug
      reads ++= mm.findStrays
      if (Pilon.debug) mm.printDebug
    }
    reads      
  }
  
  def flankRegion(region: Region) = {
    val flank = maxInsertSize
    new Region(region.name, 
              (region.start - flank) max 1, 
              (region.stop + flank))
  }
  
  def recruitBadMates(region: Region) = {
    val midpoint = region.midpoint
    val flanks = flankRegion(region)
    val mateMap = new MateMap(readsInRegion(flanks)) 

    if (Pilon.strays) {
      if (Pilon.debug) mateMap.printDebug
      mateMap.addStrays
      if (Pilon.debug) mateMap.printDebug 
    }
    
    val mm = mateMap.mateMap

    // Filter to find pairs where r1 is anchored and r2 is unmapped (we'll include r2)
    val mm2 = mm filter { pair =>  
      val (r1, r2) = pair
      val r1mapped = !(r1.getReadUnmappedFlag || r1.getProperPairFlag)
      val rc = r1.getReadNegativeStrandFlag
      val start = r1.getAlignmentStart
      val end = r1.getAlignmentEnd
      val before = start < midpoint
      val after = end > midpoint
      val r1dir = (before && !rc) || (after && rc)
      val inRegion = flanks.inRegion(start) || flanks.inRegion(end)
      r1mapped && r1dir && inRegion
    }
    if (Pilon.debug) 
      println("# Filtered jumps " + mm2.size + "/" + mm.size)
    mm2.values.toList
  }
  
  override def toString() = path
}
