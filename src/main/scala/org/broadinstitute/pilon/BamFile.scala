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

import collection.mutable.{HashMap,Map}
import java.io.File
import scala.collection.JavaConverters._
import htsjdk.samtools._

object BamFile {
  val indexSuffix = ".bai"
  val maxInsertSizes = Map(('frags -> 500), ('jumps -> 10000), ('unpaired -> 10000), ('bam -> 10000))
  val minOrientationPct = 10
  val maxFragInsertSize = 700
  // long read type codes
  val notLongRead = 0
  val nanoporeLongRead = 1
  val pacbioLongRead = 2
}

class BamFile(val bamFile: File, var bamType: Symbol, val subType: Symbol = 'none) {
  import BamFile._
  val path = bamFile.getCanonicalPath()
  var baseCount: Long = 0

  val longReadType = if (bamType == 'unpaired) {
    if (subType == 'nanopore) nanoporeLongRead
    else if (subType == 'pacbio) pacbioLongRead
    else 0
  } else 0

  def reader = {
    //val r = new SAMFileReader(bamFile, new File(path + BamFile.indexSuffix))
    val r = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile)
    r
  }

  {
    // validate at object creation time by opening file and reader
    val r = reader
    require(r.hasIndex, path + " must be indexed BAM")
    r.close
  }

  def getSeqs = {
    // returns an array of sequence records indexed by bam seq index
    val r = reader
    val seqs = r.getFileHeader.getSequenceDictionary.getSequences.asScala
    r.close
    val seqArray = new Array[SAMSequenceRecord](seqs.length)
    for (s <- seqs) seqArray(s.getSequenceIndex) = s
    seqArray
  }

  def getSeqNames = {
    val seqs = getSeqs
    seqs.map({ _.getSequenceName }).toSet
  }

  def validateSeqs(seqsOfInterest: Set[String]) = {
    require(!getSeqNames.intersect(seqsOfInterest).isEmpty, bamFile + " doesn't have sequence for any of " + seqsOfInterest.mkString(", "))
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
    val reads = r.queryUnmapped.asScala.filter(validateRead(_)).toList
    r.close
    reads
  }

  def validateRead(read: SAMRecord) = {
    (Pilon.nonPf || !read.getReadFailsVendorQualityCheckFlag) &&
      (Pilon.duplicates || !read.getDuplicateReadFlag) &&
    !read.isSecondaryAlignment
  }

  
  def process(region: GenomeRegion, printInterval: Int = 100000) = {
    val pileUpRegion = region.pileUpRegion

    // Ugh. Is there a way for samtools to do the right thing here and get
    // the pairs for which only one end overlaps?  Adding 10k slop until
    // that's figured out.
    //pileUpRegion.addReads(bam.reader.queryOverlapping(name, start, stop))
    val longRead = (bamType == 'unpaired) && (subType == 'nanopore || subType == 'pacbio)
    val r = reader
    val reads = r.queryOverlapping(region.name,
      (region.start-10000) max 0, (region.stop+10000) min region.contig.length).asScala
    val readsBefore = pileUpRegion.readCount
    val baseCountBefore = pileUpRegion.baseCount
    val covBefore = pileUpRegion.coverage
    var lastLoc = 0
    val huge = 10 * BamFile.maxInsertSizes(bamType)

    for (read <- reads) {
      val loc = read.getAlignmentStart
      if (Pilon.verbose && printInterval > 0 && loc > lastLoc + printInterval) {
        lastLoc = printInterval * (loc / printInterval)
        print("..." + lastLoc)
      }
      if (validateRead(read)) {
        val insertSize = pileUpRegion.addRead(read, region.contigBases, longReadType)
        if (insertSize > huge) {
          //if (Pilon.debug) println("WARNING: huge insert " + insertSize + " " + read)
        }
        addInsert(insertSize, read.getReadNegativeStrandFlag)
      }
    }
    r.close

    val meanCoverage = pileUpRegion.coverage - covBefore
    val nReads = pileUpRegion.readCount - readsBefore

    // Track baseCount for global coverage
    baseCount += pileUpRegion.baseCount - baseCountBefore
    meanCoverage
  }

  var mapped = 0
  var unmapped = 0
  var filtered = 0
  var secondary = 0
  var proper = 0

  class InsertSizeStats {
    var count: Long = 0
    var sum: Long = 0
    var sumSq: Long = 0

    def add(size: Int) = {
      sum += size
      sumSq += size * size
      count += 1
    }

    def mean = if (count > 0) (sum.toDouble / count.toDouble) else 0.0

    def sigma = {
      if (sum > 0) {
        scala.math.sqrt(((sumSq.toDouble / count.toDouble) - (mean * mean)).abs)
      } else 0.0
    }

    def reset = {
      count = 0
      sum = 0
      sumSq = 0
    }

    def maxInsertSize = {
      if (count > 1000) (mean + 3.0 * sigma).round.toInt else 0
    }

    override def toString = "%.0f+/-%.0f".format(mean, sigma)
  }

  val insertStatsFR = new InsertSizeStats()
  val insertStatsRF = new InsertSizeStats()
  val insertStatsUnpaired = new InsertSizeStats()

  val huge = 5 * BamFile.maxInsertSizes(bamType)

  def addInsert(insertSize: Int, rc: Boolean = false, unpaired: Boolean = false) = {
    val fr = (insertSize > 0) ^ rc
    //if (Pilon.debug) println("i=" + insertSize + " rc=" + rc)
    if (insertSize > 0 && insertSize < huge) {
      if (unpaired) insertStatsUnpaired.add(insertSize)
      else if (fr) insertStatsFR.add(insertSize)
      else insertStatsRF.add(insertSize)
    }
  }

  def pctFR = Utils.pct(insertStatsFR.count, insertStatsFR.count + insertStatsRF.count + insertStatsUnpaired.count)
  def pctRF = Utils.pct(insertStatsRF.count, insertStatsFR.count + insertStatsRF.count + insertStatsUnpaired.count)
  def pctUnpaired = Utils.pct(insertStatsUnpaired.count, insertStatsFR.count + insertStatsRF.count + insertStatsUnpaired.count)

  def insertSizeMean = {
    if (pctFR >= 50) insertStatsFR.mean else insertStatsRF.mean
  }

  def insertSizeSigma = {
    if (pctFR >= 50) insertStatsFR.sigma else insertStatsRF.sigma
  }

  def maxInsertSize = {
    // Returns reasonable max insert size for this bam, either computed or defaulted
    val maxFromStats = insertStatsFR.maxInsertSize max insertStatsRF.maxInsertSize max insertStatsUnpaired.maxInsertSize
    if (maxFromStats > 0) maxFromStats else BamFile.maxInsertSizes(bamType)
  }

  class MateMap(reads: Seq[SAMRecord] = Nil) {
    type ReadMap = HashMap[String, SAMRecord]
    //val readMap: ReadMap = HashMap()
    //val mateMap = Map[SAMRecord, SAMRecord]()
    val readMap1: ReadMap = HashMap()
    val readMap2: ReadMap = HashMap()
    var n = 0
    
    addReads(reads)
    
    def addReads(reads: Seq[SAMRecord]) = for (r <- reads) addRead(r)
    
    def addRead(read: SAMRecord) = {
      val readName = read.getReadName
      val readMap = if (read.getFirstOfPairFlag) readMap1 else readMap2
      readMap(readName) = read
    }

    def addStrays = {
      for (r <- findStrays) addRead(r)
    }
    
    def findStrays = {
      var strays: List[SAMRecord] = Nil

      def findMates(rm1: ReadMap, rm2: ReadMap) = {
        for ((name, read) <- rm1; if (!(rm2 contains name))) {
          val mate = strayMateMap.lookupMate(read)
          if (mate != null) strays ::= mate
        }
      }
      findMates(readMap1, readMap2)
      findMates(readMap2, readMap1)
      if (Pilon.debug) println("findStrays: " + strays.length)
      strays
    }
    
    def lookupMate(read: SAMRecord): SAMRecord = {
      // Look for the mate int the other readMap
      val mateReadMap = if (read.getFirstOfPairFlag) readMap2 else readMap1
      mateReadMap.getOrElse(read.getReadName, null)
    }

    def pairs() = {
      var pairList = List[(SAMRecord, SAMRecord)]()
      for ((name, read1) <- readMap1; if (readMap2 contains name)) {
        val read2 = readMap2(name)
        pairList ::= (read1, read2)
        pairList ::= (read2, read1)
      }
      pairList
    }


    def nStrays = pairs.length

    def printDebug = println("mm: " + readMap1.size + "+" + readMap2.size + "=" + nStrays/2)
  }
  
  val strayMateMap = new MateMap()

  def autoBam(): Symbol = {
    val fr = pctFR
    val rf = pctRF
    val un = pctUnpaired

    if (un >= fr && un >= rf) 'unpaired
    else {
      val insertSize = if (rf > fr) insertStatsRF.mean else insertStatsFR.mean
      if (insertSize >= maxFragInsertSize) 'jumps else 'frags
    }
  }
  
  def scan(seqsOfInterest: Set[String]) = {
    val r = reader
    for (read <- r.iterator.asScala) {
      if (!validateRead(read)) filtered += 1
      else if (read.getReadUnmappedFlag) unmapped += 1
      else {
        mapped += 1
        if (read.getReadPairedFlag) {
          val pp = read.getProperPairFlag
          if (pp) {
            proper += 1
            addInsert(read.getInferredInsertSize, read.getReadNegativeStrandFlag)
          } else if (Pilon.strays && !read.getMateUnmappedFlag) {
            // If it's not a proper pair but both ends are mapped, it's a stray mate.
            // We'll track them iff they are among the seqs we are processing.
            if ((seqsOfInterest contains read.getReferenceName)
              || (seqsOfInterest contains read.getMateReferenceName))
              strayMateMap.addRead(read)
          }
        } else {
          addInsert(read.getReadLength, read.getReadNegativeStrandFlag, true)
        }
      }
    }
    r.close

    val totalReads = mapped + unmapped + filtered
    var summary = bamFile + ": %d reads, %d filtered, %d mapped, %d proper".format(totalReads, filtered, mapped, proper)
    if (Pilon.strays) summary += ", " + strayMateMap.nStrays + " stray"
    val insertCount = insertStatsFR.count + insertStatsRF.count
    if (pctFR >= minOrientationPct)
      summary += ", FR " + pctFR + "% " + insertStatsFR
    if (pctRF >= minOrientationPct)
      summary += ", RF " + pctRF + "% " + insertStatsRF
    if (pctUnpaired >= minOrientationPct)
      summary += ", Unpaired " + pctUnpaired + "% " + insertStatsUnpaired
    summary += ", max " + maxInsertSize
    if (bamType == 'bam) {
      bamType = autoBam
      summary += " " + bamType.name
    }
    println(summary)
  }
  
  def readsInRegion(region: Region) = {
    val r = reader
    val reads = r.queryOverlapping(region.name, region.start, region.stop).asScala.filter(validateRead(_)).toList
    r.close()
    reads
  }

  def recruitFlankReads(region: Region) = {
    val flanks = flankRegion(region)
    var reads = readsInRegion(flanks) 
    if (Pilon.debug) println("readsInRegion flanks: " + flanks + " " + reads.length + " reads")
    if (Pilon.strays && bamType != 'unpaired) {
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
    
    // Filter to find pairs where r1 is anchored and r2 is unmapped (we'll include r2)
    val frPct = pctFR

    var mates: List[SAMRecord] = Nil
    for ((r1, r2) <- mateMap.pairs) {
      if (!(r1.getReadUnmappedFlag || r1.getProperPairFlag)) {
        val rc = r1.getReadNegativeStrandFlag
        val start = r1.getAlignmentStart
        val end = r1.getAlignmentEnd
        val before = start < midpoint
        val after = end > midpoint
        val goodOrientation = {
          val frOrientation = (before && !rc) || (after && rc)
          // If we are almost all FR or RF, use only appropriately mapped anchors
          if (frPct > 100 - minOrientationPct) frOrientation
          else if (frPct < minOrientationPct) !frOrientation
          // Otherwise, either orientation could be useful
          else true
        }
        if (goodOrientation && (flanks.inRegion(start) || flanks.inRegion(end)))
        mates ::= r2
      }
    }
    if (Pilon.debug)
      println("# Filtered jumps " + mates.size + "/" + mateMap.nStrays)
    mates
  }
  
  override def toString() = path
}

