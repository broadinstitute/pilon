package org.broadinstitute.pilon

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
  
  lazy val header = reader.getFileHeader()
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
    
    for (r <- reads) {
    	val loc = r.getAlignmentStart
    	if (Pilon.verbose && printInterval > 0 && loc > lastLoc + printInterval) {
    	  lastLoc = printInterval * (loc / printInterval)
    	  print("..." + lastLoc)
    	}
    	val insertSize = pileUpRegion.addRead(r, region.bases)
    	if (Pilon.debug && insertSize > huge)
          println("WARNING: huge insert " + insertSize + " " + r)
        if (insertSize > 0) addInsert(insertSize)
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
  
  def insertSizeMean = if (insertSizeSum > 0) (insertSizeSum / insertSizeCount) else 0
  
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
  
  override def toString() = path
}
