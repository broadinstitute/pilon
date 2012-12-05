package org.broadinstitute.pilon

import scala.collection.JavaConversions._
import net.sf.samtools._

object PileUpRegion {
  val defaultQual: Byte = 15
}

class PileUpRegion(name: String, start: Int, stop: Int)
  extends Region(name, start, stop) {

  val pileups = new Array[PileUp](size)
  for (i <- 0 to size - 1) pileups(i) = new PileUp

  var baseCount: Long = 0
  var readCount = 0
  val trustedFlank = Pilon.flank

  def coverage = baseCount / size

  def add(locus: Int, base: Char, qual: Int, mq: Int,
    pair: Boolean) = {
    if (inRegion(locus)) {
      if (pair) {
        pileups(index(locus)).add(base, qual, mq)
        baseCount += 1
      } else {
        if (!pair) pileups(index(locus)).badPair += 1
      }
    }
  }

  def remove(locus: Int, base: Char, qual: Int, mq: Int,
    pair: Boolean) = {
    if (inRegion(locus)) {
      if (pair) {
        pileups(index(locus)).remove(base, qual, mq)
        baseCount -= 1
      } else {
        if (!pair) pileups(index(locus)).badPair -= 1
      }
    }
  }
  var physCovStart = 0
  var insertSizeStart = 0

  def physCovIncr(aStart: Int, aEnd: Int, iSize: Int, paired: Boolean, valid: Boolean) = {
    if (!valid || (paired && iSize <= 0)) 0
    else {
      val (start, end) =
        if (!paired) {
          (aStart min aEnd, aStart max aEnd)
        } else if (iSize > 0) {
          (aStart, aStart + iSize)
        } else {
          (aEnd + 1 + iSize, aEnd + 1) // iSize negative!
        }
      val insertSize = end - start
      //println("as=" + aStart + " ae=" + aEnd + " iSize=" + iSize + " p=" + paired + " v=" + valid + " s=" + start + " e=" + end)
      if (inRegion(start)) {
        pileups(index(start)).physCov += 1
        pileups(index(start)).insertSize += insertSize
      } else if (beforeRegion(start) && !beforeRegion(end)) {
        physCovStart += 1
        insertSizeStart += insertSize
      }
      if (inRegion(end)) {
        pileups(index(end)).physCov -= 1
        pileups(index(end)).insertSize -= insertSize
      }
      insertSize
    }
  }

  def computePhysCov = {
    pileups(0).physCov += physCovStart
    pileups(0).insertSize += insertSizeStart
    for (i <- 1 until pileups.length) {
      pileups(i).physCov += pileups(i - 1).physCov
      pileups(i).insertSize += pileups(i - 1).insertSize
    }
    for (i <- 0 until pileups.length)
      if (pileups(i).physCov > 0)
        pileups(i).insertSize /= pileups(i).physCov
  }

  def addRead(r: SAMRecord, refBases: Array[Byte]) = {
    val length = r.getReadLength
    val bases = r.getReadBases
    val paired = r.getReadPairedFlag
    val valid = (!paired || r.getProperPairFlag)
    val insert = r.getInferredInsertSize
    val aStart = r.getAlignmentStart
    val aEnd = r.getAlignmentEnd
    val mq = r.getMappingQuality
    val cigar = r.getCigar
    var readOffset = 0
    var refOffset = 0
    val quals = if (r.getBaseQualities.size > 0) r.getBaseQualities
    else {
      val q = new Array[Byte](length)
      for (i <- 0 until length) q(i) = PileUpRegion.defaultQual
      q
    }

    def baseString(bases: Array[Byte]) = bases map { _.toChar } mkString ""
    def trusted(offset: Int) = offset >= trustedFlank && length - trustedFlank > offset

    // parse read alignment and add to pileups
    for (i <- 0 until cigar.numCigarElements) {
      val ele = cigar.getCigarElement(i)
      val len = ele.getLength
      val op = ele.getOperator
      val locus = aStart + refOffset
      op match {
        case CigarOperator.I =>
          var insertion = bases.slice(readOffset, readOffset + len)
          val istr = baseString(insertion)
          var iloc = locus
          var rloc = readOffset
          if (valid && trusted(readOffset) && inRegion(iloc)) {
            while (iloc > 1 && refBases(iloc - 2) == insertion(len - 1)) {
              iloc -= 1
              insertion = insertion.slice(len - 1, len) ++ insertion.slice(0, len - 1)
            }
            pileups(index(iloc)).addInsertion(insertion, quals(readOffset))
          }
        case CigarOperator.D =>
          var dloc = locus
          var rloc = readOffset
          if (valid && trusted(readOffset) && inRegion(dloc) && inRegion(dloc+len-1)) {
            while (dloc > 1 && rloc > 0 && refBases(dloc - 2) == refBases(dloc + len - 2)) {
              dloc -= 1
              rloc -= 1
              val base = bases(rloc).toChar
              val qual = quals(rloc)
              // as we slide the deletion, remove old base from pileup and
              // add to end
              if (trusted(rloc) && inRegion(dloc)) {
                remove(dloc, base, qual, mq, valid)
                if (inRegion(dloc + len))
                  add(dloc + len, base, qual, mq, valid)
              }
            }
            pileups(index(dloc)).addDeletion(refBases.slice(dloc - 1, dloc + len - 1), quals(readOffset))
          }
        case CigarOperator.M =>
          for (i <- 0 until len) {
            val rOff = readOffset + i
            if (trusted(rOff)) {
              val locusPlus = locus + i
              val base = bases(rOff).toChar
              val qual = quals(rOff)
              if (inRegion(locusPlus))
                add(locusPlus, base, qual, mq, valid)
            }
          }
        case CigarOperator.S =>
          val clipStart = if (readOffset == 0) locus - len else locus
          val clipEnd = clipStart + len - 1
          if (inRegion(clipStart)) pileups(index(clipStart)).clips += 1
          if (inRegion(clipEnd)) pileups(index(clipEnd)).clips += 1
          for (i <- 0 until len) {
            val rOff = readOffset + i
            val locusPlus = clipStart + i
            val base = bases(rOff).toChar
            val qual = quals(rOff)
            if (inRegion(locusPlus))
              add(locusPlus, base, qual, mq, false)
          }
        case _ =>
          println("unknown cigar op=" + op + " in " + r.getCigarString)
      }
      if (op.consumesReadBases()) readOffset += len
      if (op.consumesReferenceBases()) refOffset += len
    }

    readCount += 1
    physCovIncr(aStart, aEnd, insert, paired, valid)
  }

  def addReads(reads: SAMRecordIterator, refBases: Array[Byte], printInterval: Int = 100000) = {
    var lastLoc = 0
    for (r <- reads) {
      val loc = addRead(r, refBases)
      if (Pilon.verbose && printInterval > 0 && loc > lastLoc + printInterval) {
        lastLoc = printInterval * (loc / printInterval)
        print("..." + lastLoc)
      }
    }
  }

  def dump = {
    for (i <- 0 to size - 1) println(locus(i) + ": " + pileups(i))
  }

  def callRegion = {
    for (i <- 0 until size) {
      var bc = pileups(i).baseCall
      if (bc.q < 20 && !bc.homo) println(locus(i) + ": bc=" + bc + " pu=" + pileups(i))
    }
  }

  def postProcess = {
    """To be called after all reads have been added to pileUps"""
    computePhysCov
  }

  def meanReadLength = (baseCount / readCount) toInt

  def apply(i: Int) = pileups(i)
}