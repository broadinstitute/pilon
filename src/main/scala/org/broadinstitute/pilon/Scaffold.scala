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

/**
 * Created with IntelliJ IDEA.
 * User: bruce
 * Date: 10/17/13
 * Time: 2:35 PM
 * To change this template use File | Settings | File Templates.
 */


import scala.collection.JavaConversions._
import collection.mutable.{ Map, HashMap, Set, HashSet }
import htsjdk.samtools._


class MatePair(r1: SAMRecord, r2: SAMRecord) {
  // coordinates of pair, ordered lowest first
  val (scaffold1, coord1, scaffold2, coord2) = {
    val s1 = scaffold(r1)
    val s2 = scaffold(r2)
    val c1 = coord(r1)
    val c2 = coord(r2)
    if ((s1 < s2) || (s1 == s2 && c1 < c2)) (s1, c1, s2, c2)
    else (s2, c2, s1, c1)
  }
  val mq = r1.getMappingQuality min r2.getMappingQuality

  def sameScaffold = (scaffold1 == scaffold2)

  def scaffold(read: SAMRecord): Int = {
    // make 1-based, and note negative means rc, positive means fw
    //val s = read.getReferenceIndex + 1
    //if (read.getReadNegativeStrandFlag) -s else s
    read.getReferenceIndex
  }

  def coord(read: SAMRecord) = {
    // Note that fw coords are negative, rc are postitive
    // so fw + rc gives insertSize in normal orientation
    if (read.getReadNegativeStrandFlag) read.getAlignmentEnd
    else -read.getAlignmentStart
  }

  // 64bit coord: High order word is scaffold (negative if rc, low order is signed coord
  def mask32bit: Long = (1.toLong << 32) - 1
  def mask16bit: Long = (1.toLong << 16) - 1

  def longCoord(s: Int, c: Int): Long = (s.toLong << 32) | (c.toLong & mask32bit)

  def longCoord1 = longCoord(scaffold1, coord1)
  def longCoord2 = longCoord(scaffold2, coord2)


  def distance : Long = {
    // A simple expression designed to yield insertSize in the case of
    // same-scaffold FR orientation, small negative for RF "innies",
    // but something big and different
    // for different intra-scaffold combinations, still keeping neighborhoods
    // close together.
    // We muck with high bits for intra-scaffold
    //((scaffold1 ^ scaffold2) << 24) ^ (coord1 + coord2)
    if (sameScaffold)
      (coord1 + coord2)
    else
      ((scaffold1 & mask16bit).toLong << 48) | ((scaffold2 & mask16bit) << 32) | ((coord1 + coord2) & mask32bit)
  }

  def isRunt(maxInsert : Int = 10000) = {
    if (sameScaffold) {
      val d = distance
      (d < maxInsert) && (d > -maxInsert)
    } else false
  }

  def ambiguousPlacement = mq < 4

  override def toString = "<MatePair %d:%d %d:%d %d %d".format(scaffold1, coord1, scaffold2, coord2,
    mq, distance)
}

class LinkCluster(val matePairs: Array[MatePair], scaffolds: Array[SAMSequenceRecord], sigma: Int) {
  val nLinks = matePairs.size
  val scaffold1 = matePairs.map(_.scaffold1).min
  val minCoord1 = matePairs.map(_.coord1).min
  val maxCoord1 = matePairs.map(_.coord1).max
  val scaffold2 = matePairs.map(_.scaffold2).min
  val minCoord2 = matePairs.map(_.coord2).min
  val maxCoord2 = matePairs.map(_.coord2).max
  val mq = matePairs.map(_.mq).sum / matePairs.length
  val seq1 = scaffolds(scaffold1)
  val seq2 = scaffolds(scaffold2)
  val size1 = seq1.getSequenceLength
  val size2 = seq2.getSequenceLength

  //println("LinkCluster " + this)
  //dumpCoords(matePairs)
  require(scaffold1 == matePairs.map(_.scaffold1).max, "Not all scaffold1 the same!")
  require(scaffold2 == matePairs.map(_.scaffold2).max, "Not all scaffold2 the same!")

  def sameScaffold = (scaffold1 == scaffold2)

  def sameOrientation = (minCoord1 < 0 && minCoord2 < 0) || (minCoord1 > 0 && minCoord2 > 0)

  def spread1 = maxCoord1 - minCoord1
  def spread2 = maxCoord2 - minCoord2

  def valid = (spread1 > sigma) && (spread2 > sigma)

  def circular = {
    valid && sameScaffold && nearEnds
  }

  def reportCoord(name: String, coord: Int) = {
    val dir = if (coord < 0) "+" else "-"
    "%s:%d%s".format(name, coord.abs, dir)
  }

  def reportCoord1 = reportCoord(seq1.getSequenceName, minCoord1)
  def reportCoord2 = reportCoord(seq2.getSequenceName, minCoord2)

  def reportCoords = {
    (reportCoord1, reportCoord2)
  }

  def reportCircular = {
    assert(valid && circular, "don't call me unless you now I'm circular")
    val name = seq1.getSequenceName
    ""
  }

  def nearEnd(coord: Int, size: Int) = {
    (coord < sigma) || (-coord > size - sigma)
  }

  def nearEnds = nearEnd(minCoord1, size1) && nearEnd(minCoord2, size2)

  def scaffoldLink = {
     valid && nearEnds && !sameScaffold
  }

  def rearrangement = {
    valid && !circular && !scaffoldLink
  }

  override def toString = {
    //var str = "<LinkCluster %d %d:%d+%d %d:%d+%d %d ".format(nLinks, scaffold1, minCoord1, spread1,
    //  scaffold2, minCoord2, spread2, mq)
    var str = "<LinkCluster %d %d:%d%+d %d:%d%+d %d".format(nLinks, scaffold1, minCoord1, maxCoord1,
      scaffold2, minCoord2, maxCoord2, mq)
    if (valid) str += " valid"
    if (circular) str += " circular"
    str + ">"
  }


  def dumpCoords(coords: Array[MatePair]) = coords foreach println
}

object Scaffold {

  def findClusters(coords: Array[MatePair], scaffolds: Array[SAMSequenceRecord],
                   sigma: Int, minCluster: Int = 10) = {
    val windowSize = 4 * sigma
    val clustersByDistance = findClustersInternal(coords, windowSize, minCluster, {_.distance})
    val clusters = clustersByDistance map {findClustersInternal(_, windowSize, minCluster, {_.longCoord1})}

    val linkClusters = clusters.flatten.map({new LinkCluster(_, scaffolds, sigma)}).sortWith({_.nLinks > _.nLinks})
    if (Pilon.debug) linkClusters foreach println
    linkClusters
  }

  // This is used to find clusters within a neighborhood along some dimension (mpFunc)
  def findClustersInternal(unsortedCoords: Array[MatePair],
                           window: Int,
                           minCluster: Int = 10,
                           mpFunc: MatePair => Long) = {
    if (Pilon.debug)
      println("finding clusters: size=%d window=%d min=%d func=%s".format(unsortedCoords.size,
        window, minCluster, mpFunc))

    val coords = unsortedCoords.sortWith({mpFunc(_) < mpFunc(_)})
    if (Pilon.debug) dumpCoords(coords)
    var best = (0, 0)
    var clusterRanges : List[(Int, Int)] = Nil
    for (tail <- 0 until coords.size) {
      val windowLimit = mpFunc(coords(tail)) + window
      var head = tail + 1
      while (head < coords.size && mpFunc(coords(head)) < windowLimit)
        head += 1
      if (head - tail > minCluster) {
        if (head - tail > best._2 - best._1)
          best = (tail, head)
      } else {
        if (best._2 > 0) clusterRanges ::= best
        best = (0, 0)
      }
    }
    if (best._2 > 0) clusterRanges ::= best
    //clusterRanges foreach {c => println("cluster " + (c._2-c._1) + " " + c._1 + "-" + c._2)}
    clusterRanges map {c => coords.slice(c._1, c._2)}
  }

  def analyzeStrays(bam: BamFile) = {
    val mm = bam.strayMateMap.mateMap
    val sigma = bam.insertSizeSigma
    val scaffolds = bam.getSeqs
    val scaffoldSizes = scaffolds.map({_.getSequenceLength})
    val genomeSize = scaffoldSizes.sum
    println("Analyzing large-scale structure using " + bam)
    if (Pilon.verbose) {
      println("analyzing strays in " + bam)
      println("genome size " + genomeSize)
      println("max imsert " + bam.maxInsertSize)
    }

    var innie = 0
    var intra = 0
    var ambig = 0
    var links: List[MatePair] = Nil
    for ((r1,r2) <- mm if (r2.getSecondOfPairFlag)) {
      val mp = new MatePair(r1, r2)
      if (mp.isRunt(/* bam.maxInsertSize */10000)) innie += 1
      else if (mp.ambiguousPlacement) ambig += 1
      else {
        if (mp.sameScaffold) intra += 1
        links ::= mp
      }
    }
    val nLinks = links.length
    val backgroundRate = (nLinks.toFloat * 4 * sigma / genomeSize.toFloat).toInt
    if (Pilon.debug)
      println("ambig=" + ambig + " same=" + intra + " innies=" + innie + " links=" + nLinks + " background=" + backgroundRate)
    val coords = links.toArray

    val clusters = findClusters(links.toArray, scaffolds, sigma.toInt, /*TODO: 1 * backgroundRate*/ 25)
    for (c <- clusters if c.valid) {
      val (c1, c2) = c.reportCoords
      if (c.circular)
        print("Circular element " + scaffolds(c.scaffold1).getSequenceName)
      else if (c.scaffoldLink) {
        print("Candidate scaffold link " + c1 + " to " + c2)
      }
      else if (c.rearrangement) {
        print("Candidate rearrangement " + c1 + " connects to " + c2)
      }
      else print("What is this? " + c)
      if (c.sameOrientation) print(" reversed")
      println(" (%d supporting links)".format(c.nLinks))
    }
  }


  def dumpCoords(coords: Array[MatePair]) = coords foreach println
}
