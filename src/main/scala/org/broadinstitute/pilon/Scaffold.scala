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
import net.sf.samtools._


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


  // convert back to (scaffold, coord) pair of ints
  def decodeLongCoord(c: Long) = ((c >> 32).toInt, (c & mask32bit).toInt)

  def coordsAndDistance= (longCoord1, longCoord2, distance)

  override def toString = "<MatePair %d:%d %d:%d %d %d".format(scaffold1, coord1, scaffold2, coord2,
    mq, distance)
}

class LinkCluster(val matePairs: Array[MatePair]) {
  val nLinks = matePairs.size
  val scaffold1 = matePairs.map(_.scaffold1).min
  val minCoord1 = matePairs.map(_.coord1).min
  val maxCoord1 = matePairs.map(_.coord1).max
  val scaffold2 = matePairs.map(_.scaffold2).min
  val minCoord2 = matePairs.map(_.coord2).min
  val maxCoord2 = matePairs.map(_.coord2).max
  val mq = matePairs.map(_.mq).sum / matePairs.length

  //println("LinkCluster " + this)
  //dumpCoords(matePairs)
  assert(scaffold1 == matePairs.map(_.scaffold1).max, "Not all scaffold1 the same!")
  assert(scaffold2 == matePairs.map(_.scaffold2).max, "Not all scaffold2 the same!")

  override def toString = {
    "<LinkCluster %d %d:%d+%d %d:%d+%d %d>".format(nLinks, scaffold1, minCoord1, maxCoord1-minCoord1,
      scaffold2, minCoord2, maxCoord2-minCoord2, mq)
  }


  def dumpCoords(coords: Array[MatePair]) = coords foreach println
}

object Scaffold {

  def findClusters(coords: Array[MatePair], windowSize: Int, minCluster: Int = 10) = {
    val clustersByDistance = findClustersInternal(coords, windowSize, minCluster, {_.distance})
    val clusters = clustersByDistance map {findClustersInternal(_, windowSize, minCluster, {_.longCoord1})}

    val linkClusters = clusters.flatten.map({new LinkCluster(_)}).sortWith({_.nLinks > _.nLinks})
    linkClusters foreach println
    linkClusters
  }

  // This is used to find clusters within a neighborhood along some dimension (mpFunc)
  def findClustersInternal(unsortedCoords: Array[MatePair],
                           window: Int,
                           minCluster: Int = 10,
                           mpFunc: MatePair => Long) = {
    //println("finding clusters: size=%d window=%d min=%d func=%s".format(unsortedCoords.size,
    //  window, minCluster, mpFunc))

    val coords = unsortedCoords.sortWith({mpFunc(_) < mpFunc(_)})
    //dumpCoords(coords)
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
    println("analyzing strays in " + bam)
    val mm = bam.strayMateMap.mateMap
    val window = (bam.insertSizeSigma * 4).toInt
    val genomeSize = bam.getSeqs.map({_.getSequenceLength}).sum
    println("genome size " + genomeSize)
    println("max imsert " + bam.maxInsertSize)

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
    val backgroundRate = (nLinks.toFloat * window.toFloat / genomeSize.toFloat).toInt
    println("ambig=" + ambig + " same=" + intra + " innies=" + innie + " links=" + nLinks + " background=" + backgroundRate)
    val coords = links.toArray

    findClusters(links.toArray, window, /*1 * backgroundRate*/ 25)
  }


  def dumpCoords(coords: Array[MatePair]) = coords foreach println
}
