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
import scala.annotation.tailrec
import collection.JavaConversions._
import collection.mutable.{ Map, HashMap, Set, HashSet }
import net.sf.samtools._

object Assembler {
  val K = 47
  val minDepth = 5
  val minGap = 10
  val minExtend = 20
  val minNovel = 200
  val maxBranches = 5
  type Kmer = String
  type KmerPileup = HashMap[Kmer, PileUp]
  type KmerGraph = HashMap[Kmer, (Kmer,Int)]
}

class Assembler(val minDepth: Int = Assembler.minDepth) {
  import Assembler._
  var pileups: KmerPileup = HashMap()
  val kGraph: KmerGraph = HashMap()
  val altGraph: KmerGraph = HashMap()
  var nReads: Long = 0
  var nBases: Long = 0

  def addReads(reads: List[SAMRecord]) = {
    reads foreach addRead
  }

  def addPair(r1: SAMRecord, r2: SAMRecord) = {
    addRead(r1)
    addRead(r2)
  }

  def addRead(r: SAMRecord) = {
    val bases = r.getReadString
    if (bases.size > K) {
      val quals = r.getBaseQualities
      val mq = r.getMappingQuality
      addToPileups(bases, quals, mq)

      val rcBases = Bases.reverseComplement(bases)
      val rcQuals = quals.reverse
      addToPileups(rcBases, rcQuals, mq)

      nReads += 1
      nBases += bases.length
      if (Pilon.verbose && nReads % 10000 == 0) print("..." + nReads)
      if (nReads % 100000 == 0) prunePileups(2)
    }
  }

  def addToPileups(bases: String, quals: Array[Byte], mq: Int) = {
    val length = bases.length
    for (offset <- 0 to length - K - 1) {
      val kmer = bases.slice(offset, offset + K)
      val qmer = quals.slice(offset, offset + K)
      if (!(pileups contains kmer))
        pileups(kmer) = new PileUp
      pileups(kmer).add(bases(offset + K), quals(offset + K), mq)
    }
  }

  def buildGraph = {
    val targets = Set[Kmer]()
    val multiTargets = Set[Kmer]()

    def addLink(g: KmerGraph, k1: Kmer, k2: Kmer, weight: Int) = {
      g(k1) = (k2, weight)
      if (targets contains k2) multiTargets += k2
      else targets += k2
    }


    if (Pilon.debug) println("building kmer Graph")

    for ((k, pu) <- pileups.iterator) {
      if (pu.depth >= Assembler.minDepth) {
        val bc = pu.baseCall
        val prefix = k.substring(1)
        addLink(kGraph, k, prefix + bc.base, pu.baseCount.sums(bc.baseIndex).toInt)
        if (!bc.homo)
          addLink(altGraph, k, prefix + bc.altBase, pu.baseCount.sums(bc.altBaseIndex).toInt)
      }
    }
    // We don't need  or targets any more, free up the memory
    //pileups = null

    if (Pilon.debug) println("kmer graph: t=" + targets.size + " mt=" + multiTargets.size)

    val sources = kGraph.keys filter {k => (multiTargets contains k) || !(targets contains k)}
    //if (Pilon.debug) println("kmer Graph built")
    //if (Pilon.debug) println("sources: " + sources.size, " " + sources)
    for (s <- sources) {

    }
  }


  def prunePileups(minCount: Int = Assembler.minDepth) = {
    if (Pilon.debug) print("[prune " + pileups.size)
    for ((k, pu) <- pileups.iterator)
      if (pu.count < minCount) pileups -= k
    if (Pilon.debug) print("->" + pileups.size + "]")
  }
  
  def kmerPathString(kmers: List[String], prependLength: Boolean = false) = {
    val path = kmers.reverse
    val pathStr = path.head + (path.tail map {_.last}).mkString("")
    if (prependLength)
      "(" + pathStr.length + ")" + pathStr
    else
      pathStr
  }


  def kmerPathsForward(kmersIn: List[String], branches: Int = 0): List[List[String]] = {
    var kmers = kmersIn
    while (true) {
      val kmer = kmers.head
      if (!(kGraph contains kmer)) {
        // end of the line
        return List(kmers)
      } else if (altGraph contains kmer) {
        // two choices forward
        val next1 = kGraph(kmer)._1
        val next2 = altGraph(kmer)._1
        val seen1 = kmers.tail contains next1
        val seen2 = kmers.tail contains next2

        // explore branches not already taken
        if (seen1 && seen2) return List(kmers)
        else if (seen1 && !seen2) kmers ::= next2
        else if (seen2 && !seen1) kmers ::= next1
        else {
          //if (Pilon.debug) println("branch " + branches + " " + next1 + " " + next2)
          if (branches < maxBranches)
            return kmerPathsForward(next1 :: kmers, branches + 1) ++
              kmerPathsForward(next2 :: kmers, branches + 1)
          else return List(kmers)
        }
      } else {
        // only one way forward
        val next = kGraph(kmer)._1
        val nextCount = kmers count {_ == next}
        // punt if we're looping (we'll allow one full time around repeat)
        if (nextCount > 1)
          return List(kmers)
        kmers ::= next
      }
    }
    // shouldn't get here
    Nil
  }

  def pathsForward(startingKmer: String): List[String] = {
    require(startingKmer.length == K, "starting kmer must be size K")
    if (kGraph.isEmpty) buildGraph
    if (Pilon.debug) print("pathsForward: " + startingKmer)
    val kmerPaths = kmerPathsForward(List(startingKmer))
    val paths = (kmerPaths map { kmerPathString(_) }) sortWith {(a, b) => a.length > b.length}
    if (Pilon.debug) {
      println(": " + paths.length + " paths")
      paths foreach {p => println("  [" + p.length + "]" + p)}
    }
    paths
  }

  def pathsReverse(startingKmer: String): List[String] = {
    val paths = pathsForward(Bases.reverseComplement(startingKmer))
    paths map {Bases.reverseComplement(_)}
  }

  def tryForward(anchor: String): List[String] = {
    if (Pilon.debug) println("tryForward: [" + anchor.length + "]" + anchor)
    if (anchor.length < K) List(anchor)
    else {
      val startingKmer = anchor.substring(0, K)
      val paths = pathsForward(startingKmer) filter {_.length > anchor.length + minExtend}
      if (paths.length > 0)
        paths
      else
        tryForward(anchor.substring(K)) map {p => startingKmer + p}
    }
  }

  def tryReverse(anchor: String) = {
    if (Pilon.debug) println("tryReverse")
    val rcAnchor = Bases.reverseComplement(anchor)
    val paths = tryForward(rcAnchor) map { Bases.reverseComplement(_) }
    if (Pilon.debug) {
      println("RC " + paths.length + " paths")
      paths foreach {p => println("  [" + p.length + "]" + p)}
    }
    paths
  }

  def multiBridge(left: String, right: String) = {
    val pathsForward = tryForward(left)
    val pathsReverse = tryReverse(right)
    (pathsForward, pathsReverse)
  }



  def novel: List[String] = {
    println("Assembling novel sequence: " + this)
    prunePileups()
    var usedKmers = HashSet[String]()
    var paths: List[String] = List()
    var n = 0

    for (kmer <- pileups.keysIterator) {
      if (!(usedKmers contains kmer)) {
        var dup = false
        val forwards = pathsForward(kmer)
        val reverses = pathsReverse(kmer)
        val forward = if (forwards.isEmpty) ""  else forwards.head
        val reverse = if (reverses.isEmpty) ""  else reverses.head
        val path = reverse + forward.substring(K)
        path.sliding(K) foreach { k =>
          if (usedKmers contains kmer)
            dup = true
          usedKmers += k
          usedKmers += Bases.reverseComplement(k)
        }
        if (!dup && path.length > Assembler.minNovel) {
          paths ::= path
          println("novel " + path.length + " " + path)
        }
      }
      n += 1
      //if (n % 100000 == 0) print("..." + n)
    }
    //println
    paths
  }

  override def toString =
    "<assembler K=" + K + " nReads=" + nReads + " nBases=" +
      nBases + " nKmers=" + pileups.size + ">"
}
