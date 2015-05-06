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
import scala.annotation.tailrec
import collection.JavaConversions._
import collection.mutable.{ Map, HashMap, Set, HashSet }
import htsjdk.samtools._

object Assembler {
  var K = 47
  val minDepth = 5
  val minExtend = 20
  val maxBranches = 5
  val minNovel = 200
  val minNovelPct = 50
  type Kmer = String
  type KmerPileup = HashMap[Kmer, PileUp]
  type KmerGraph = HashMap[Kmer, Kmer]
}

class Assembler(val minDepth: Int = Assembler.minDepth) {
  import Assembler._
  var pileups: KmerPileup = HashMap()
  val kGraph: KmerGraph = HashMap()
  val altGraph: KmerGraph = HashMap()
  var nReads: Long = 0
  var nBases: Long = 0
  var loopLength = 0
  var loopSequence = ""

  def addReads(reads: List[SAMRecord]) = {
    reads foreach addRead
  }

  def addPair(r1: SAMRecord, r2: SAMRecord) = {
    addRead(r1)
    addRead(r2)
  }

  def addRead(r: SAMRecord) = {
    val bases = r.getReadString
    val length = bases.size
    if (length > K) {
      val quals = {
        val q = r.getBaseQualities
        if (q.length > 0) q
        else Array.fill[Byte](length)(Pilon.defaultQual)
      }
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
      if (!(pileups contains kmer))
        pileups(kmer) = new PileUp
      //if (Pilon.debug) println("K=" + K + " " + kmer + " offset=" + offset + " len=" + quals.length + " " + quals)
      pileups(kmer).add(bases(offset + K), quals(offset + K), mq)
    }
  }

  def dumpReads(prefix: String, reads: List[SAMRecord]): Unit = {

  }

  // Used to create an assembly graph from sequence, e.g., contigs.
  // Uses fake base and mapping qualities.
  def addSeq(bases: String) = {
    val mq = 10
    val quals = Array.fill(bases.length) { 10.toByte }
    addToPileups(bases, quals, mq)
    val rcBases = Bases.reverseComplement(bases)
    addToPileups(rcBases, quals, mq)
  }

  // Used to create an assembly graph directly from sequence, e.g., contigs.
  def addGraphSeq(bases: String) = {
    graphSeq(bases)
    val rcBases = Bases.reverseComplement(bases)
    graphSeq(rcBases)
  }
  
  def addGraphSeqs(seqs: List[String]) = seqs foreach addGraphSeq

  // Used to create an assembly graph directly from sequence, 
  // e.g., contigs.
  def graphSeq(bases: String) = {
    val length = bases.length
    for (offset <- 0 to length - K - 1) {
      val k = bases.slice(offset, offset + K)
      val nextK = k.substring(1) + bases(offset + K)

      if ((kGraph contains k) && (kGraph(k) != nextK))
    	addLink(altGraph, k, nextK, 1)
      else
    	addLink(kGraph, k, nextK, 1)
    }
  }
  
  def addLink(g: KmerGraph, k1: Kmer, k2: Kmer, weight: Int) = {
	// not using weight for now
	g(k1) = k2
  }
  
  def buildGraph = {

    if (Pilon.debug) println("building kmer Graph")

    for ((k, pu) <- pileups.iterator) {
      if (pu.depth >= minDepth) {
        val bc = pu.baseCall
        val prefix = k.substring(1)
        val nextK = prefix + bc.base
        val weight = pu.baseCount.sums(bc.baseIndex).toInt
        if ((kGraph contains k) && (kGraph(k) != nextK))
          addLink(altGraph, k, nextK, weight)
        else
          addLink(kGraph, k, nextK, weight)
        if (!bc.homo)
          addLink(altGraph, k, prefix + bc.altBase, pu.baseCount.sums(bc.altBaseIndex).toInt)
      }
    }
    // We don't need pileups any more, free up the memory
    pileups = HashMap()

    if (Pilon.debug) println("kmer graph: t=" + kGraph.size + " mt=" + altGraph.size)
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

  def noteKmerLoop(loopIndex: Int, kmers: List[String]) = {
    val length = loopIndex + 1
    if (loopLength == 0 || length < loopLength) {
      loopLength = length
      loopSequence = kmerPathString(kmers.take(loopIndex+1)).take(length)
      if (Pilon.verbose)
        println("# loop " + loopLength + ": " + loopSequence)
    }
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
        val next1 = kGraph(kmer)
        val next2 = altGraph(kmer)
        val seen1 = kmers.tail contains next1
        val seen2 = kmers.tail contains next2

        if (seen1 || seen2) {
          val loop = kmers.indexOf(next1) max kmers.indexOf(next2)
          noteKmerLoop(loop, kmers)
          // if we're being careful about tandem repeats, better punt now!
          if (Pilon.trSafe) return List(kmers)
        }

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
        val next = kGraph(kmer)
        val loop = kmers indexOf next
        if (loop >= 0) {
          noteKmerLoop(loop, kmers)
          // if we're being careful about tandem repeats, better punt now!
          if (Pilon.trSafe) return List(kmers)

          // punt if we're looping (we'll allow one full time around repeat)
          val nextCount = kmers count {_ == next}
          if (nextCount > 1)
            return List(kmers)
        }
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
      //val paths = pathsForward(startingKmer) filter {_.length > anchor.length + minExtend}
      //if (paths.length > 0)
      val paths = pathsForward(startingKmer)
      if (paths.length > 0 && paths(0).length > anchor.length + minExtend)
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
    (pathsForward, pathsReverse, loopSequence)
  }

  def novel(ref: Assembler): List[String] = {
    def novelKmers(seq: String) = seq.sliding(K) count { !ref.kGraph.contains(_)  }

    if (Pilon.verbose) println("Assembling novel sequence: " + this + ", ref " + ref)
    prunePileups()
    var usedKmers = HashSet[String]()
    var paths: List[String] = List()
    var n = 0

    for (kmer <- pileups.keysIterator) {
      if (!(usedKmers contains kmer)) {
        val forwards = pathsForward(kmer)
        val reverses = pathsReverse(kmer)
        val forward = if (forwards.isEmpty) ""  else forwards.head
        val reverse = if (reverses.isEmpty) ""  else reverses.head
        val path = reverse + forward.substring(K)
        path.sliding(K) foreach { k =>
          usedKmers += k
          usedKmers += Bases.reverseComplement(k)
        }
        if (path.length >= minNovel)
          paths ::= path
      }
      n += 1
      //if (n % 100000 == 0) print("..." + n)
    }
    paths = paths sortWith {(a, b) => a.length > b.length}
    paths filter { path =>
      val kLength = path.length - (K - 1)
      val novel = novelKmers(path)
      val novelPct = Utils.pct(novel, kLength)
      //println("L=" + path.length + " K=" + kLength + " N=" + novel + " P=" + novelPct)
      if (novel >= minNovel && novelPct >= minNovelPct) {
        // side effect...if we're accepting this one, add it to ref to avoid future dupes
        ref.addGraphSeq(path)
        //ref.addSeq(path)
        //ref.buildGraph
        if (Pilon.verbose) println("novel " + path.length + " " + novelPct + "% " + path)
        true
      } else false
    }
  }

  override def toString =
    "<assembler K=" + K + " nReads=" + nReads + " nBases=" +
      nBases + " nKmers=" + pileups.size + ">"
}
