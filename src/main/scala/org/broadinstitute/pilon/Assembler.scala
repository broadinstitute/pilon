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


  def pathsForward(kmersIn: List[String], branches: Int = 0): List[List[String]] = {
    if (kGraph.isEmpty) buildGraph
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
          return pathsForward(next1 :: kmers, branches + 1) ++
            pathsForward(next2 :: kmers, branches + 1)
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


  def pathForward(kmersIn: List[String], target: String, forks: Int): List[String] = {
    var kmers = kmersIn
    def debug(s: String, kmerPath: List[String] = kmers) = {
      if (Pilon.debug) {
        val ks = kmerPathString(kmerPath)
        for (i <- 1 to forks) print("  ")
        println(s + "[(" + ks.length + ")" + ks + "]")
      }
    }
    while (true) {
      val kmer = kmers.head
      //if (target.indexOf(kmer) >= 99999999) {

      /*
      if (kmer == target) {
        debug("pFw:target ")
        return kmers
      }
      */
      if (!(pileups contains kmer)) {
	    // we're off the graph, so punt!
	    debug("pFw:off ")
        return kmers
      }
      val seen0 = kmers.tail count {_ == kmer}
      val pu = pileups(kmer)
      val bc = pu.baseCall
      val prefix = kmer.substring(1)
      if (seen0 > 1 || forks > 10) {
    	// chop off most recent kmer, as we don't want to go around again!
        val repeatStart = kmers.tail.indexOf(kmer)
    	if (Pilon.debug) {
    	  val repeatKmers = kmers.slice(0, repeatStart)
    	  val repeatStr = kmerPathString(repeatKmers)
    	  debug("pFw: loop " + seen0 + " " + forks +
    			  	" (" + repeatStr.length + ")" + repeatStr)          
    	}

    	return kmers.tail //tail drop repeatStart+1
      } else if (pu.depth < Assembler.minDepth) { // TODO: fixed depth or computed?
        // not enough depth to move forward
        debug("pFw: " + pu)
        return kmers
      } else if (bc.homo) {
        // no-brainer: only one path forward
        val newKmer = prefix + bc.base
        kmers ::= newKmer
      } else {
        // two choices forward: we choose based on where we've been before
        val newKmer1 = prefix + bc.base
        val newKmer2 = prefix + bc.altBase
        if (seen0 > 0) {
          // if we've been here before and we have two ways forward, let's see
          // what we've explored
          val seen1 = kmers contains newKmer1
          val seen2 = kmers contains newKmer2
          if (Pilon.debug) {
            val repeatStart = kmers.tail.indexOf(kmer)
            val repeatKmers = kmers.slice(0, repeatStart)
            val repeatStr = kmerPathString(repeatKmers)
            debug("pFw: loop " + seen0 + "," + seen1 + "," + seen2 + " " + forks +
               " (" + repeatStr.length + ")" + repeatStr)
          }
          if (seen1) {
            // we've already taken 1st branch, so try extending with 2nd 
            // if we haven't already
            if (seen2) {
       	      debug("pFw:s1+2 " + pu)
              return kmers
            } else
              kmers ::= newKmer2
          } else if (seen2) {
            // likewise, if we've been through 2nd, try 1st if we haven't
            kmers ::= newKmer1
          } else {
            // shouldn't happen; if we've seen this kmer, we should have moved forward
            assert(false, "shouldn't happen")
            return kmers
          }
        } else {
          // we haven't been here, but two ways forward: try both and take longest extension
          debug("pFw:forkA " + pu + " " + forks)
          val path1 = pathForward(newKmer1 :: kmers, target, forks + 1)
          debug("pFw:forkB")
          val path2 = pathForward(newKmer2 :: kmers, target, forks + 1)
          debug("A:", path1)
          debug("B:", path2)
          val hit1 = target.indexOf(path1.head) >= 0
          val hit2 = target.indexOf(path2.head) >= 0
          if (hit1 && !hit2) return path1 
          else if (hit2 && !hit1) return path2
          else
            if (path1.length >= path2.length) return path1
          else return path2
        }
      }
    }
    kmers
  }
  
  def pathForward(startingKmer: String, target: String = "") : String = {
    require(startingKmer.length == K, "starting kmer must be size K")
    val path = pathForward(List(startingKmer), target, 0)
    val bases = kmerPathString(path)
    if (Pilon.debug) println("pFw:" + bases.length + " " + bases)
    if (Pilon.debug) print("pathsForward:")
    val paths = pathsForward(List(startingKmer))
    val pathStrs = (paths map { kmerPathString(_) }) sortWith {(a, b) => a.length > b.length}
    if (Pilon.debug) {
      println(paths.length + " paths")
      pathStrs foreach {p => println("  [" + p.length + "]" + p)}
      assert(pathStrs.indexOf(bases) == 0)
    }
    bases
  }

  def pathForwardOld(anchor: String, kmersVisited: HashSet[String] = HashSet()): String = {
    val startingKmer = anchor.slice(anchor.length - K, anchor.length)
    if (pileups contains startingKmer) {
      if (kmersVisited contains startingKmer) {
        anchor
      } else {
        val pu = pileups(startingKmer)
        val bc = pu.baseCall
        //if (bc.homo && !bc.indel && pu.depth >= GapFiller.minDepth) {
        if (pu.depth >= Assembler.minDepth && // TODO: fixed depth or computed?
          (bc.homo || (Pilon.diploid && bc.majority))) {
          kmersVisited += startingKmer
          pathForwardOld(anchor + bc.base, kmersVisited)
        } else {
          anchor
        }
      }
    } else {
      anchor
    }
  }

  def pathReverse(anchor: String, target: String = ""): String = {
    val rcAnchor = Bases.reverseComplement(anchor)
    //val path = pathForward(rcAnchor, kmersVisited)
    val path = pathForward(rcAnchor)
    Bases.reverseComplement(path)
  }

  def tryForward(anchor: String, target: String = ""): String = {
    var longest = ""
    var longestFull = ""
    val anchorK = anchor.length / K
    for (offsetK <- 0 until (anchorK min 2)) {
      val offset = K * offsetK
      val start = anchor.slice(offset, offset + K)
      if (Pilon.debug) println("tryForward o=" + offset + " " + start)
      val path = pathForward(start, target)
      if (Pilon.debug) {
    	  val pathOld = pathForwardOld(start)
    	  if (pathOld.length != path.length)
    	    println("pFw: new " + path.length + " old " + pathOld.length)
      } 
      val fullPath = anchor.slice(0, offset) + path
      // If we make it sufficiently beyond our starting point, call it good
      if (fullPath.length > anchor.length + Assembler.minExtend)
        return fullPath
      if (path.length > longest.length) {
        longest = path
        longestFull = fullPath
      }
    }
    //anchor
    longestFull
  }

  def tryReverse(anchor: String, target: String = "") = {
    if (Pilon.debug) println("tryReverse")
    val rcAnchor = Bases.reverseComplement(anchor)
    val rcTarget = Bases.reverseComplement(target)
    val path = tryForward(rcAnchor, rcTarget)
    Bases.reverseComplement(path)
  }
  
  def bridge(left: String, right: String) = {
    val pathForward = tryForward(left, right.substring(right.length - K))
    if (Pilon.debug) println("bridgeF:(" + pathForward.length + ")" + pathForward)
    val pathReverse = tryReverse(right, left.substring(0, K))
    if (Pilon.debug) println("bridgeR:(" + pathReverse.length + ")" + pathReverse)
    (pathForward, pathReverse)
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
        val forward = pathForward(kmer)
        val reverse = pathReverse(kmer)
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
