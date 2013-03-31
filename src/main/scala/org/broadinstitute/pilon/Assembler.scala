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

import collection.JavaConversions._
import collection.mutable.{ Map, HashMap, Set, HashSet }
import net.sf.samtools._

object Assembler {
  val K = 47
  val minDepth = 5
  val minGap = 10
  val minExtend = 20
  val minNovel = 200
  type Kmer = String
  type KmerGraph = HashMap[Kmer, PileUp]
}

class Assembler(val minDepth: Int = Assembler.minDepth) {
  val K = Assembler.K
  val graph: Assembler.KmerGraph = HashMap()
  var nReads: Long = 0
  var nBases: Long = 0

  def addReads(reads: List[SAMRecord]) = {
    reads foreach addRead
  }

  def addRead(r: SAMRecord) = {
    val bases = r.getReadString
    if (bases.size > K) {
      val quals = r.getBaseQualities
      val mq = r.getMappingQuality
      addToGraph(bases, quals, mq)

      val rcBases = Bases.reverseComplement(bases)
      val rcQuals = quals.reverse
      addToGraph(rcBases, rcQuals, mq)

      nReads += 1
      nBases += bases.length
      if (Pilon.verbose && nReads % 10000 == 0) print("..." + nReads)
      if (nReads % 100000 == 0) pruneGraph(2)
    }
  }

  def addToGraph(bases: String, quals: Array[Byte], mq: Int) = {
    val length = bases.length
    for (offset <- 0 to length - K - 1) {
      val kmer = bases.slice(offset, offset + K)
      val qmer = quals.slice(offset, offset + K)
      if (!(graph contains kmer))
        graph(kmer) = new PileUp
      graph(kmer).add(bases(offset + K), quals(offset + K), mq)
    }
  }

  def pruneGraph(minCount: Int = Assembler.minDepth) = {
    if (Pilon.debug) print("[prune " + graph.size)
    for ((k, pu) <- graph.iterator)
      if (pu.count < minCount) graph -= k
    if (Pilon.debug) print("->" + graph.size + "]")
  }
  
  def kmerPathString(kmers: List[String]) = {
    val path = kmers.reverse
    path.head + (path.tail map {_.last}).mkString("")
  }

  def pathForwardRecurse(kmers: List[String], target: String, forks: Int): List[String] = {
    val kmer = kmers.head
    if (target.indexOf(kmer) >= 99999999) {
      kmers
    } else if (graph contains kmer) {
      val seen0 = kmers.tail count {_ == kmer}
      val pu = graph(kmer)
      val bc = pu.baseCall
      val prefix = kmer.substring(1)
      //if (bc.homo && !bc.indel && pu.depth >= GapFiller.minDepth) {
      //if (Pilon.debug) println("pFw:" + kmer + " " + pu) 
      if (seen0 > 1 || forks > 10) {
    	  if (Pilon.debug) println("pFw:twice " + pu + " " + forks) 
          kmers
      } else if (pu.depth < Assembler.minDepth) { // TODO: fixed depth or computed?
        // not enough depth to move forward
    	if (Pilon.debug) println("pFw: " + pu) 
        kmers
      } else if (bc.homo) {
        val newKmer = prefix + bc.base
        pathForward(newKmer :: kmers, target, forks)
      } else {
        val newKmer1 = prefix + bc.base
        val newKmer2 = prefix + bc.altBase
        if (seen0 > 0) {
          // if we've been here before and we have two ways forward, let's see
          // what we've explored
          val seen1 = kmers contains newKmer1
          val seen2 = kmers contains newKmer2
          if (Pilon.debug) println("pFw: " + seen0 + " " + seen1 + " " + seen2 + " " + forks)
          if (seen1) {
            // we've already taken 1st branch, so try extending with 2nd 
            // if we haven't already
            if (seen2) {
       	      if (Pilon.debug) println("pFw:s1+2 " + pu) 
              kmers
            } else 
              pathForward(newKmer2 :: kmers, target, forks)
          } else if (seen2) {
            // likewise, if we've been through 2nd, try 1st if we haven't
            pathForward(newKmer1 :: kmers, target, forks)
          } else {
            // shouldn't happen; if we've seen this kmer, we should have moved forward
            assert(false, "shouldn't happen")
            kmers
          }
        } else {
          // we haven't been here, but two ways forward: try both and take longest extention
       	  if (Pilon.debug) println("pFw:fork " + pu + " " + forks) 
          val path1 = pathForward(newKmer1 :: kmers, target, forks + 1)
          val path2 = pathForward(newKmer2 :: kmers, target, forks + 1)
          val hit1 = target.indexOf(path1.head) >= 0
          val hit2 = target.indexOf(path2.head) >= 0
          if (hit1 && !hit2) path1 
          else if (hit2 && !hit1) path2
          else if (path1.length >= path2.length) path1
          else path2
        }
      }
    } else {
      if (Pilon.debug) println("pFw:off ")
      kmers	// we're off the graph, so punt! 
    }
  }
  
  def pathForward(kmersIn: List[String], target: String, forks: Int): List[String] = {
    var kmers = kmersIn
    while (true) {
      val kmer = kmers.head
      if (target.indexOf(kmer) >= 99999999) {
        if (Pilon.debug) println("pFw:target ")
        return kmers
      }
      if (!(graph contains kmer)) {
        if (Pilon.debug) println("pFw:off ")
        return kmers	// we're off the graph, so punt!
      }
      val seen0 = kmers.tail count {_ == kmer}
      val pu = graph(kmer)
      val bc = pu.baseCall
      val prefix = kmer.substring(1)
      //if (bc.homo && !bc.indel && pu.depth >= GapFiller.minDepth) {
      //if (Pilon.debug) println("pFw:" + kmer + " " + pu) 
      if (seen0 > 1 || forks > 10) {
    	  if (Pilon.debug) println("pFw:twice " + pu + " " + forks) 
          return kmers
      } else if (pu.depth < Assembler.minDepth) { // TODO: fixed depth or computed?
        // not enough depth to move forward
    	if (Pilon.debug) println("pFw: " + pu) 
        return kmers
      } else if (bc.homo) {
        val newKmer = prefix + bc.base
        kmers ::= newKmer
        //pathForward(newKmer :: kmers, target, forks)
      } else {
        val newKmer1 = prefix + bc.base
        val newKmer2 = prefix + bc.altBase
        if (seen0 > 0) {
          // if we've been here before and we have two ways forward, let's see
          // what we've explored
          val seen1 = kmers contains newKmer1
          val seen2 = kmers contains newKmer2
          if (Pilon.debug) println("pFw: " + seen0 + " " + seen1 + " " + seen2 + " " + forks)
          if (seen1) {
            // we've already taken 1st branch, so try extending with 2nd 
            // if we haven't already
            if (seen2) {
       	      if (Pilon.debug) println("pFw:s1+2 " + pu) 
              return kmers
            } else 
              kmers ::= newKmer2
              //pathForward(newKmer2 :: kmers, target, forks)
          } else if (seen2) {
            // likewise, if we've been through 2nd, try 1st if we haven't
            kmers ::= newKmer1
            //pathForward(newKmer1 :: kmers, target, forks)
          } else {
            // shouldn't happen; if we've seen this kmer, we should have moved forward
            assert(false, "shouldn't happen")
            return kmers
          }
        } else {
          // we haven't been here, but two ways forward: try both and take longest extension
       	  if (Pilon.debug) println("pFw:fork " + pu + " " + forks) 
          val path1 = pathForward(newKmer1 :: kmers, target, forks + 1)
          val path2 = pathForward(newKmer2 :: kmers, target, forks + 1)
          val hit1 = target.indexOf(path1.head) >= 0
          val hit2 = target.indexOf(path2.head) >= 0
          if (hit1 && !hit2) return path1 
          else if (hit2 && !hit1) return path2
          else if (path1.length >= path2.length) return path1
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
    bases
  }

  def pathForwardOld(anchor: String, kmersVisited: HashSet[String] = HashSet()): String = {
    val startingKmer = anchor.slice(anchor.length - K, anchor.length)
    if (graph contains startingKmer) {
      if (kmersVisited contains startingKmer) {
        anchor
      } else {
        val pu = graph(startingKmer)
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
      if (Pilon.debug) print("tryForward o=" + offset + " ")
      val start = anchor.slice(offset, offset + K)
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
    if (Pilon.debug) print("tryReverse ")
    val rcAnchor = Bases.reverseComplement(anchor)
    val rcTarget = Bases.reverseComplement(target)
    val path = tryForward(rcAnchor, rcTarget)
    Bases.reverseComplement(path)
  }
  
  def bridge(left: String, right: String) = {
    val pathForward = tryForward(left, right.substring(right.length - K))
    val pathReverse = tryReverse(right, left.substring(0, K))
    (pathForward, pathReverse)
  }

  def novel: List[String] = {
    println("Assembling novel sequence: " + this)
    pruneGraph()
    var usedKmers = HashSet[String]()
    var paths: List[String] = List()
    var n = 0

    for (kmer <- graph.keysIterator) {
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
    "<assembler K=" + K + " nReads=" + nReads + " nBases=" + nBases + " nKmers=" + graph.size + ">"
}
