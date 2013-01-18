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

  def pathForward(kmers: List[String]): List[String] = {
    val kmer = kmers.head
    if (graph contains kmer) {
      val seen0 = kmers.tail count {_ == kmer}
      val pu = graph(kmer)
      val bc = pu.baseCall
      val prefix = kmer.substring(1)
      //if (bc.homo && !bc.indel && pu.depth >= GapFiller.minDepth) {
      //if (Pilon.debug) println("pFw:" + kmer + " " + pu) 
      if (seen0 > 1) {
    	  if (Pilon.debug) println("pFw:twice " + pu) 
          kmers
      } else if (pu.depth < Assembler.minDepth) { // TODO: fixed depth or computed?
        // not enough depth to move forward
    	if (Pilon.debug) println("pFw: " + pu) 
        kmers
      } else if (bc.homo) {
        val newKmer = prefix + bc.base
        pathForward(newKmer :: kmers)
      } else {
        val newKmer1 = prefix + bc.base
        val newKmer2 = prefix + bc.altBase
        if (seen0 > 0) {
          // if we've been here before and we have two ways forward, let's see
          // what we've explored
          val seen1 = kmers contains newKmer1
          val seen2 = kmers contains newKmer2
          if (Pilon.debug) println("pFw: " + seen0 + " " + seen1 + " " + seen2)
          if (seen1) {
            // we've already taken 1st branch, so try extending with 2nd 
            // if we haven't already
            if (seen2) {
       	      if (Pilon.debug) println("pFw:s1+2 " + pu) 
              kmers
            } else 
              pathForward(newKmer2 :: kmers)
          } else if (seen2) {
            // likewise, if we've been through 2nd, try 1st if we haven't
            pathForward(newKmer1 :: kmers)
          } else {
            // shouldn't happen; if we've seen this kmer, we should have moved forward
            assert(false, "shouldn't happen")
            kmers
          }
        } else {
          // we haven't been here, but two ways forward: try both and take longest extention
       	  if (Pilon.debug) println("pFw:fork " + pu) 
          val path1 = pathForward(newKmer1 :: kmers)
          val path2 = pathForward(newKmer2 :: kmers)
          if (path1.length >= path2.length) path1 else path2
        }
      }
    } else {
      if (Pilon.debug) println("pFw:off ")
      kmers	// we're off the graph, so punt! 
    }
  }
  
  def pathForward(startingKmer: String) : String = {
    require(startingKmer.length == K, "starting kmer must be size K")
    val path = pathForward(List(startingKmer))
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

  def pathReverse(anchor: String, kmersVisited: HashSet[String] = HashSet()): String = {
    val rcAnchor = Bases.reverseComplement(anchor)
    //val path = pathForward(rcAnchor, kmersVisited)
    val path = pathForward(rcAnchor)
    Bases.reverseComplement(path)
  }

  def tryForward(anchor: String): String = {
    var longest = ""
    var longestFull = ""
    for (offset <- 0 until (anchor.length, K)) {
      if (Pilon.debug) print("tryForward o=" + offset + " ")
      val start = anchor.slice(offset, offset + K)
      val path = pathForward(start)
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

  def tryReverse(anchor: String) = {
    if (Pilon.debug) print("tryReverse ")
    val rcAnchor = Bases.reverseComplement(anchor)
    val path = tryForward(rcAnchor)
    Bases.reverseComplement(path)
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
