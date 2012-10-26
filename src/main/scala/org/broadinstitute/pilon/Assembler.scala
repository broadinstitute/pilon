package org.broadinstitute.pilon

import collection.JavaConversions._
import collection.mutable.{Map, HashMap, Set, HashSet}
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

class Assembler {
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
        graph(kmer).add(bases(offset+K), quals(offset+K), mq)
    }
  }
  
  def pruneGraph(minCount: Int = Assembler.minDepth) = {
    print("[prune " + graph.size)
    for ((k, pu) <- graph.iterator)
      if (pu.count < minCount) graph -= k
    print("->" + graph.size + "]")
  }

  def pathForward(anchor: String, kmersVisited: HashSet[String] = HashSet()): String = {
    val startingKmer = anchor.slice(anchor.length - K, anchor.length)
    if (graph contains startingKmer) {
      if (kmersVisited contains startingKmer) {
        if (Pilon.debug) println(anchor.length + ":" + anchor + " <loop>")
        anchor
      } else {
        val pu = graph(startingKmer)
        val bc = pu.baseCall
        if (bc.homo && !bc.indel && pu.depth >= GapFiller.minDepth) {
          kmersVisited += startingKmer
          pathForward(anchor + bc.base, kmersVisited)
        } else {
          if (Pilon.debug) println(anchor.length + ":" + anchor + " " + pu)
          anchor
        }
      }
    } else {
      if (Pilon.debug) println(anchor.length + ":" + anchor + " <off graph>")
      anchor
    }
  }

  def pathReverse(anchor: String, kmersVisited: HashSet[String] = HashSet()): String = {
    val rcAnchor = Bases.reverseComplement(anchor)
    val path = pathForward(rcAnchor, kmersVisited)
    Bases.reverseComplement(path)
  }
  
  def tryForward(anchor: String): String = {
    var longest = ""
    var longestFull = ""
    for (offset <- 0 until (anchor.length, K)) {
      if (Pilon.debug) print("tryForward o=" + offset + " ")
      val start = anchor.slice(offset, offset+K)
      val path = pathForward(start)
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
  
  def novel : List[String] = {
    println("Assembling novel sequence: " + this)
    pruneGraph()
    var usedKmers = HashSet[String]()
    var paths : List[String] = List()
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
