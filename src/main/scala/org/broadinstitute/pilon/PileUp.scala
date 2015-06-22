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

import scala.collection.mutable.Map
import Utils._

class PileUp {
  var baseCount = new BaseSum
  var qualSum = new BaseSum
  //var mqSum = new BaseSum
  //var qmqSum = new BaseSum
  var mqSum = 0
  var qSum = 0
  var physCov = 0
  var insertSize = 0
  var badPair = 0
  var deletions = 0
  var delQual = 0
  var insertions = 0
  var insQual = 0
  var clips = 0
  var insertionList = List[Array[Byte]]()
  var deletionList = List[Array[Byte]]()

  def count = baseCount.sum
  def depth = baseCount.sum + deletions
    
  def baseIndex(c: Char) : Int = c match {
    case 'A' => 0
    case 'C' => 1
    case 'G' => 2
    case 'T' => 3
    case _ => -1
  }
  
  def indexBase(i: Int) : Char = "ACGT"(i)
  
  def weightedMq = {
    roundDiv(qualSum.sum, qSum)
  }

  def weightedQual = {
    roundDiv(qualSum.sum, mqSum)
  }
  
  def meanQual = {
    // scale mqSum by count/depth because mqSum includes deletions
    roundDiv(qualSum.sum, roundDiv(mqSum * count, depth))
  }

  // In all the mq handling, note that we use mq+1 internally to avoid ignoring mq0
  def meanMq = {
    roundDiv(mqSum - depth, depth)
  }
 
  // add a base to the stack
  def add(base: Char, qual: Int, mq: Int) = {
    val bi = baseIndex(base)
    if (bi >= 0 && qual >= Pilon.minQual) {
    	val mq1 = mq + 1
    	baseCount.add(bi)
    	qualSum.add(bi, qual * mq1)
    	mqSum += mq1
    	qSum += qual
    }
  }

  // remove a base from the stack
  def remove(base: Char, qual: Int, mq: Int) = {
    val bi = baseIndex(base)
    if (bi >= 0) {
    	val mq1 = mq + 1
    	baseCount.remove(bi)
    	qualSum.remove(bi, qual * mq1)
    	mqSum -= mq1
    	qSum -= qual
    }
  }
  
  def addInsertion(insertion: Array[Byte], qual: Int, mq: Int) = {
	  val mq1 = mq + 1
    insQual += mq1
    //mqSum += mq1 //don't add twice...already have a base here!
    qSum += qual
    insertionList ::= insertion
    insertions += 1
  }
  
  def addDeletion(deletion: Array[Byte], qual: Int, mq: Int) = {
	  val mq1 = mq + 1
    mqSum += mq1
	  delQual += mq1
    qSum += qual
	  deletionList ::= deletion
    deletions += 1
  }
  
  def totalQSum = qualSum.sum + insQual + delQual
  
  //def insPct = pct(insertions, count)
  //def delPct = pct(deletions, count + deletions)
  def insPct = pct(insQual, mqSum)
  def delPct = pct(delQual, mqSum)
  def indelPct = insPct + delPct
  def qualPct = {
    val (base, max, sum) = qualSum.maxBase
    pct(max, sum)
  }  
  def clipPct = pct(clips, depth)
  
  class BaseCall {
    val n = count
    val qs = qualSum
    val total = qs.sum //+ insQual + delQual
    val order = qs.order
    val baseIndex = order(0)
    val base = if (n > 0) indexBase(baseIndex) else 'N'
    val baseSum = qs.sums(baseIndex)
    val altBaseIndex = order(1)
    val altBase = indexBase(altBaseIndex)
    val altBaseSum = qs.sums(altBaseIndex)
    val homoScore = baseSum - (total - baseSum)
    val halfTotal = total / 2
    val heteroScore = total - (halfTotal - baseSum).abs - (halfTotal - altBaseSum).abs
    var homo = homoScore >= heteroScore
    val score = if (mqSum > 0) (homoScore - heteroScore).abs  * n / mqSum else 0
    val insertion = insertCall != ""
    val deletion = !insertion && deletionCall != ""
    val indel = insertion || deletion
    val called = (base != 'N' || indel)
    val q = if (n > 0) score / n else 0
    val highConfidence = q >= 10
    def callString(indelOk : Boolean = true) = {
      if (indelOk && insertion) insertCall
      else if (indelOk && deletion) deletionCall
      else base.toString //+ (if (!homo) "/" + altBase else "")
    }
    def baseMatch(refBase: Char) {
      refBase == base	// TODO: handle IUPAC codes
    }
    override def toString = {
      if (insertion) "bc=i" + insertCall + ",cq=" + insQual / insertions
      else if (deletion) "bc=d" + deletionCall + ",cq=" + delQual / deletions
      else 
    	  "bc=" + base + (if (!homo) "/" + altBase else "") +",cq=" + q
    }
  }
  
  def baseCall = new BaseCall
  
  def insertCall = indelCall(insertionList, insPct)
  def deletionCall = indelCall(deletionList, delPct)
  //def insertCall = indelCall(insertionList)
  //def deletionCall = indelCall(deletionList)
  
  def indelCall(indelList: List[Array[Byte]], pct: Int): String = {
    val map = Map.empty[String, Int]
    if (depth < Pilon.minMinDepth || pct < 33 || indelList.isEmpty) return ""
    for (indel <- indelList) {
      val indelStr = indel.toSeq map {_.toChar} mkString  ""
      map(indelStr) = map.getOrElse(indelStr, 0) + 1
    }
    val winner = map.toSeq.sortBy({ _._2 }).last
    if (winner._2 < indelList.length / 2) return ""
    val winStr = winner._1
    if (pct >= 50 - winStr.length)
      winStr
    else
        ""
  }


  override def toString = {
    "<PileUp " + (new BaseCall).toString + ",b=" + baseCount + "/" + qualSum.toStringPct + ",c=" + depth + "/" + (depth + badPair) + 
   	",i=" + insertions + ",d=" + deletions + ",q=" + weightedQual + ",mq=" + weightedMq +
    ",p=" + physCov + ",s=" + insertSize + ",x=" + clips + ">"
  }
  

}

