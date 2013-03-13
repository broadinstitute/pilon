/*
 * Copyright (c) 2013. The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc.
 * All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package org.broadinstitute.pilon

import collection.mutable.ArrayBuffer
import math.pow

// Calculates the moments of a sample distribution.
class NormalDistribution(values: Array[Double], nMoments: Int) {
  require(values.size > 0, "can't compute moments of empty distribution")
  require(nMoments > 1, "assumes at least two moments (mean, stddev)")

  val mean = values.sum / values.size
  val median = {
    val sorted = values.clone.sorted
    val n = sorted.size
    if (n % 2 == 0)
      (values(n/2) + values(n/2 - 1)) / 2.0
    else
      values(n/2)
  } 
  
  val moments = new Array[Double](nMoments)
  moments(0) = mean
  for (n <- 1 until nMoments) {
	for (v <- values) moments(n) += pow(v - mean, n + 1)
    moments(n) = pow(moments(n) / values.size, 1.0 / (n + 1))
  }

  def toSigma(value: Double, maxSigma: Double = 5.0) = {
    val sigma = (value - moments(0)) / moments(1)
    sigma
  }
  
  def toSigma10x(value: Double) = (toSigma(value) * 10.0).round.toInt
  
  def fromSigma(sigma: Double) = moments(0) + sigma * moments(1)
    
  def this(ivalues: Array[Int], nMoments: Int) =
   this(for { v<-ivalues } yield v.toDouble, nMoments) 
  
  def this(ivalues: Array[Short], nMoments: Int) =
   this(for { v<-ivalues } yield v.toDouble, nMoments) 
  
  def this(ivalues: Array[Byte], nMoments: Int) =
   this(for { v<-ivalues } yield v.toDouble, nMoments) 
  
  override def toString = "<moments: n=" + values.size + ",moments=" + (moments mkString ",") + ">"
}