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

object Bases {
  // values chosen such that x^3 is the complementary base
  val A = 0
  val C = 1
  val G = 2
  val T = 3

  def toIndex(c: Char) : Int = c match {
    case 'A' => A
    case 'C' => C
    case 'G' => G
    case 'T' => T
    case _ => -1
  }  

  def toBase(c: Int) : Char = c match {
    case A => 'A'
    case C => 'C'
    case G => 'G'
    case T => 'T'
    case _ => 'N'
  }

  def toBase(c1: Int, c2: Int) = 'x'

  def complement(c: Char) = c match {
    case 'A' => 'T'
    case 'C' => 'G'
    case 'G' => 'C'
    case 'T' => 'A'
    case _ => c
  }
  
  
  def complementIndex(i: Int) = i ^ 3
  
  def reverseComplement(s: String) = {
    s map complement reverse    
  }

  def bit(i: Int) = 1 << i
  val Abit = bit(A)
  val Cbit = bit(C)
  val Gbit = bit(G)
  val Tbit = bit(T)
  val mapToBits =
    Map('A' -> Abit,
        'C' -> Cbit,
        'G' -> Gbit,
        'T' -> Tbit,
        'R' -> (Abit | Gbit), 
        'Y' -> (Cbit | Tbit),
        'S' -> (Gbit | Cbit),
        'W' -> (Abit | Tbit),
        'K' -> (Gbit | Tbit),
        'M' -> (Abit | Cbit),
        'B' -> (Cbit | Gbit | Tbit),
        'D' -> (Abit | Gbit | Tbit),
        'H' -> (Abit | Cbit | Tbit),
        'V' -> (Abit | Cbit | Gbit),
        'N' -> (Abit | Cbit | Gbit | Tbit))
  val mapToIUPAC = mapToBits map {x => (x._2 -> x._1)}
  
  def toIUPAC(base1: Char, base2: Char) = {
    val bits = bit(toIndex(base1)) | bit(toIndex(base2))
    mapToIUPAC(bits)
  }
    
}