/*
 * Copyright 2012-2015 Broad Institute, Inc.
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
    s.map(complement).reverse
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

  def baseMatch(iupac: Char, base: Char) = {
    // returns true if "base" is a subset of iupac
    (iupac == base) || {
      val iuBits = mapToBits(iupac)
      val baseBits = mapToBits(base)
      (iuBits & baseBits) == baseBits
    }
  }
    
}


