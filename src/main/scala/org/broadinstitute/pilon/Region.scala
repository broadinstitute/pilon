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

class Region(val name: String, val start: Int, val stop: Int) {
  def inRegion(locus: Int) = locus >= start && locus <= stop
  def beforeRegion(locus: Int) = locus < start
  def afterRegion(locus: Int) = locus > stop
  def index(locus: Int) = locus - start
  def locus(index: Int) = start + index
  def size = stop + 1 - start
  def midpoint = (start + stop) / 2

  def equals(other: Region) = 
    start == other.start && stop == other.stop && name == other.name
  def overlaps(other: Region) =
    other.start <= stop && other.stop >= start && other.name == name
  def contains(other: Region) = 
    other.start >= start && other.stop <= stop && other.name == name
    
  def near(other: Region, distance: Int = 100) =
    other.name == name && 
      (overlaps(other) || (other.stop - start).abs <= distance || (other.start - stop).abs <= distance)
  def nearAny(others: List[Region], distance: Int = 100) = others.exists({_.near(this, distance)})
  def regionString = name + ":" + start + (if (size < 2) "" else "-" + stop)
  //override def toString = "<Region " + regionString + ">"
  override def toString = regionString
}

