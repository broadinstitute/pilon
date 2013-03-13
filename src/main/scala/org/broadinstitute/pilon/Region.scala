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
  def nearAny(others: List[Region]) = others.exists({_.near(this)})
  override def toString = 
    "<Region " + name + ":" + start + (if (size == 1) "" else "-" + stop + ">")
}

