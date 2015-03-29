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

/**
 * Created by IntelliJ IDEA.
 * User: bruce
 * Date: 11/20/11
 * Time: 5:19 PM
 * To change this template use File | Settings | File Templates.
 */
import java.io._

import scala.collection.JavaConversions._
import htsjdk.samtools._

class Tracks(val reference: GenomeFile, val prefix : String = "") {
  def standardTracks = {
    makeBedTrack("Pilon.bed", "Pilon")
    changesTrack("Changes.wig")
    unconfirmedTrack("Unconfirmed.wig")
    copyNumberTrack("CopyNumber.wig")
    coverageTrack("Coverage.wig")
    badCoverageTrack("BadCoverage.wig")
    pctBadTrack("PctBad.wig")
    //coverageTrackSD("CoverageSD.wig")
    //badCoverageTrackSD("BadCoverageSD.wig")
    deltaCoverageTrack("DeltaCoverage.wig")
    dipCoverageTrack("DipCoverage.wig")
    //fragCoverageTrack("FragCoverage.wig")
    physicalCoverageTrack("PhysicalCoverage.wig")
    //physicalCoverageTrackSD("PhysicalCoverageSD.wig")
    //insertSizeTrack("InsertSize.wig")
    //insertSizeTrackSD("InsertSizeSD.wig")
    clippedAlignmentTrack("ClippedAlignments.wig")
    weightedQualTrack("WeightedQual.wig")
    weightedMqTrack("WeightedMq.wig")
    gcTrack("GC.wig")
    //kmerCopyNumberTrack("KmerCopyNumber.wig")
  }
  
  def changesTrack(file: String) = {
    makeTrack(file, "Changes", 
        { (r: GenomeRegion, i: Int) => if (r.changed(i)) 1 else 0 })
  }

  def unconfirmedTrack(file: String) = {
    makeTrack(file, "Unconfirmed", 
        { (r: GenomeRegion, i: Int) => if (r.confirmed(i)) 0 else 1 })
  }
  
  def copyNumberTrack(file: String) = {
    makeTrack(file, "Copy Number", 
        { (r: GenomeRegion, i: Int) => r.copyNumber(i) - 1 })
  }
  
  def kmerCopyNumberTrack(file: String) = {
    makeTrack(file, "Kmer Copy Number", 
        { (r: GenomeRegion, i: Int) => (r.kmerCopyNumber(i) - 1) max 0})
  }
  
  def coverageTrack(file: String) = {
    makeTrack(file, "Coverage",
        { (r: GenomeRegion, i: Int) => r.coverage(i) })
  }

  def fragCoverageTrack(file: String) = {
    makeTrack(file, "Frag Coverage",
        { (r: GenomeRegion, i: Int) => r.fragCoverage(i) })
  }

  def coverageTrackSD(file: String) = {
    makeTrack(file, "Coverage SD",
        { (r: GenomeRegion, i: Int) => r.coverageDist.toSigma10x(r.coverage(i)) },
        "viewLimits=-30:30")
  }

  def badCoverageTrack(file: String) = {
    makeTrack(file, "Bad Coverage", 
        { (r: GenomeRegion, i: Int) => r.badCoverage(i) })
  }

  def badCoverageTrackSD(file: String) = {
    makeTrack(file, "Bad Coverage SD",
        { (r: GenomeRegion, i: Int) => r.badCoverageDist.toSigma10x(r.badCoverage(i)) },
        "viewLimits=-30:30")
  }

  def deltaCoverageTrack(file: String, radius: Int = 100) = {
    makeTrack(file, "Delta Coverage", 
        { (r: GenomeRegion, i: Int) => r.deltaCoverage(i, radius) })
  }

  def dipCoverageTrack(file: String, radius: Int = 100) = {
    makeTrack(file, "Dip Coverage", 
        { (r: GenomeRegion, i: Int) => r.dipCoverage(i, radius) })
  }
    
  def physicalCoverageTrack(file: String) = {
    makeTrack(file, "Physical Coverage", 
    		{ (r: GenomeRegion, i: Int) => r.physCoverage(i) })
  }
  def physicalCoverageTrackSD(file: String) = {
    makeTrack(file, "Physical Coverage SD",
    		{ (r: GenomeRegion, i: Int) => r.physCoverageDist.toSigma10x(r.physCoverage(i)) },
        "viewLimits=-30:30")
  }

  def gcTrack(file: String) = {
    makeTrack(file, "GC", 
        { (r: GenomeRegion, i: Int) => r.gc(i) },
        "graphType=heatmap midRange=35:65 midColor=0,255,0")
  }
  
  def insertSizeTrack(file: String) = {
    makeTrack(file, "Insert Size", 
        { (r: GenomeRegion, i: Int) => r.insertSize(i) })
  }

  def insertSizeTrackSD(file: String) = {
    makeTrack(file, "Insert Size SD",
        { (r: GenomeRegion, i: Int) => r.insertSizeDist.toSigma10x(r.insertSize(i)) },
        "viewLimits=-30:30")
  }

  def pctBadTrack(file: String) = {
    makeTrack(file, "Pct Bad", 
        { (r: GenomeRegion, i: Int) =>
          val good = r.coverage(i)
          val bad = r.badCoverage(i)
          if (good+bad > 0) bad * 100 / (good + bad)
          else 0 
        })
  }


  def weightedMqTrack(file: String) = {
    makeTrack(file, "Weighted MQ", 
        //{ (r: GenomeRegion, i: Int) => r.weightedMqDist.toSigma10x(r.weightedMq(i)) },
        //"viewLimits=-30:30")
    { (r: GenomeRegion, i: Int) => r.weightedMq(i) })
  }

  def weightedQualTrack(file: String) = {
    makeTrack(file, "Weighted Qual", 
    { (r: GenomeRegion, i: Int) => r.weightedQual(i) })
  }


  def clippedAlignmentTrack(file: String) = {
    makeTrack(file, "Clipped Alignments",
    { (r: GenomeRegion, i: Int) => r.clips(i) })
  }

  def prefixedFile(fileName: String) = {
    val prefixedFileName = if (prefix == "") fileName else prefix + fileName
    new File(prefixedFileName)
  }

  def makeTrack(fileName: String, name: String, func: (GenomeRegion, Int) => Int, options: String = "") = {
    val file = prefixedFile(fileName)
    println ("Creating " + name + " track in file " + file.getPath())
    val writer = new PrintWriter(file)
    var headLine = "track type=wiggle_0 graphType=line color=0,0,255 altColor=255,0,0 name=\"" + name + "\""
    if (options != "") headLine += " " + options
    writer.println(headLine)
    for ((cName, regions) <- reference.regions) {
    	writer.println("fixedStep chrom=" + cName + " start=1 step=1")
    	regions foreach { region =>
    	  for (rIndex <- 0 to region.size-1) {
    	    val value = func(region, rIndex)
    	    writer.println(value)
    	  }
    	}
    }
    writer.close()   
  }
  
  def regionsToBed(regions: List[Region], name: String, writer: PrintWriter, rgb: String = "0,255,0") = {
    regions map {region: Region =>
      List(region.name, region.start, region.stop, name, "0", "+", region.start-1, region.stop, rgb) mkString("\t")
    } foreach { writer.println(_) }
  }
  
  def makeBedTrack(fileName: String, name: String, options: String = "") = {
    val file = prefixedFile(fileName)
    println ("Creating " + name + " track in file " + file.getPath())
    val writer = new PrintWriter(file)
    var headLine = "track description=\"Issues found by Pilon\" name=\"" + name + "\""
    if (options != "") headLine += " " + options
    writer.println(headLine)
    
    for ((cName, regions) <- reference.regions) {
    	regions foreach { r: GenomeRegion =>
    	  regionsToBed(r.unConfirmedRegions, "?", writer, "255,0,0") }
    	regions foreach { r: GenomeRegion =>
    	  regionsToBed(r.possibleCollapsedRepeats, "#", writer, "0,0,255") }
    	regions foreach { r: GenomeRegion =>
    	  regionsToBed(r.changeRegions, "X", writer, "255,0,0") }
    	regions foreach { r: GenomeRegion =>
    	  regionsToBed(r.possibleInsertions, "I", writer, "255,255,0") }
        regions foreach { r: GenomeRegion =>
          regionsToBed(r.possibleDeletions, "D", writer, "255,0,255") }
        regions foreach { r: GenomeRegion =>
          regionsToBed(r.gaps, "G", writer, "32,32,32") }
        regions foreach { r: GenomeRegion =>
          regionsToBed(r.possibleBreaks, "B", writer, "0,255,255") }
    }
    writer.close()       
  }
}
