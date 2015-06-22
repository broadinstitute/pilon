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

import java.io.{ File, PrintWriter, FileWriter, BufferedWriter }

class Vcf(val file: File, val contigsWithSizes: List[(String, Int)] = Nil) {
  val writer = new PrintWriter(new BufferedWriter(new FileWriter(file)))
  val tab = "\t"

  def writeHeader = {
    val date = (new java.text.SimpleDateFormat("yyyyMMdd")).format(new java.util.Date())
    val ref = (new File(Pilon.genomePath)).toURI
    writer.println("##fileformat=VCFv4.1")
    writer.println("##fileDate=" + date)
    writer.println("##source=\"" + Version.version + "\"")
    writer.println("##PILON=\"" + Pilon.commandArgs.mkString(" ") + "\"")
    writer.println("##reference=" + ref)
    for ((c, s) <- contigsWithSizes)
      writer.println("##contig=<ID=" + c + ",length=" + s + ">")
    //writer.println("##FILTER=<ID=LowConf,Description=\"Low Confidence Call\">")
    writer.println("##FILTER=<ID=LowCov,Description=\"Low Coverage of good reads at location\">")
    //writer.println("##FILTER=<ID=LowMQ,Description=\"Low mean mapping quality at location\">")
    writer.println("##FILTER=<ID=Amb,Description=\"Ambiguous evidence in haploid genome\">")
    writer.println("##FILTER=<ID=Del,Description=\"This base is in a deletion or change event from another record\">")
    writer.println("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Valid read depth; some reads may have been filtered\">")
    writer.println("##INFO=<ID=TD,Number=1,Type=Integer,Description=\"Total read depth including bad pairs\">")
    writer.println("##INFO=<ID=PC,Number=1,Type=Integer,Description=\"Physical coverage of valid inserts across locus\">")
    writer.println("##INFO=<ID=BQ,Number=1,Type=Integer,Description=\"Mean base quality at locus\">")
    writer.println("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mean read mapping quality at locus\">")
    writer.println("##INFO=<ID=QD,Number=1,Type=Integer,Description=\"Variant confidence/quality by depth\">")
    writer.println("##INFO=<ID=BC,Number=4,Type=Integer,Description=\"Count of As, Cs, Gs, Ts at locus\">")
    if (Pilon.vcfQE)
      writer.println("##INFO=<ID=QE,Number=4,Type=Integer,Description=\"Evidence for As, Cs, Gs, Ts weighted by Q & MQ at locus\">")
    else
      writer.println("##INFO=<ID=QP,Number=4,Type=Integer,Description=\"Percentage of As, Cs, Gs, Ts weighted by Q & MQ at locus\">")
    writer.println("##INFO=<ID=IC,Number=1,Type=Integer,Description=\"Number of reads with insertion here\">")
    writer.println("##INFO=<ID=DC,Number=1,Type=Integer,Description=\"Number of reads with deletion here\">")
    writer.println("##INFO=<ID=XC,Number=1,Type=Integer,Description=\"Number of reads clipped here\">")
    writer.println("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">")
    writer.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Fraction of evidence in support of alternate allele(s)\">")
    writer.println("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
    writer.println("##INFO=<ID=SVLEN,Number=.,Type=String,Description=\"Difference in length between REF and ALT alleles\">")
    writer.println("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">")
    writer.println("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise change from local reassembly (ALT contains Ns)\">")
    writer.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    writer.println("##FORMAT=<ID=AD,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed\">")
    writer.println("##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Approximate read depth; some reads may have been filtered\">")
    writer.println("##ALT=<ID=DUP,Description=\"Possible segmental duplication\">")
    writer.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE")
  }

  writeHeader

  def close = writer.close

  def writeRecord(region: GenomeRegion, index: Int,
      embedded: Boolean = false, indelOk: Boolean = true): Unit = {
    val locus = region.locus(index)
    val pileUp = region.pileUpRegion(index)
    val bc = pileUp.baseCall
    val bcString = bc.callString(indelOk)
    val baseDP = bc.baseSum.toInt //pileUp.baseCount.sums(bc.baseIndex)
    val altBaseDP = bc.altBaseSum.toInt //pileUp.baseCount.sums(bc.altBaseIndex)
    val depth = pileUp.depth.toInt
    var loc = locus
    val (rBase, cBase, callType, refDP, altDP) = {
      if (indelOk && !embedded && bc.deletion) {
        loc -= 1
        val rBase = region.refBase(loc)
        //(rBase + bcString, rBase.toString, "1/1", depth - pileUp.deletions, pileUp.deletions)
        (rBase + bcString, rBase.toString, "1/1", pileUp.mqSum - pileUp.delQual, pileUp.delQual)
      } else if (indelOk && !embedded && bc.insertion) {
        loc -= 1
        val rBase = region.refBase(loc)
        //(rBase.toString, rBase + bcString, "1/1", depth - pileUp.insertions, pileUp.insertions)
        (rBase.toString, rBase + bcString, "1/1", pileUp.mqSum - pileUp.insQual, pileUp.insQual)
      } else if (bc.homo) {
        val rBase = region.refBase(loc)
        if (rBase == bc.base || bcString == "N")
          (rBase.toString, bc.base.toString, "0/0", baseDP, altBaseDP)
        else {
          (rBase.toString, bc.base.toString, "1/1", altBaseDP, baseDP)
        }
      } else {
        val rBase = region.refBase(loc)
        if (rBase == bc.base) {
          (rBase.toString, bc.altBase.toString, "0/1", baseDP, altBaseDP)
        } else {
          (rBase.toString, bc.base.toString, "0/1", altBaseDP, baseDP)
        }
      }
    }
    var filters = List[String]()
    if (depth < region.minDepth) filters ::= "LowCov"
    //if (!bc.highConfidence && !bc.indel) filters ::= "LowConf"
    if (!Pilon.diploid && !bc.homo && !(indelOk && bc.indel)) filters ::= "Amb"
    if (embedded) filters ::= "Del"
    if (filters.isEmpty) filters ::= "PASS"
    val cBaseVcf = if (cBase == "N" || cBase == rBase) "." else cBase
    val filter = filters.mkString(";")

    val ac = callType match {
      case "0/0" => 0
      case "0/1" => 1
      case "1/1" => 2
    }
    val af = if (refDP + altDP > 0 && cBaseVcf != ".")
      (altDP.toFloat / (refDP + altDP).toFloat)
    else 0.0

    val info = "DP=" + (if (!embedded) pileUp.depth else pileUp.count) +
    	";TD=" + (pileUp.depth + pileUp.badPair) +
    	";BQ=" + pileUp.meanQual +
    	";MQ=" + pileUp.meanMq +
    	";QD=" + bc.q +
    	";BC=" + pileUp.baseCount +
      (if (Pilon.vcfQE) ";QE=" + pileUp.qualSum.toString else ";QP=" + pileUp.qualSum.toStringPct ) +
    	";PC=" + pileUp.physCov +
      ";IC=" + pileUp.insertions +
      //";IF=" + pileUp.insPct +
      ";DC=" + pileUp.deletions +
      //";DF=" + pileUp.delPct +
      ";XC=" + pileUp.clips +
    	";AC=" + ac +
    	";AF=" + ("%.2f".format(af))

    val gt = "GT"
    val gtInfo = callType
    //AD removed
    //if (ac > 0) {
    //  gt += ":AD"
    //  gtInfo += ":" + refDP + "," + altDP 
    //}
    
    val line = region.name +
    	tab + loc + 
    	tab + "." + 
    	tab + rBase + 
    	tab + cBaseVcf +
    	tab + (if (indelOk && bc.deletion) "." else bc.score.toString) +
    	tab + filter +
    	tab + info +
    	tab + gt + 
    	tab + gtInfo

    writer.println(line)
    
    if (indelOk && bc.indel && !embedded) writeRecord(region, index, bc.deletion, false)
  }

  def writeFixRecord(region: GenomeRegion, fix: GenomeRegion.Fix) = {
    val loc = fix._1 - 1
    val rBase = region.refBase(loc)
    val ref = rBase + fix._2
    val alt = rBase + fix._3
    val svlen = alt.length - ref.length
    val svend = loc + ref.length - 1
    val svtype = if (svlen < 0) "DEL" else "INS"
    var line = region.name + tab + loc + tab + "." + tab
    line += ref + tab + alt + tab + "." + tab + "PASS" + tab
    line += "SVTYPE=" + svtype + ";SVLEN=" + svlen + ";END=" + svend
    if (fix._3 contains 'N') line += ";IMPRECISE"
    line += tab + "GT" + tab + "1/1"
    writer.println(line)
  }

  def writeDup(region: GenomeRegion, dup: Region) = {
    var loc = dup.start - 1
    val rBase = region.refBase(loc)
    var line = region.name + tab + loc + tab + "." + tab
    line += rBase + tab + "<DUP>" + tab + "." + tab + "PASS" + tab
    line += "SVTYPE=DUP;SVLEN=" + dup.size + ";END=" + dup.stop + ";IMPRECISE"
    line += tab + "GT" + tab + "./."
    writer.println(line)
  }
}
