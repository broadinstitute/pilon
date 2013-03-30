/*
 * Copyright (c) 2012, 2013. The Broad Institute, Inc.
 *
 * This file is part of Pilon.
 *
 * Pilon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2
 * as published by the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License
 * along with Pilon.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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

