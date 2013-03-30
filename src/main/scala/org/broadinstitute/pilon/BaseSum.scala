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


class BaseSum {
  val sums = new Array[Long](4)

  def add(base: Int, n: Int = 1) = sums(base) += n

  def remove(base: Int, n: Int = 1) = sums(base) -= n

  def sum = sums.sum

  def maxBase = {
    var maxIndex = 0
    var max = sums(maxIndex)
    var sum = max

    for (b <- 1 to 3) {
      var n = sums(b)
      sum += n
      if (n > max) {
        max = n
        maxIndex = b
      }
    }
    (maxIndex, max, sum)
  }

  def excessBase = {
    val (base, max, total) = maxBase
    var excess = max
    for (b <- 0 to 3) {
    	if (b != base) excess -= sums(b)
    }
    (base, excess, total)
  }
  
  def order = {
    val indicies = List(0,1,2,3).toArray
    indicies.sortWith({(a,b) => sums(a) > sums(b)})
  }
  
  def /(divisor : Int) = {
    val bs = new BaseSum
    for (i <- 0 to 3) if (divisor != 0) bs.sums(i) = this.sums(i) / divisor else 0
    bs
  }
  
  def toStringPct = {
    val div = sum
    sums map ({x => if (div == 0) 0 else (100*x + div/2) / div}) mkString(",")
  }

  override def toString = sums.mkString(",")
}