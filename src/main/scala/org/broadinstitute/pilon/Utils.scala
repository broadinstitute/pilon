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

object Utils {
  def roundDiv(n: Long, d: Long) = if (d > 0) (n + d/2) / d else 0.toLong
  def roundDiv(n: Int, d: Int) = if (d > 0) (n + d/2) / d else 0
  def pct(n: Long, d: Long) = roundDiv(100 * n, d)
  def pct(n: Int, d: Int) = roundDiv(100 * n, d)
}
