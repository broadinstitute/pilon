package org.broadinstitute.pilon

object Utils {
  def roundDiv(n: Long, d: Long) = if (d > 0) (n + d/2) / d else 0
  def pct(n: Long, d: Long) = roundDiv(100 * n, d)
}