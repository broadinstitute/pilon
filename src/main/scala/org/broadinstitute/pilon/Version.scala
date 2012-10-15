package org.broadinstitute.pilon


//
// This is a generic version file which is updated by a build script.
//

object Version {
	val date = "[DATE NOT SET]"
	val svn = "[VERSION NOT SET]"
	def version = "Pilon version " + svn + " compiled " + date
}
