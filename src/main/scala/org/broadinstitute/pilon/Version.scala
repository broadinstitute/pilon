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


//
// This is a generic version file which is updated by a build script.
//

object Version {
	val date = "[DATE NOT SET]"
	val svn = "[VERSION NOT SET]"
	def version = "Pilon version " + svn + " " + date
}
