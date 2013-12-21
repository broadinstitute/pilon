name := "pilon"

version := "1.6"

scalaVersion := "2.10.2"

scalacOptions += "-deprecation"

scalacOptions += "-feature"

seq(com.github.retronym.SbtOneJar.oneJarSettings: _*)

libraryDependencies += "commons-lang" % "commons-lang" % "2.6"
