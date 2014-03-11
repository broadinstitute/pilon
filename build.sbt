name := "pilon"

version := "1.7"

scalaVersion := "2.10.3"

scalacOptions += "-deprecation"

scalacOptions += "-feature"

seq(com.github.retronym.SbtOneJar.oneJarSettings: _*)

libraryDependencies += "commons-lang" % "commons-lang" % "2.6"
