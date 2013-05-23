name := "pilon"

version := "1.0rc"

scalaVersion := "2.9.2"

scalacOptions += "-deprecation"

seq(com.github.retronym.SbtOneJar.oneJarSettings: _*)

libraryDependencies += "commons-lang" % "commons-lang" % "2.6"
