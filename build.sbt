name := "pilon"

version := "1.5"

scalaVersion := "2.10.2"

scalacOptions += "-deprecation"

seq(com.github.retronym.SbtOneJar.oneJarSettings: _*)

libraryDependencies += "commons-lang" % "commons-lang" % "2.6"
