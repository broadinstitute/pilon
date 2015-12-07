name := "pilon"

version := "1.16"

scalaVersion := "2.11.6"

scalacOptions += "-deprecation"

scalacOptions += "-feature"

Seq(com.github.retronym.SbtOneJar.oneJarSettings: _*)

libraryDependencies += "commons-lang" % "commons-lang" % "2.6"
