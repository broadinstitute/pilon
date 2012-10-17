#!/bin/sh

tmp=Version.scala.tmp
date=`date`
svn=`svnversion`
f=`find . -name Version.scala`
cp -p $f $tmp
sed -e "s/\(date.*=\).*/\\1 \"$date\"/"  -e "s/\(svn.*=\).*/\\1 \"$svn\"/" <$tmp >$f
#sbt $* package
sbt $* one-jar
cp -p target/scala-2.9.2/pilon_2.9.2-1.0-one-jar.jar ~/lib/pilon/pilon.jar
mv $tmp $f
