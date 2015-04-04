#!/bin/sh
version=`grep ^version build.sbt |sed -e 's/.*\"\(.*\)\"/\\1/'`
date=`date`
commit=`git describe | sed -e 's/^v//'`
commitdate=`git log -n1 | grep '^Date' | sed -e 's/Date: *\(.*\)/\\1/'`
tmp=Version.scala.tmp
f=`find . -name Version.scala`
cp -p $f $tmp
sed -e "s/\(date.*=\).*/\\1 \"$commitdate\"/" \
    -e "s/\(commit.*=\).*/\\1 \"$commit\"/" \
    -e "s/\(sbt.*=\).*/\\1 \"$version\"/" \
    <$tmp >$f
#sbt $* package
sbt $* one-jar
ln -sf `pwd`/target/scala-2.11/pilon_2.11-$version-one-jar.jar ~/lib/pilon/pilon-dev.jar
mv $tmp $f
