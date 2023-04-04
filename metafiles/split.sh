#!/bin/sh

file=$1
prefix=$2
n=`cat $file | wc -l`
line_num=`seq 1 1 $n`

cat $file | csplit - `echo $line_num` -f ${prefix}

rm ${prefix}00
ls ${prefix}0* | while read file
do mv $file ${prefix}${file#${prefix}0}
done

