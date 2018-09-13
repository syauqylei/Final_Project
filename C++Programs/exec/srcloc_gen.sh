#!/bin/bash

touch setfile{10..940..20}.txt
for f in {10..940..20}
do
	cat fwdset.txt > "setfile$f.txt"
done
for i in {10..940..20}
do
	echo "$i" >> "setfile$i.txt"
	echo "0" >> "setfile$i.txt"
done

for i in {10..940..20}
do
	par1="./rec$i"
	par2="./marm.txt"
	par3="./setfile$i.txt"
	./wve_nacdgpu $par1 $par2 $par3
done
