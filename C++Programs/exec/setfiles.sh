#!/bin/bash

touch setfile{20..440..20}.txt
for f in {20..440..20}
do
	cat fwdset.txt > "setfile$f.txt"
done
for i in {20..440..20}
do
	echo "$i" >> "setfile$i.txt"
	echo "0" >> "setfile$i.txt"
done

