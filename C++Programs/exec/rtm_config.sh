#!/bin/bash

for i in {20..440..20}
do
	par1="./image$i"
	par2="./marm.txt"
	par3="./setfile$i.txt"
	par4="./rec$i.txt"
	./Rtm_nacdol $par1 $par2 $par3 $par4
done
