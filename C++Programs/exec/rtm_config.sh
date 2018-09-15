#!/bin/bash

for i in {10..940..20}
do
	par1="/scratch/syauqy/result/NACD/image$i"
	par2="./velbp2004new.txt"
	par3="./setfile$i.txt"
	par4="/scratch/syauqy/record/SynBP2004/rec$i.txt"
	./rtm_nacd $par1 $par2 $par3 $par4
done
