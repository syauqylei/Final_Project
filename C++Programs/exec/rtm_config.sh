#!/bin/bash

for i in {20..440..20}
do
	par1="./scratch/syauqy/result/NACD/image$i"
	par2="./velbp2004.txt"
	par3="./setfile$i.txt"
	par4="./rec$i.txt"
	./Rtm_fd4 $par1 $par2 $par3 $par4
done
