#!/bin/bash

for i in {20..440..20}
do
	par1="/scratch/syauqy/synmarmosi/image$i.txt"
	par2="./marm.txt"
	par3="./setfile$i.txt"
	par4="/scratch/syauqy/synmarmosi/rec$i.txt"
	./Rtm_fd $par1 $par2 $par3 $par4
done
