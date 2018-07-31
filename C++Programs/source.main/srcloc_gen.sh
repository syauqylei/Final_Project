#!/bin/bash

for i in {1..5}
do
	par1="./rec$i"
	par2="./vel$i.txt"
	par3="./fwdset.txt"
	./wve_nacdgpu $par1 $par2 $par3
done
