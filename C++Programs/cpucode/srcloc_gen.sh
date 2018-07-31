#!/bin/bash

for i in {5..7}
do
	par1="./rec$i"
	par2="./vel$i.txt"
	par3="./fwdset.txt"
	./wve_nacdgpu $par1 $par2 $par3
done
