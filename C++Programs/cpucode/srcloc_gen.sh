#!/bin/bash
for number in {10..1110..50}
do
echo "$number" >> fwdset.txt
done
for num in {10..1110..50}
do
echo '0' >> fwdset.txt
done
