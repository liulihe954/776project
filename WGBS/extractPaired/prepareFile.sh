#!/bin/bash

sample=(M6125 M6126 M6152 M6153 M6157 M6197 M6226 M6254 M6412 M6415 M6419 M6438 M6453 M6472 M6485 M6520)

for i in "${sample[@]}" 
do 
       cat $i*deduplicated.bismark.cov > /ufrc/penagaricano/fpenagaricano/Methylation/fpenagaricano/fpenagaricano.GSE117194/$i."final.txt"
done

