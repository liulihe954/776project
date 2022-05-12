#!/bin/bash

sample=(8252 8867 6334 8862 8966 8897 8832 8971 8514 8983 8896 8841)
for i in "${sample[@]}" 
do 
        sbatch counts.sbatch $i
	sleep 1
done

