#!/bin/bash
if [ "$#" -ne 6 ]; then
	echo "Usage: $0 <m> <n> <k> <threads> <w> <h>"
	exit 1
fi

M_NUM=$1
N_NUM=$2
K_NUM=$3
THREADS=$4
W_NUM=$5
H_NUM=$6

~/openmp/openmp_exercise/openmp_code $M_NUM $N_NUM $K_NUM $THREADS $W_NUM $H_NUM

exit 0
