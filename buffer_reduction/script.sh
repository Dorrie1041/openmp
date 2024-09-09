#!/bin/bash

if [ "$#" -ne 1 ]; then
   echo "Usage: $0 <file_number>"
   exit 1
fi
NUMBER=$1

f="../dCNS_test/dCNS_setaria_maize_V23/"$NUMBER
c=100
for i in 1, 2, 4, 8, 16
do
    srun -n 1 -c 128 --cpu_bind=cores numactl --interleave=all ./diag_buffer_compared ${c} ${i} ${f}
done

    
	 



