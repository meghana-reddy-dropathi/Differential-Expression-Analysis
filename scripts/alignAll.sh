#!/usr/bin/env bash
# alignAll.sh
outDir='quant/'
sample='Aip02' # TODO: update to loop over all Aip## samples
function align {
    salmon quant -l IU \
        -1 /work/courses/BINF6309/AiptasiaMiSeq/fastq/${sample}.R1.fastq \
        -2 /work/courses/BINF6309/AiptasiaMiSeq/fastq/${sample}.R2.fastq \
        -i AipIndex \
        --validateMappings \
        -o ${outDir}${sample}
}

align 
