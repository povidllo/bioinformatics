REF=ref/GCF_000005845.2_ASM584v2_genomic.fna

fastqc data/SRR9075672_1.fastq data/SRR9075672_2.fastq

bwa index $REF
samtools faidx $REF

bwa mem $REF data/SRR9075672_1.fastq data/SRR9075672_2.fastq > results/alignment.sam

samtools view -bS results/alignment.sam > results/alignment.bam

samtools flagstat results/alignment.bam > results/flagstat.txt 

mapped=$(grep "mapped (" results/flagstat.txt | head -n 1 | awk -F '[()%]' '{print $2}')
mapped_int=${mapped%.*}

echo "Mapped: ${mapped}%" > results/parsed.txt

THRESHOLD=90

if [ "$mapped_int" -ge "$THRESHOLD" ]; then
    echo "OK" >> results/parsed.txt

    samtools sort results/alignment.bam -o results/alignment.sorted.bam

    freebayes -f $REF results/alignment.sorted.bam > results/variants.vcf

    echo "Pipeline finished successfully"
else
    echo "NOT OK" >> results/parsed.txt
    echo "Mapping quality too low (<90%)"
fi
