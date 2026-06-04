nextflow.enable.dsl=2

params.ref      = 'ref/GCF_000005845.2_ASM584v2_genomic.fna'
params.read1    = 'data/SRR9075672_1.fastq'
params.read2    = 'data/SRR9075672_2.fastq'
params.sample   = 'SRR9075672'
params.threshold = 90
params.outdir   = 'results'


process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "*.html", emit: html
    path "*.zip",  emit: zip

    script:
    """
    fastqc -o . $read1 $read2
    """
}


process ALIGN {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), path(ref), path(read1), path(read2)

    output:
    tuple val(sample), path(ref), path("${sample}.sam"), path("${sample}.bam"), emit: aligned

    script:
    """
    bwa index $ref
    samtools faidx $ref

    bwa mem $ref $read1 $read2 > ${sample}.sam
    samtools view -bS ${sample}.sam > ${sample}.bam
    """
}


process FLAGSTAT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), path(ref), path(sam), path(bam)

    output:
    tuple val(sample), path(ref), path(bam), path("${sample}.flagstat.txt"), emit: stats

    script:
    """
    samtools flagstat $bam > ${sample}.flagstat.txt
    """
}


process PARSE_FLAGSTAT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), path(ref), path(bam), path(flagstat)
    val threshold

    output:
    tuple val(sample), path(ref), path(bam), path("${sample}.parsed.txt"), path("${sample}.mapped.int"), emit: parsed

    script:
    """
    mapped=\$(grep "mapped (" $flagstat | head -n 1 | awk -F '[()%]' '{print \$2}')
    mapped_int=\${mapped%.*}

    if [ "\$mapped_int" -ge "$threshold" ]; then
        status="OK"
    else
        status="NOT OK"
    fi

    {
      echo "Mapped: \${mapped}%"
      echo "\${status}"
    } > ${sample}.parsed.txt

    echo "\${mapped_int}" > ${sample}.mapped.int
    """
}


process SORT_BAM {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), path(ref), path(bam)

    output:
    tuple val(sample), path(ref), path("${sample}.sorted.bam"), emit: sorted

    script:
    """
    samtools sort $bam -o ${sample}.sorted.bam
    """
}


process FREEBAYES {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), path(ref), path(sorted_bam)

    output:
    path "${sample}.variants.vcf", emit: vcf

    script:
    """
    freebayes -f $ref $sorted_bam > ${sample}.variants.vcf
    """
}


workflow {
    def ref = file(params.ref)

    def reads_ch = Channel.of(
        tuple(params.sample, file(params.read1), file(params.read2))
    )

    FASTQC(reads_ch)

    def align_in = reads_ch.map { sample, r1, r2 ->
        tuple(sample, ref, r1, r2)
    }

    ALIGN(align_in)

    FLAGSTAT(ALIGN.out.aligned)

    PARSE_FLAGSTAT(FLAGSTAT.out.stats, params.threshold)

    def qc = PARSE_FLAGSTAT.out.parsed.map { sample, ref2, bam, parsed_file, mapped_file ->
        def pct = mapped_file.text.trim().toInteger()
        tuple(sample, ref2, bam, parsed_file, pct)
    }

    def split = qc.branch {
        ok:     it[4] >= params.threshold
        not_ok: it[4] < params.threshold
    }

    SORT_BAM(
        split.ok.map { sample, ref2, bam, parsed_file, pct ->
            tuple(sample, ref2, bam)
        }
    )

    FREEBAYES(SORT_BAM.out.sorted)

    split.not_ok.view { row ->
        println "NOT OK for ${row[0]}: ${row[4]}% mapped"
    }
}
