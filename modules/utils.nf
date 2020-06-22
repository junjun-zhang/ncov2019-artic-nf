process prepareFastqPair {
    input:
        val sampleId
        path fastq_files

    output:
        tuple sampleId, path("${sampleId}_R1.fastq.gz"), path("${sampleId}_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        mkdir indir
        mv *.fastq.gz indir
        ls indir/*_{,R}1{,_001}.fastq.gz 2>/dev/null |xargs -I {} ln -s {} ${sampleId}_R1.fastq.gz || true
        ls indir/*_{,R}2{,_001}.fastq.gz 2>/dev/null |xargs -I {} ln -s {} ${sampleId}_R2.fastq.gz || true
        """
}

process performHostFilter {
    input:
        tuple(val(sampleId), path(forward), path(reverse))
    output:
        tuple sampleId, path("${sampleId}_hostfiltered_R1.fastq.gz"), path("${sampleId}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.human_ref} ${forward} ${reverse} | samtools view -q 30 -U ${sampleId}.host_unaligned.bam -Sb > ${sampleId}.host_aligned.bam
        samtools sort -n ${sampleId}.host_unaligned.bam | \
             samtools fastq -1 ${sampleId}_hostfiltered_R1.fastq.gz -2 ${sampleId}_hostfiltered_R2.fastq.gz -s ${sampleId}_singletons.fastq.gz -
        """
}


process getSampleId {
    input:
        path song_metadata

    output:
        stdout()

    """
    set -euxo pipefail
    cat ${song_metadata} | jq -er .sample.sample_id | tr -d '\\n'
    """
}
