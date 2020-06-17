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
        ls indir/*_R1{,_001}.fastq.gz 2>/dev/null |xargs -I {} ln -s {} ${sampleId}_R1.fastq.gz || true
        ls indir/*_R2{,_001}.fastq.gz 2>/dev/null |xargs -I {} ln -s {} ${sampleId}_R2.fastq.gz || true
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