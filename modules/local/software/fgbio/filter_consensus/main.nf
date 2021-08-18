// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process FILTER_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
//    container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), file(bam) 
    path fasta
    val min_reads

    output:
    tuple val(meta), file("*_filt.bam")

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g FilterConsensusReads \\
    -i ${bam} \\
    -o ${meta.id}_cons_filt.bam \\
    -r ${fasta} \\
    --min-reads ${min_reads} \\
    --max-read-error-rate 0.05 \\
    --min-base-quality 30 \\
    --max-base-error-rate 0.1 \\
    --max-no-call-fraction 0.1 \\
    --reverse-per-base-tags true
    """
    stub:
    """
    touch ${meta.id}_cons_filt.bam
    """
}
