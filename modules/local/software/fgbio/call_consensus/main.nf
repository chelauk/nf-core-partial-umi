// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process CALL_CONSENSUS {
    tag "$meta.id"
    label 'CALL_CONSENSUS'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    //container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*consensus.bam")

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g -XX:+AggressiveOpts -XX:+AggressiveHeap CallMolecularConsensusReads \\
    -i $bam \\
    -o ${meta.id}_consensus.bam \\
    --min-reads 1 \\
    --min-input-base-quality 30 \\
    --tag MI
    """
    stub:
    """
    touch ${meta.id}_consensus.bam
    """
}
