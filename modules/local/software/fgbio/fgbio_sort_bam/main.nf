// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../../../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process FGBIO_SORT_BAM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    //container "quay.io/biocontainers/fgbio:1.3.0--0"

    input:
    tuple val(meta), file(bam)

    output:
    tuple val(meta), file("*sort*.bam")

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g SortBam \\
    -i $bam \\
    -o ${meta.id}_sort.bam \\
    --sort-order TemplateCoordinate \\
    --max-records-in-ram 4000000
    """
    stub:
    """
    touch ${meta.id}_sort.bam
    """
    }