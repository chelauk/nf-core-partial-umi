// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FILTER_UMIS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*trimmed_?.fq.gz"), emit: reads

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ln -s ${reads[0]} ${prefix}_trimmed_1.fq.gz
    ln -s ${reads[2]} ${prefix}_trimmed_3.fq.gz
    filter_umis.py -u ${prefix}_2.fq.gz -v ${prefix}_1_val_1.fq.gz
    """
    
    stub:
    """
    touch fastq_trimmed_1.fq.gz
    touch fastq_trimmed_2.fq.gz
    touch fastq_trimmed_3.fq.gz
    """
}
