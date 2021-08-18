/*
===================================
          UMI stage one
===================================
*/

include { FGBIO_SORT_BAM }             from '../../software/fgbio/fgbio_sort_bam/main'          addParams(options: params.fgbio_sort_mapping_options)
include { CALL_CONSENSUS }             from '../../software/fgbio/call_consensus/main'          addParams(options: params.fgbio_call_consensus_mapping_options)
include { FILTER_CONSENSUS }           from '../../software/fgbio/filter_consensus/main'        addParams(options: params.fgbio_filter_mapping_options)

workflow UMI_STAGE_ONE_MINI {
    take:
    input_samples

    main:
    FGBIO_SORT_BAM(input_samples)
    CALL_CONSENSUS(FGBIO_SORT_BAM.out)
    FILTER_CONSENSUS(CALL_CONSENSUS.out,fasta, min_reads)

    emit:
    filtered_bam  = FILTER_CONSENSUS.out
}