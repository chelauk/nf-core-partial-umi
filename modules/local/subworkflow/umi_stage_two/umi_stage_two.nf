params.second_sort = true

include { BED_TO_INTERVAL_LIST }       from '../../software/picard/bed_to_interval_list/main'   addParams(options: params.bed_to_intervals_options)
include { BAM_TO_FASTQ }                    from '../../software/picard/bam_to_fastq/main'      addParams(options: params.bam_to_fastq_options)
include { BWA_ALN }                         from '../../software/bwa/bwa_aln/main'              addParams(options: params.bwamem1_mem_options)
include { PICARD_SORT_BAM }                 from '../../software/picard/picard_sort_bam/main'   addParams(options: params.picard_sort_mapping_options)
include { PICARD_SORT_BAM as SORT_BAM_TWO } from '../../software/picard/picard_sort_bam/main'   addParams(options: params.picard_sort_mapping_options, second_file: true )
include { PICARD_MERGE_BAMS }               from '../../software/picard/picard_merge_bams/main' addParams(options: params.picard_merge_bams_options)
include { PICARD_COLLECT_HS_METRICS_2 }     from '../../software/picard/collect_hs_metrics_2/main'     addParams(options: params.collect_hs_metrics_options)
include { MARK_DUPLICATES }                 from '../../software/gatk/markduplicates/main'      addParams(options: params.gatk_mark_duplicates_options)
include { ERRORRATE_BY_READ_POSITION }      from '../../software/fgbio/error_rate/main'         addParams(options: params.error_rate_options, second_file: true)

workflow UMI_STAGE_TWO {

    take:
    filtered_bam
    bwa_index
    fasta
    fasta_fai
    target_bed
    dict
    dbsnp
    dbsnp_index

    main:
    BED_TO_INTERVAL_LIST(target_bed, dict)
    BAM_TO_FASTQ(filtered_bam)
    BWA_ALN(BAM_TO_FASTQ.out,bwa_index,fasta,fasta_fai)
    PICARD_SORT_BAM(BWA_ALN.out.bam,fasta,dict)
    SORT_BAM_TWO(filtered_bam,fasta,dict)
    PICARD_MERGE_BAMS(PICARD_SORT_BAM.out.join(SORT_BAM_TWO.out),fasta,dict)
    MARK_DUPLICATES(PICARD_MERGE_BAMS.out.merged_bam)
    PICARD_COLLECT_HS_METRICS_2( MARK_DUPLICATES.out.md_bam, MARK_DUPLICATES.out.md_bai,fasta, fasta_fai, dict, BED_TO_INTERVAL_LIST.out.interval_list)
    md_bam = MARK_DUPLICATES.out.md_bam
    md_tsv = MARK_DUPLICATES.out.tsv

    md_tsv.collectFile(storeDir: "${params.outdir}/preprocessing/tsv") { meta ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam = "${launchDir}/${params.outdir}/preprocessing/mark_duplicates/${patient}_${sample}.md.bam"
            bai = "${launchDir}/${params.outdir}/preprocessing/mark_duplicates/${patient}_${sample}.md.bam.bai"
            ["md_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\n"]
    }

    md_tsv.map { meta ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam = "${launchDir}/${params.outdir}/preprocessing/mark_duplicates/${patient}_${sample}.md.bam"
            bai = "${launchDir}/${params.outdir}/preprocessing/mark_duplicates/${patient}_${sample}.md.bam.bai"
            "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\n"
    }.collectFile(name: 'md.tsv', sort: true, storeDir: "${params.outdir}/preprocessing/tsv") 
    ERRORRATE_BY_READ_POSITION(MARK_DUPLICATES.out.md_bam,fasta,dict,dbsnp,dbsnp_index,BED_TO_INTERVAL_LIST.out.interval_list)

    emit:
    error_rate_2  = ERRORRATE_BY_READ_POSITION.out
    md_hs_metrics = PICARD_COLLECT_HS_METRICS_2.out.md_hs_metrics
    md_report     = MARK_DUPLICATES.out.report

}
