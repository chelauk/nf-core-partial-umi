nextflow.enable.dsl=2

/*
=================================
          PRINT HELP
=================================
*/

def json_schema = "$baseDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

/*
=================================
       PARAMETER SUMMARY
=================================
*/

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

/*
=================================
       PARAMETER CHECKS
=================================
*/

Checks.aws_batch(workflow, params) // Check AWS batch settings
Checks.hostname(workflow, params, log)  // Check the hostnames against configured profiles

// MultiQC - Stage config files

multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)
/*
================================================================================
                     UPDATE MODULES OPTIONS BASED ON PARAMS
================================================================================
*/
modules = params.modules
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
params.bwa                     = params.genome ? params.genomes[params.genome].bwa         ?: false : false
params.dbsnp                   = params.genome ? params.genomes[params.genome].dbsnp       ?: false : false
params.dbsnp_index             = params.genome ? params.genomes[params.genome].dbsnp_index ?: false : false
params.dict                    = params.genome ? params.genomes[params.genome].dict        ?: false : false
params.fasta                   = params.genome ? params.genomes[params.genome].fasta       ?: false : false
params.fasta_fai               = params.genome ? params.genomes[params.genome].fasta_fai   ?: false : false
file("${params.outdir}/no_file").text = "no_file\n"

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
bwa_index         = params.bwa               ? file(params.bwa)               : file("${params.outdir}/no_file")
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : file("${params.outdir}/no_file")
dbsnp_index       = params.dbsnp_index       ? file(params.dbsnp_index)       : file("${params.outdir}/no_file")
dict              = params.dict              ? file(params.dict)              : file("${params.outdir}/no_file")
fasta             = params.fasta             ? file(params.fasta)             : file("${params.outdir}/no_file")
fasta_fai         = params.fasta_fai         ? file(params.fasta_fai)         : file("${params.outdir}/no_file")
target_bed        = params.target_bed        ? file(params.target_bed)        : file("${params.outdir}/no_file")

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
min_reads         = params.min_reads         ?: Channel.empty()  // this the minimum reads parameter passed to FilterConsensusReads
read_structure    = params.read_structure    ?: Channel.empty()
params.enable_conda = false
params.second_file  = false

//  include functions
include {
    extract_fastq;
    extract_bam;
    has_extension
} from './modules/functions'


// handle bam input for mini
bam_tsv_path = null
if (params.input && (has_extension(params.input, "tsv")) && (params.stage == "mini") ) bam_tsv_path = params.input
if (bam_tsv_path) {
    bam_tsv_file = file(bam_tsv_path)
    input_samples = extract_bam(bam_tsv_file)
    }

// params summary for MultiQC
workflow_summary = Schema.params_summary_multiqc(workflow, summary_params)
workflow_summary = Channel.value(workflow_summary)


include { UMI_STAGE_ONE_MINI } from './modules/local/subworkflow/umi_stage_one_mini/umi_stage_one_mini' addParams(
    fgbio_sort_mapping_options:           modules['fgbio_sort_mapping'],
    fgbio_call_consensus_mapping_options: modules['fgbio_call_consensus_mapping'],
    fgbio_filter_mapping_options:         modules['fgbio_filter_mapping']
)


workflow {

    if ( params.stage == 'mini' ) { bam_grouped_by_umi = input_samples }
    UMI_STAGE_ONE_MINI(bam_grouped_by_umi,fasta, fasta_fai, dict)

}
