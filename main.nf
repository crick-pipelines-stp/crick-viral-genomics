#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { params_summary_map   } from './modules/local/util/logging/main'
include { summary_log          } from './modules/local/util/logging/main'
include { multiqc_summary      } from './modules/local/util/logging/main'
include { get_genome_attribute } from './modules/local/util/references/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = file("$projectDir/assets/multiqc_illumina_config.yml", checkIfExists: true)
ch_multiqc_logo = file("$projectDir/assets/The_Francis_Crick_Institute_logo.png", checkIfExists: true)
ch_seq_sim_config = file(params.seq_sim_config, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def host_fasta = get_genome_attribute(params, 'fasta')
def host_bwa   = get_genome_attribute(params, 'bwa'  )
if(params.host_fasta) {
    host_fasta = params.host_fasta
}
if(params.host_bwa) {
    host_bwa = params.host_bwa
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info summary_log(workflow, params, params.debug, params.monochrome_logs)
def summary_params = params_summary_map(workflow, params, params.debug)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    viral_fasta: params.viral_fasta
]
for (param in check_param_list) {
    if (!param.value) {
        exit 1, "Required parameter not specified: ${param.key}"
    }
    else {
        file(param.value, checkIfExists: true)
    }
}

// If no data being generated, samplesheet is manditory
if(!params.generate_reads && params.samplesheet == null) {
     exit 1, "Required parameter not specified: samplesheet"
}

// Check non-manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    params.samplesheet,
    params.viral_gff,
    params.host_fasta,
    params.host_bwa,
    params.seq_sim_ref_dir,
    params.seq_sim_config,
    params.primers_fasta,
    params.primers_csv
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LINUX_COMMAND as MERGE_REFS           } from './modules/local/linux/command/main'
include { SEQ_SIMULATOR                         } from './modules/local/seq_simulator/main'
include { SAMPLESHEET_CHECK                     } from './modules/local/samplesheet/check/main'
include { CAT_FASTQ                             } from './modules/nf-core/cat/fastq/main'
include { GFF_FLU                               } from './modules/local/gff_flu/main'
include { ITERATIVE_ALIGNMENT                   } from './modules/local/iterative_alignment/main'
include { MINIMAP2_INDEX                        } from './modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                        } from './modules/nf-core/minimap2/align/main'
include { SAMTOOLS_FAIDX                        } from './modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES                 } from './modules/nf-core/picard/markduplicates/main'
include { ARTIC_ALIGN_TRIM                      } from './modules/local/artic/align_trim/main'
include { SAMTOOLS_INDEX as INDEX_TRIMED        } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_PRIMER_TRIMED } from './modules/nf-core/samtools/index/main'
include { QUAST                                 } from './modules/nf-core/quast/main'
include { SNPEFF_BUILD                          } from './modules/local/snpeff/build/main'
include { SNPEFF_ANN                            } from './modules/local/snpeff/ann/main'
include { MOSDEPTH                              } from './modules/nf-core/mosdepth/main'
include { LINUX_COMMAND as MERGE_CONSENSUS_REF  } from './modules/local/linux/command/main'
include { LINUX_COMMAND as MERGE_CONSENSUS      } from './modules/local/linux/command/main'
include { MUSCLE                                } from './modules/nf-core/muscle/main'
include { PANGOLIN                              } from './modules/nf-core/pangolin/main'
include { NEXTCLADE_DATASETGET                  } from './modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                         } from './modules/nf-core/nextclade/run/main'
include { VCF_REPORT                            } from './modules/local/vcf_report/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from './modules/local/custom_dumpsoftwareversions.nf'
include { MULTIQC                               } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_TRIM_FASTP_FASTQC                         } from './subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { FASTQ_NANOPORE_QC_TRIM                          } from './subworkflows/local/fastq_nanopore_qc_trim/main'
include { ILLUMINA_REMOVE_HOST                            } from './subworkflows/local/illumina_remove_host/main'
include { ASSEMBLE_REFERENCE                              } from './subworkflows/local/assemble_reference/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_VIRAL_SORT_STATS } from './subworkflows/nf-core/bam_sort_stats_samtools/main'
// include { PREPARE_PRIMERS                                 } from './subworkflows/local/prepare_primers/main'
include { NANOPORE_VARCALL                                } from './subworkflows/local/nanopore_varcall/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Init persistant channels
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Init single file channels
    ch_host_fasta     = []
    if(host_fasta) {
        ch_host_fasta = file(host_fasta, checkIfExists: true)
    }
    ch_host_bwa_index = []
    if(host_bwa) {
        ch_host_bwa_index = file(host_bwa, checkIfExists: true)
    }
    ch_viral_gff = []
    if(params.viral_gff) {
        ch_viral_gff = file(params.viral_gff, checkIfExists: true)
    }

    //
    // MODULE: Concat the reference files into one file
    //
    ch_viral_fasta_merge = Channel.fromPath(params.viral_fasta).toSortedList().map{[[id:"viral_reference"], it]}
    .branch {
        meta, fasta ->
            single  : fasta.size() == 1
                return [ meta, fasta.flatten() ]
            multiple: fasta.size() > 1
                return [ meta, fasta.flatten() ]
    }
    MERGE_REFS (
        ch_viral_fasta_merge.multiple,
        [],
        true,
        "merged"
    )
    ch_viral_fasta = MERGE_REFS.out.file.mix(ch_viral_fasta_merge.single)

    //
    // MODULE: Generate fake reads if required
    //
    ch_fastq = Channel.empty()
    if(params.generate_reads) {
        ch_seq_sim_refs   = Channel.from(file(params.seq_sim_ref_dir, checkIfExists: true))
        ch_seq_sim_config = file(params.seq_sim_config, checkIfExists: true)
        SEQ_SIMULATOR (
            ch_seq_sim_refs.map{[ [id: "${params.seq_sim_profile}_test"], it ]},
            ch_seq_sim_config,
            params.seq_sim_profile,
            params.seq_sim_num_reads
        )
        ch_fastq = SEQ_SIMULATOR.out.fastq
    }

    //
    // SECTION: Samplesheet parsing and meta data parsing
    //
    ch_samplesheet   = Channel.empty()
    if (params.samplesheet) {
        //
        // MODULE: Load samplesheet
        //
        ch_samplesheet = file(params.samplesheet, checkIfExists: true)
        SAMPLESHEET_CHECK (
            ch_samplesheet
        )

        //
        // CHANNEL: Construct meta and fastq channel
        //
        ch_fastq = SAMPLESHEET_CHECK.out.csv
        .splitCsv (header:true, sep:",")
        .map {
            it.single_end = true
            def read1 = file(it.read1, checkIfExists: true)
            it.remove("read1")
            it.remove("read2")
            def read2 = null
            if(it.read2) {
                read2 = file(it.read2, checkIfExists: true)
                it.single_end = false
            }

            if (it.read2) {
                [it, [read1, read2]]
            }
            else {
                [it, [read1]]
            }
        }
    }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    ch_fastq_merge = ch_fastq
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    CAT_FASTQ (
        ch_fastq_merge.multiple
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    ch_fastq    = CAT_FASTQ.out.reads.mix(ch_fastq_merge.single)

    //
    // SECTION: Read QC and preprocessing
    //
    if(params.run_illumina_qc_trim) {
        //
        // SUBWORKFLOW: Fastqc and trimming
        //
        ch_trim_primers  = []
        if (params.trim_primers_from_reads && params.primers_fasta) {
            ch_trim_primers = Channel.from(file(params.primers_fasta, checkIfExists: true)).collect()
        }
        FASTQ_TRIM_FASTP_FASTQC (
            ch_fastq,        // ch_reads
            ch_trim_primers, // ch_adapter_fasta
            false,           // val_save_trimmed_fail
            false,           // val_save_merged
            false,           // val_skip_fastp
            false,           // val_skip_fastqc
        )
        ch_versions      = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.trim_json.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip.collect{it[1]})
        ch_fastq         = FASTQ_TRIM_FASTP_FASTQC.out.reads
    }
    if(params.run_nanopore_qc_trim) {
        //
        // SUBWORKFLOW: Nanopore qc and trimming
        //
        FASTQ_NANOPORE_QC_TRIM (
            ch_fastq
        )
        ch_versions = ch_versions.mix(FASTQ_NANOPORE_QC_TRIM.out.versions)
    }

    //
    // SUBWORKFLOW: Remove host reads
    //
    if(params.remove_host_reads) {
        ILLUMINA_REMOVE_HOST (
            ch_fastq,
            ch_host_fasta,
            ch_host_bwa_index
        )
        ch_versions      = ch_versions.mix(ILLUMINA_REMOVE_HOST.out.versions)
        ch_fastq         = ILLUMINA_REMOVE_HOST.out.viral_fastq
        ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA_REMOVE_HOST.out.host_bam_stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA_REMOVE_HOST.out.host_bam_flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA_REMOVE_HOST.out.host_bam_idxstats.collect{it[1]})
    }

    //
    // SECTION: Assemble reference or assign and index
    //
    ch_viral_ref   = Channel.empty()
    ch_fastq_fasta = Channel.empty()
    if(params.assemble_ref) {
        //
        // SUBWORKFLOW: Assemble reference from a list of possible references in the viral_fasta
        //
        ASSEMBLE_REFERENCE (
            ch_fastq,
            ch_viral_fasta
        )
        ch_versions    = ch_versions.mix(ASSEMBLE_REFERENCE.out.versions)
        ch_viral_ref   = ASSEMBLE_REFERENCE.out.viral_ref
        ch_fastq_fasta = ASSEMBLE_REFERENCE.out.fastq_fasta
    } else {
        ch_viral_ref = ch_viral_fasta
    }

    //
    // MODULE: Index ref
    //
    SAMTOOLS_FAIDX (
        ch_viral_ref,
        [[],[]]
    )
    ch_versions      = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_viral_ref_fai = SAMTOOLS_FAIDX.out.fai

    //
    // CHANNEL: Join ref to fai
    //
    ch_viral_ref_fasta_fai = ch_viral_ref
    .map { [it[0].id, it ]}
    .join ( ch_viral_ref_fai.map { [it[0].id, it[1]] })
    .map{ [it[1][0], it[1][1], it[2]] }

    //
    // SECTION: Annotate ref if able and required
    //

    //
    // MODULE: Annotate flu ref
    //
    if(!params.viral_gff && params.annotate_flu_ref) {
        GFF_FLU (
            ch_viral_ref
        )
    }

    //
    // SECTION: Alignment
    //
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    if(params.run_iterative_align) {
        //
        // MODULE: Run iterative alignment
        //
        ITERATIVE_ALIGNMENT (
            ch_fastq_fasta
        )
        ch_bam            = ITERATIVE_ALIGNMENT.out.bam
        ch_bai            = ITERATIVE_ALIGNMENT.out.bai
        ch_consensus_wref = ITERATIVE_ALIGNMENT.out.consensus_wref
        ch_consensus_wn   = ITERATIVE_ALIGNMENT.out.consensus_wn
        ch_final_ref      = ITERATIVE_ALIGNMENT.out.final_ref
        ch_align_metrics  = ITERATIVE_ALIGNMENT.out.metrics
    }
    else if(params.run_bwa_align) {

    }
    else if(params.run_minimap_align) {
        //
        // MODULE: Minimap index
        //
        MINIMAP2_INDEX (
            ch_viral_ref
        )
        ch_versions    = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_mm2_index   = MINIMAP2_INDEX.out.index

        //
        // MODULE: Minimap align
        //
        MINIMAP2_ALIGN (
            ch_fastq,
            ch_mm2_index.collect(),
            true,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
        ch_bam      = MINIMAP2_ALIGN.out.bam
    }

    //
    // MODULE: Mark duplicates
    //
    if(params.run_illumina_mark_dups) {
        PICARD_MARKDUPLICATES (
            ch_bam,
            [[],[]],
            [[],[]]
        )
        ch_versions      = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})
        ch_bam           = PICARD_MARKDUPLICATES.out.bam
    }

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_VIRAL_SORT_STATS (
        ch_bam,
        [[],[]]
    )
    ch_versions      = ch_versions.mix(BAM_VIRAL_SORT_STATS.out.versions)
    ch_bam           = BAM_VIRAL_SORT_STATS.out.bam
    ch_bai           = BAM_VIRAL_SORT_STATS.out.bai
    ch_multiqc_files = ch_multiqc_files.mix(BAM_VIRAL_SORT_STATS.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_VIRAL_SORT_STATS.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_VIRAL_SORT_STATS.out.idxstats.collect{it[1]})

    //
    // SECTION: Primer pre-processing
    //
    ch_primer_bed = Channel.empty()
    if(params.primers_fasta && params.primers_csv) {
        // PREPARE_PRIMERS (
        //     ch_viral_ref,
        //     file(params.primers_fasta),
        //     file(params.primers_csv)
        // )
        // ch_versions = ch_versions.mix(PREPARE_PRIMERS.out.versions)
    } else if(params.primers_bed) {
        ch_primer_bed = Channel.from(file(params.primers_bed, checkIfExists: true)).collect()
    }

    //
    // SECTION: Primer trimming
    //
    ch_trimmed_bam            = Channel.empty()
    ch_primer_trimmed_bam     = Channel.empty()
    ch_trimmed_bai        = Channel.empty()
    ch_primer_trimmed_bai = Channel.empty()
    if(params.run_artic_primer_trim) {
        //
        // MODULE: Trim primers from reads and assig read group to primer pool
        //
        ARTIC_ALIGN_TRIM (
            ch_bam,
            ch_primer_bed
        )
        ch_trimmed_bam        = ARTIC_ALIGN_TRIM.out.trimmed_bam
        ch_primer_trimmed_bam = ARTIC_ALIGN_TRIM.out.primer_trimmed_bam

        //
        // MODULE: Index the trimmed reads
        //
        INDEX_TRIMED ( ch_trimmed_bam )
        INDEX_PRIMER_TRIMED ( ch_primer_trimmed_bam )
        ch_trimmed_bai = INDEX_TRIMED.out.bai
        ch_primer_trimmed_bai = INDEX_PRIMER_TRIMED.out.bai
    }

    //
    // CHANNEL: Join bam to bai
    //
    ch_trimmed_bam_bai = ch_trimmed_bam
    .map { [it[0].id, it ]}
    .join ( ch_trimmed_bai.map { [it[0].id, it[1]] })
    .map{ [it[1][0], it[1][1], it[2]] }

    ch_primer_trimmed_bam_bai = ch_primer_trimmed_bam
    .map { [it[0].id, it ]}
    .join ( ch_primer_trimmed_bai.map { [it[0].id, it[1]] })
    .map{ [it[1][0], it[1][1], it[2]] }

    //
    // SECTION: Variant and consensus calling
    //
    ch_consensus = Channel.empty()
    ch_variants  = Channel.empty()
    ch_vcf_files = Channel.empty()
    if(params.run_nanopore_varcall) {
        NANOPORE_VARCALL (
            ch_trimmed_bam_bai,
            ch_primer_trimmed_bam_bai,
            ch_primer_bed,
            params.pool_primer_reads,
            ch_viral_ref_fasta_fai,
            ch_viral_gff,
            params.clair3_model,
            params.clair3_platform
        )
        ch_versions  = ch_versions.mix(NANOPORE_VARCALL.out.versions)
        ch_consensus = NANOPORE_VARCALL.out.consensus
        ch_variants  = NANOPORE_VARCALL.out.clair3_vcf_tbi
        ch_vcf_files = NANOPORE_VARCALL.out.vcf_files
    } else if(params.run_illumina_varcall) {

    }

    //
    // SECTION: Post consensus analysis
    //

    //
    // MODULE: Quast assembly QC
    //
    QUAST (
        ch_consensus,
        ch_viral_ref.collect(),
        Channel.of(ch_viral_gff).map{[[], it]}.collect()
    )
    ch_versions      = ch_versions.mix(QUAST.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect{it[1]})

    //
    // MODULE: Variant annotation (only run if annotation provided)
    //
    ch_annotation_vcf = Channel.empty()
    if(ch_viral_gff) {
        SNPEFF_BUILD (
            ch_viral_ref.collect{it[1]},
            ch_viral_gff
        )
        ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions)
        SNPEFF_ANN (
            ch_variants.map{[it[0], it[1]]},
            SNPEFF_BUILD.out.db.collect(),
            SNPEFF_BUILD.out.config.collect(),
            ch_viral_ref.collect{it[1]}
        )
        ch_versions      = ch_versions.mix(SNPEFF_ANN.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_ANN.out.csv.collect{it[1]})
        ch_vcf_files     = ch_vcf_files.mix(SNPEFF_ANN.out.vcf.map{[it[0], it[1], "snpeff", 4]})
    }

    //
    // CHANNEL: Join bam to bai and ref
    //
    ch_bam_bai_fasta_fai = Channel.empty()
    if(params.assemble_ref) {
        ch_bam_bai_fasta_fai = ch_primer_trimmed_bam_bai
            .map { [it[0].id, it ]}
            .join ( ch_viral_ref_fasta_fai.map { [it[0].id, it[1], it[2]] })
            .map{ [it[1][0], it[1][1], it[1][2], it[3], it[4]] }
    } else {
        ch_bam_bai_fasta_fai = ch_primer_trimmed_bam_bai
            .combine(ch_viral_ref_fasta_fai)
            .map{ [it[0], it[1], it[2], it[4][0], it[5]] }
    }

    //
    // MODULE: Genome-wide coverage
    //
    MOSDEPTH (
        ch_bam_bai_fasta_fai.map{[it[0], it[1], it[2], []]},
        ch_bam_bai_fasta_fai.map{[it[0], it[3]]},
    )
    ch_versions      = ch_versions.mix(MOSDEPTH.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})

    //
    // MODULE: Merge ref and conesensus seq into one file
    //
    ch_consensus_fasta_merge_ref = ch_viral_ref
        .map{it[1]}
        .mix(ch_consensus.map{it[1]})
        .flatten()
        .toSortedList()
        .map{[[id:"consensus"], it]}
    MERGE_CONSENSUS_REF (
        ch_consensus_fasta_merge_ref,
        [],
        true,
        "merged"
    )
    ch_merged_consensus_ref = MERGE_CONSENSUS_REF.out.file

    //
    // MODULE: conesensus seqs into one file
    //
    ch_consensus_fasta_merge = ch_consensus
        .map{it[1]}
        .flatten()
        .toSortedList()
        .map{[[id:"consensus"], it]}
    MERGE_CONSENSUS (
        ch_consensus_fasta_merge,
        [],
        true,
        "merged"
    )
    ch_merged_consensus = MERGE_CONSENSUS_REF.out.file

    //
    // MODULE: MSA
    //
    MUSCLE (
        ch_merged_consensus_ref
    )
    ch_versions = ch_versions.mix(MUSCLE.out.versions)

    //
    // MODULE: Pangolin
    //
    if(params.run_panglolin) {
        PANGOLIN (
            ch_consensus
        )
        ch_versions = ch_versions.mix(PANGOLIN.out.versions)
    }

    if(params.run_nextclade && params.nextclade_dataset_name) {
        //
        // MODULE: Get nextclade dataset
        //
        NEXTCLADE_DATASETGET (
            params.nextclade_dataset_name,
            params.nextclade_dataset_tag ?: []
        )
        ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

        //
        // MODULE: Run nextclade
        //
        NEXTCLADE_RUN (
            ch_merged_consensus,
            NEXTCLADE_DATASETGET.out.dataset
        )
        ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)
    }

    //
    // CHANNEL: Prepare VCF files for report
    //
    ch_vcf_files = ch_vcf_files
        .groupTuple(by: [0])
        .map { meta, files, callers, order ->
                def sorted_files_and_callers = [files, callers].transpose().sort { a, b ->
                order[callers.indexOf(a[1])] <=> order[callers.indexOf(b[1])]
            }.transpose()
            [meta, sorted_files_and_callers[0], sorted_files_and_callers[1]]
        }

    //
    // MODULE: Generate VCF report
    //
    VCF_REPORT (
        ch_vcf_files
    )

    //
    // MODULE: Track software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile()
    )

    //
    // MODULE: MULTIQC
    //
    workflow_summary = multiqc_summary(workflow, params)
    ch_workflow_summary = Channel.value(workflow_summary)
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        [],
        ch_multiqc_logo,
        [],
        []
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
