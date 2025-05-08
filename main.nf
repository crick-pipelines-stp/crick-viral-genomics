#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { params_summary_map      } from './modules/local/util/logging/main'
include { params_summary_map_json } from './modules/local/util/logging/main'
include { summary_log             } from './modules/local/util/logging/main'
include { multiqc_summary         } from './modules/local/util/logging/main'
include { get_genome_attribute    } from './modules/local/util/references/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW INIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create workflow summary
log.info summary_log(workflow, params, params.debug, params.monochrome_logs)
def summary_params = params_summary_map(workflow, params, params.debug)
def json_summary = params_summary_map_json(workflow, params, false)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LINUX_COMMAND as MERGE_REFS            } from './modules/local/linux/command/main'
include { SEQ_SIMULATOR                          } from './modules/local/seq_simulator/main'
include { SAMPLESHEET_CHECK                      } from './modules/local/samplesheet/check/main'
include { SNAPGENE_TO_GFF                        } from './modules/local/snapgene_to_gff/main'
include { LINUX_COMMAND as FORCE_REF_UPPER       } from './modules/local/linux/command/main'
include { CAT_FASTQ                              } from './modules/nf-core/cat/fastq/main'
include { LINUX_COMMAND as SPLIT_REF             } from './modules/local/linux/command/main'
include { GFF_FLU                                } from './modules/local/gff_flu/main'
include { LINUX_COMMAND as MERGE_GFF             } from './modules/local/linux/command/main'
include { ITERATIVE_ALIGNMENT                    } from './modules/local/iterative_alignment/main'
include { MINIMAP2_INDEX                         } from './modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                         } from './modules/nf-core/minimap2/align/main'
include { BWA_INDEX as BWA_INDEX_VIRUS           } from './modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_ALIGN_VIRUS             } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FAIDX                         } from './modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES                  } from './modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VIRUS } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_VIRUS } from './modules/nf-core/samtools/view/main'
include { SEQKIT_BAM                             } from './modules/local/seqkit/bam/main'
include { ARTIC_ALIGN_TRIM                       } from './modules/local/artic/align_trim/main'
include { SAMTOOLS_SORT as SORT_PRIMER_TRIMMED   } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as INDEX_TRIMMED        } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_PRIMER_TRIMMED } from './modules/nf-core/samtools/index/main'
include { IVAR_TRIM                              } from './modules/nf-core/ivar/trim/main'
include { QUAST                                  } from './modules/nf-core/quast/main'
include { SAMTOOLS_MPILEUP                       } from './modules/nf-core/samtools/mpileup/main'
include { GEN_COUNT_TABLE                        } from './modules/local/gen_count_table/main'
include { SNPEFF_BUILD                           } from './modules/local/snpeff/build/main'
include { SNPEFF_ANN                             } from './modules/local/snpeff/ann/main'
include { MOSDEPTH                               } from './modules/nf-core/mosdepth/main'
include { LINUX_COMMAND as MERGE_CONSENSUS_REF   } from './modules/local/linux/command/main'
include { LINUX_COMMAND as MERGE_CONSENSUS       } from './modules/local/linux/command/main'
include { MUSCLE                                 } from './modules/nf-core/muscle/main'
include { PANGOLIN                               } from './modules/nf-core/pangolin/main'
include { NEXTCLADE_DATASETGET                   } from './modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                          } from './modules/nf-core/nextclade/run/main'
include { EXPORT_REPORT_DATA                     } from './modules/local/export_report_data/main'
include { VCF_REPORT                             } from './modules/local/vcf_report/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from './modules/local/custom_dumpsoftwareversions.nf'
include { MULTIQC                                } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_TRIM_FASTP_FASTQC                         } from './subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { FASTQ_NANOPORE_QC_TRIM                          } from './subworkflows/local/fastq_nanopore_qc_trim/main'
include { REMOVE_HOST                                     } from './subworkflows/local/remove_host/main'
include { REMOVE_CONTAMINANTS                             } from './subworkflows/local/remove_contaminants/main'
include { ASSEMBLE_REFERENCE                              } from './subworkflows/local/assemble_reference/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_VIRAL_SORT_STATS } from './subworkflows/nf-core/bam_sort_stats_samtools/main'
include { PREPARE_PRIMERS                                 } from './subworkflows/local/prepare_primers/main'
include { NANOPORE_VARCALL                                } from './subworkflows/local/nanopore_varcall/main'
include { ILLUMINA_VARCALL                                } from './subworkflows/local/illumina_varcall/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //////////////////////////////////////
    // INIT:
    //////////////////////////////////////

    // Load static config
    ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_logo   = file("$projectDir/assets/The_Francis_Crick_Institute_logo.png", checkIfExists: true)
    ch_seq_sim_config = file(params.seq_sim_config, checkIfExists: true)

    // Configure run_id
    def run_id = params.run_id
    if(params.run_id == null) {
        run_id = "viral_genomics_pipeline"
    }

    // Resolve references, load from genome if possible but then override from params if supplied
    def host_fasta = get_genome_attribute(params, 'fasta')
    def host_bwa   = get_genome_attribute(params, 'bwa'  )
    if(params.host_fasta) {
        host_fasta = params.host_fasta
    }
    if(params.host_bwa) {
        host_bwa = params.host_bwa
    }

    // Init persistant channels
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //////////////////////////////////////
    // CHECK PARAMS:
    //////////////////////////////////////

    // If using ref per sample, load from samplesheet, otherewise we need a ref dir or single fasta
    if(!params.use_independant_refs && params.viral_fasta == null) {
        exit 1, "Required parameter not specified: viral_fasta"
    }

    // If no data being generated, samplesheet is manditory
    if(!params.run_generate_reads && params.samplesheet == null) {
        exit 1, "Required parameter not specified: samplesheet"
    }

    // Check non-manditory input parameters to see if the files exist if they have been specified
    check_param_list = [
        params.samplesheet,
        params.viral_gff,
        params.host_fasta,
        params.host_bwa,
        params.contaminant_fasta,
        params.seq_sim_ref_dir,
        params.seq_sim_config,
        params.primers_bed,
        params.primers_fasta,
        params.primers_csv
    ]
    check_param_list.each { param -> if (param) { file(param, checkIfExists: true) } }

    // Init variables
    def multi_ref = params.run_assemble_ref || params.use_independant_refs || params.run_iterative_align
    def is_gff    = params.viral_gff != null || params.annotate_flu_ref || params.use_independant_gff

    // Init single file channels
    ch_host_fasta = []
    if(host_fasta) {
        ch_host_fasta = file(host_fasta, checkIfExists: true)
    }
    ch_host_bwa_index = []
    if(host_bwa) {
        ch_host_bwa_index = file(host_bwa, checkIfExists: true)
    }
    ch_contaminent_fasta = []
    if(params.contaminant_fasta) {
        ch_contaminent_fasta = Channel.of(file(params.contaminant_fasta, checkIfExists: true)).map{[[id:"contaminent"], it]}.collect()
    }
    ch_viral_gff = [[],[]]
    if(params.viral_gff) {
        // ch_viral_gff = file(params.viral_gff, checkIfExists: true)
        ch_viral_gff = Channel.of(file(params.viral_gff, checkIfExists: true)).map{[[id:"gff"], it]}
    }

    //////////////////////////////////////
    // MAIN WORKFLOW:
    //////////////////////////////////////

    //
    // MODULE: Concat the reference files into one file
    //
    ch_viral_fasta = Channel.empty()
    if(!params.use_independant_refs) {
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
        ch_viral_fasta = ch_viral_fasta.map{ [ it[0], (it[1] instanceof List ? it[1][0] : it[1] ) ] }
    }

    //
    // MODULE: Generate fake reads if required
    //
    ch_fastq = Channel.empty()
    if(params.run_generate_reads) {
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
            if (!it.containsKey("id")) {
                exit 1, "The header 'id' is not in the samplesheet"
            }
            if (!it.containsKey("read1")) {
                exit 1, "The header 'read1' is not in the samplesheet"
            }

            it.single_end = true
            def read1 = file(it.read1, checkIfExists: true)
            it.remove("read1")
            def read2 = null
            if(it.read2) {
                read2 = file(it.read2, checkIfExists: true)
                it.single_end = false
                it.remove("read2")
            }
            if (read2) {
                [it, [read1, read2]]
            }
            else {
                it.remove("read2")
                [it, [read1]]
            }
        }

        //
        // CHANNEL: Check for folders and expand
        //
        ch_fastq_folder = ch_fastq.branch {
            meta, fastq ->
                single  : !fastq[0].isDirectory()
                    return [ meta, fastq ]
                folder: fastq[0].isDirectory()
                    return [ meta, fastq[0] ]
        }
        ch_fastq_folder_expanded = ch_fastq_folder.folder
            .flatMap {
                meta, folder ->
                    def files = folder
                        .listFiles()
                        .findAll { it.name ==~ /.*\.(fastq|fq)(\.gz)?/ }
                        .sort { it.name }
                    files.collect { file ->
                        [ meta, file ]
                    }
            }
        ch_fastq = ch_fastq_folder.single.mix(ch_fastq_folder_expanded)
    }

    //
    // MODULE: Independant reference and gff processing
    //
    if(params.use_independant_refs) {
        ch_viral_fasta = ch_fastq
            .map{
                def ref = file(it[0].ref, checkIfExists: true)
                [[id:it[0].id], ref]
            }
    }
    if(params.use_independant_gff) {
        // Load GFF
        ch_viral_gff = ch_fastq
            .map{
                def gff = file(it[0].gff, checkIfExists: true)
                [[id:it[0].id], [gff]]
            }
    }

    if(params.convert_snapgene) {
        // Extract the ref name from the fasta file
        ch_viral_fasta_ref_name = ch_viral_fasta
            .map { meta, fasta ->
                def firstLine = fasta.withReader { reader -> reader.readLine() }
                firstLine = firstLine.replaceFirst('^>', '')
                [meta, firstLine]
            }

        // Join with gff for contig naming
        if(params.use_independant_refs) {
            ch_viral_gff = ch_viral_gff.join(ch_viral_fasta_ref_name)
        }
        else {
            ch_viral_gff = ch_viral_gff
                .combine(ch_viral_fasta_ref_name)
                .map{ [it[0], it[1], it[3]]}
        }

        //
        // MODULE: Convert snapgene to gff if required
        //
        SNAPGENE_TO_GFF (
            ch_viral_gff
        )
        ch_viral_gff = SNAPGENE_TO_GFF.out.gff
    }

    //
    // MODULE: Force ref upper
    //
    if(params.force_ref_to_upper) {
        FORCE_REF_UPPER (
            ch_viral_fasta,
            [],
            false,
            "upper"
        )
        ch_viral_fasta = FORCE_REF_UPPER.out.file
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
    ch_versions   = ch_versions.mix(CAT_FASTQ.out.versions)
    ch_fastq      = CAT_FASTQ.out.reads.mix(ch_fastq_merge.single)
    ch_orig_fastq = []
    if(params.run_nanopore_qc_trim) {
        ch_orig_fastq = ch_fastq.map{it[1]}.collect()
    }

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
    ch_report_data_host = []
    if(params.run_remove_host_reads && host_fasta != null) {
        def remove_host_mode = "illumina"
        if(params.run_minimap_align) {
            remove_host_mode = "ont"
        }
        REMOVE_HOST (
            ch_fastq,
            ch_host_fasta,
            ch_host_bwa_index,
            remove_host_mode

        )
        ch_versions         = ch_versions.mix(REMOVE_HOST.out.versions)
        ch_fastq            = REMOVE_HOST.out.viral_fastq
        ch_multiqc_files    = ch_multiqc_files.mix(REMOVE_HOST.out.host_bam_stats.collect{it[1]})
        ch_multiqc_files    = ch_multiqc_files.mix(REMOVE_HOST.out.host_bam_flagstat.collect{it[1]})
        ch_multiqc_files    = ch_multiqc_files.mix(REMOVE_HOST.out.host_bam_idxstats.collect{it[1]})
        ch_report_data_host = REMOVE_HOST.out.host_bam_flagstat.collect{it[1]}
    }

    //
    // SECTION: Assemble reference or assign and index
    //
    ch_viral_ref   = Channel.empty()
    ch_fastq_fasta = Channel.empty()
    if(params.run_assemble_ref) {
        //
        // SUBWORKFLOW: Assemble reference from a list of possible references in the viral_fasta
        //
        def assemble_mode = "illumina"
        if(params.run_minimap_align) {
            assemble_mode = "ont"
        }
        ASSEMBLE_REFERENCE (
            ch_fastq,
            ch_viral_fasta,
            assemble_mode
        )
        ch_versions    = ch_versions.mix(ASSEMBLE_REFERENCE.out.versions)
        ch_viral_ref   = ASSEMBLE_REFERENCE.out.viral_ref
        ch_fastq_fasta = ASSEMBLE_REFERENCE.out.fastq_fasta
    } else {
        ch_viral_ref   = ch_viral_fasta
        ch_fastq_fasta = ch_fastq.combine(ch_viral_ref.map{it[1]})
    }

    //
    // SUBWORKFLOW: Remove contaminent reads
    //
    ch_report_data_contam = []
    if(params.run_remove_contaminant_reads && params.contaminant_fasta != null) {
        def remove_contaminent_mode = "illumina"
        if(params.run_minimap_align) {
            remove_contaminent_mode = "ont"
        }
        REMOVE_CONTAMINANTS (
            ch_fastq,
            ch_viral_ref,
            ch_contaminent_fasta,
            remove_contaminent_mode,
            multi_ref

        )
        ch_versions           = ch_versions.mix(REMOVE_CONTAMINANTS.out.versions)
        ch_fastq              = REMOVE_CONTAMINANTS.out.viral_fastq
        ch_multiqc_files      = ch_multiqc_files.mix(REMOVE_CONTAMINANTS.out.contam_bam_stats.collect{it[1]})
        ch_multiqc_files      = ch_multiqc_files.mix(REMOVE_CONTAMINANTS.out.contam_bam_flagstat.collect{it[1]})
        ch_multiqc_files      = ch_multiqc_files.mix(REMOVE_CONTAMINANTS.out.contam_bam_idxstats.collect{it[1]})
        ch_report_data_contam = REMOVE_CONTAMINANTS.out.contam_bam_idxstats.collect{it[1]}
    }

    //
    // SECTION: Annotate ref if able and required
    //

    // Annotate flu ref
    if(!params.viral_gff && params.annotate_flu_ref) {
        //
        // MODULE: Split ref into chunks
        //
        SPLIT_REF (
            ch_viral_ref,
            [],
            false,
            "split"
        )
        ch_split_refs = SPLIT_REF.out.file
        .flatMap { meta, files ->
            files.collect { file -> [meta, file] }
        }

        //
        // MODULE: Annotate chunks and group back into samples
        //
        GFF_FLU (
            ch_split_refs
        )
        ch_grouped_gff = GFF_FLU.out.snpeff_gff
            .groupTuple(by: [0])

        //
        // MODULE: Merge GFF
        //
        MERGE_GFF (
            ch_grouped_gff,
            [],
            false,
            "merged"
        )
        ch_viral_gff = MERGE_GFF.out.file
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
        ch_iter_final_ref = ITERATIVE_ALIGNMENT.out.final_ref

        //
        // CHANNEL: Assign new consensus as ref
        //
        ch_viral_ref = ch_iter_final_ref
    }
    else if(params.run_bwa_align) {
        //
        // MODULE: BWA index
        //
        BWA_INDEX_VIRUS (
            ch_viral_ref
        )
        ch_versions  = ch_versions.mix(BWA_INDEX_VIRUS.out.versions)
        ch_bwa_index = BWA_INDEX_VIRUS.out.index

        if(multi_ref) {
            //
            // CHANNEL: BWA align prep for multiple refs
            //
            ch_fastq_idx_ref = ch_fastq
                .map { [it[0].id, it ]}
                .join ( ch_bwa_index.map { [it[0].id, it[1]] })
                .join ( ch_viral_ref.map { [it[0].id, it[1]] })
                .map{ [it[1][0], it[1][1], it[2], it[3]] }

            //
            // MODULE: BWA align
            //
            BWA_ALIGN_VIRUS (
                ch_fastq_idx_ref.map{[it[0], it[1]]},
                ch_fastq_idx_ref.map{[it[0], it[2]]},
                ch_fastq_idx_ref.map{[it[0], it[3]]},
                true
            )
            ch_versions = ch_versions.mix(BWA_ALIGN_VIRUS.out.versions)
            ch_bam      = BWA_ALIGN_VIRUS.out.bam
        }
        else {
            //
            // MODULE: BWA align
            //
            BWA_ALIGN_VIRUS (
                ch_fastq,
                ch_bwa_index.collect(),
                ch_viral_ref.collect(),
                true
            )
            ch_versions = ch_versions.mix(BWA_ALIGN_VIRUS.out.versions)
            ch_bam      = BWA_ALIGN_VIRUS.out.bam
        }
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

        if(multi_ref) {
            ch_fastq_mm2_index = ch_fastq
                .map { [it[0].id, it ]}
                .join ( ch_mm2_index.map { [it[0].id, it[1]] })
                .map{ [it[1][0], it[1][1], it[2]] }

            //
            // MODULE: Minimap align
            //
            MINIMAP2_ALIGN (
                ch_fastq_mm2_index.map{[it[0], it[1]]},
                ch_fastq_mm2_index.map{[it[0], it[2]]},
                true,
                false,
                false,
                false
            )
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
            ch_bam      = MINIMAP2_ALIGN.out.bam
        }
        else {
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
    // MODULE: Index and Filter for primary alignments only
    //
    SAMTOOLS_INDEX_VIRUS (
        ch_bam,
    )
    ch_versions   = ch_versions.mix(SAMTOOLS_INDEX_VIRUS.out.versions)
    ch_bai        = SAMTOOLS_INDEX_VIRUS.out.bai
    ch_bam_bai = ch_bam
        .map { [it[0].id, it ]}
        .join ( ch_bai.map { [it[0].id, it[1]] })
        .map{ [it[1][0], it[1][1], it[2]] }
    SAMTOOLS_FILTER_VIRUS (
        ch_bam_bai,
        [[],[]],
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FILTER_VIRUS.out.versions)
    ch_bam      = SAMTOOLS_FILTER_VIRUS.out.bam

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
    ch_versions          = ch_versions.mix(BAM_VIRAL_SORT_STATS.out.versions)
    ch_bam               = BAM_VIRAL_SORT_STATS.out.bam
    ch_bai               = BAM_VIRAL_SORT_STATS.out.bai
    ch_multiqc_files     = ch_multiqc_files.mix(BAM_VIRAL_SORT_STATS.out.stats.collect{it[1]})
    ch_multiqc_files     = ch_multiqc_files.mix(BAM_VIRAL_SORT_STATS.out.flagstat.collect{it[1]})
    ch_multiqc_files     = ch_multiqc_files.mix(BAM_VIRAL_SORT_STATS.out.idxstats.collect{it[1]})
    ch_report_data_align = BAM_VIRAL_SORT_STATS.out.flagstat.collect{it[1]}

    //
    // MODULE: Calculate Seqkit stats
    //
    SEQKIT_BAM(
        ch_bam
    )
    ch_versions   = ch_versions.mix(SEQKIT_BAM.out.versions)
    ch_seqkit_tsv = SEQKIT_BAM.out.tsv

    //
    // CHANNEL: Join bam to bai
    //
    ch_bam_bai = ch_bam
    .map { [it[0].id, it ]}
    .join ( ch_bai.map { [it[0].id, it[1]] })
    .map{ [it[1][0], it[1][1], it[2]] }

    //
    // SECTION: Primer pre-processing
    //
    ch_primer_bed = Channel.empty()
    if(!params.primers_bed && params.primers_fasta && params.primers_csv) {
        ch_primer_ref = ch_viral_ref
        if(params.run_iterative_align) {
            ch_primer_ref = ch_consensus_wref
        }
        PREPARE_PRIMERS (
            ch_primer_ref,
            file(params.primers_fasta),
            file(params.primers_csv)
        )
        ch_versions   = ch_versions.mix(PREPARE_PRIMERS.out.versions)
        ch_primer_bed = PREPARE_PRIMERS.out.primer_bed
    } else if(params.primers_bed) {
        ch_primer_bed = Channel.from(file(params.primers_bed, checkIfExists: true)).collect()
    }

    //
    // SECTION: Primer trimming
    //
    ch_trimmed_bam        = ch_bam
    ch_primer_trimmed_bam = ch_bam
    ch_trimmed_bai        = ch_bai
    ch_primer_trimmed_bai = ch_bai
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
        INDEX_TRIMMED ( ch_trimmed_bam )
        INDEX_PRIMER_TRIMMED ( ch_primer_trimmed_bam )
        ch_trimmed_bai = INDEX_TRIMMED.out.bai
        ch_primer_trimmed_bai = INDEX_PRIMER_TRIMMED.out.bai
    }
    else if(params.run_ivar_primer_trim) {
        //
        // CHANNEL: merge bam/bai with primer beds if we did iterative align
        //
        ch_ivar_bam_bai = ch_bam_bai
        ch_ivar_primer_bed = ch_primer_bed
        if(params.run_iterative_align) {
            ch_ivar_input = ch_bam_bai
                .map{[it[0].id, it]}
                .join(ch_primer_bed.map{[it[0].id, it[1]]})
                .map{ [it[1][0], it[1][1], it[1][2], it[2]] }
            ch_ivar_bam_bai = ch_ivar_input.map{[it[0], it[1], it[2]]}
            ch_ivar_primer_bed = ch_ivar_input.map{[it[3]]}
        }

        //
        // MODULE: Trim primers using ivar
        //
        IVAR_TRIM (
            ch_ivar_bam_bai,
            ch_ivar_primer_bed
        )
        ch_versions           = ch_versions.mix(IVAR_TRIM.out.versions)
        ch_trimmed_bam        = IVAR_TRIM.out.bam
        ch_primer_trimmed_bam = IVAR_TRIM.out.bam

        //
        // MODULE: Sort and Index the trimmed reads
        //
        SORT_PRIMER_TRIMMED ( ch_trimmed_bam, [[],[]] )
        INDEX_PRIMER_TRIMMED ( SORT_PRIMER_TRIMMED.out.bam )
        ch_primer_trimmed_bam = SORT_PRIMER_TRIMMED.out.bam
        ch_primer_trimmed_bai = INDEX_PRIMER_TRIMMED.out.bai
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
    // CHANNEL: Join bam to bai and ref
    //
    ch_bam_bai_fasta_fai = Channel.empty()
    if(multi_ref) {
        ch_bam_bai_fasta_fai = ch_primer_trimmed_bam_bai
            .map { [it[0].id, it ]}
            .join ( ch_viral_ref_fasta_fai.map { [it[0].id, it[1], it[2]] })
            .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }
    }
    else {
        ch_bam_bai_fasta_fai = ch_primer_trimmed_bam_bai
            .combine(ch_viral_ref_fasta_fai)
            .map{ [it[0], it[1], it[2], it[4], it[5]] }
    }

    //
    // SECTION: Variant and consensus calling
    //
    ch_consensus = Channel.empty()
    ch_variants  = Channel.empty()
    ch_vcf_files = Channel.empty()
    ch_compressed_vcf = []
    if(params.run_nanopore_varcall) {
        NANOPORE_VARCALL (
            ch_trimmed_bam_bai,
            ch_primer_trimmed_bam_bai,
            ch_primer_bed,
            params.pool_primer_reads,
            ch_viral_ref_fasta_fai,
            params.clair3_model,
            params.clair3_platform,
            multi_ref
        )
        ch_versions  = ch_versions.mix(NANOPORE_VARCALL.out.versions)
        ch_consensus = NANOPORE_VARCALL.out.consensus
        ch_variants  = NANOPORE_VARCALL.out.clair3_vcf_tbi
        ch_vcf_files = NANOPORE_VARCALL.out.vcf_files
        ch_compressed_vcf = NANOPORE_VARCALL.out.medaka_vcf_tbi.map{[it[1], it[2]]}
            .mix(NANOPORE_VARCALL.out.clair3_vcf_tbi.map{[it[1], it[2]]})
            .mix(NANOPORE_VARCALL.out.lofreq_vcf_tbi.map{[it[1], it[2]]})
            .mix(NANOPORE_VARCALL.out.sniffles_vcf_tbi.map{[it[1], it[2]]})
            .flatten()
            .collect()
    } else if(params.run_illumina_varcall) {
        ILLUMINA_VARCALL(
            ch_bam_bai_fasta_fai
        )
        ch_versions  = ch_versions.mix(ILLUMINA_VARCALL.out.versions)
        ch_consensus = ILLUMINA_VARCALL.out.consensus
        ch_variants  = ILLUMINA_VARCALL.out.lofreq_vcf_tbi
        ch_vcf_files = ILLUMINA_VARCALL.out.vcf_files
    }

    //
    // SECTION: Post consensus analysis
    //

    //
    // CHANNEL: Prep channels for QC
    //
    if(!multi_ref) {
        ch_viral_ref = ch_viral_ref.collect()
    }
    else {
        ch_con_ref = ch_consensus
            .map { [it[0].id, it ]}
            .join ( ch_viral_ref.map { [it[0].id, it[1]] })
            .map{ [it[1][0], it[1][1], it[2]] }

        ch_consensus = ch_con_ref.map{[it[0], it[1]]}
        ch_viral_ref = ch_con_ref.map{[it[0], it[2]]}
    }
    if(!multi_ref && params.annotate_flu_ref) {
        ch_viral_gff = ch_viral_gff.collect()
    }
    else if(!multi_ref && is_gff) {
        ch_viral_gff = ch_viral_gff.collect()
    }
    else if (multi_ref && is_gff) {
        ch_con_ref_gff = ch_consensus
            .map { [it[0].id, it ]}
            .join ( ch_viral_ref.map { [it[0].id, it[1]] })
            .join ( ch_viral_gff.map { [it[0].id, it[1]] })
            .map{ [it[1][0], it[1][1], it[2], it[3]] }

        ch_consensus = ch_con_ref_gff.map{[it[0], it[1]]}
        ch_viral_ref = ch_con_ref_gff.map{[it[0], it[2]]}
        ch_viral_gff = ch_con_ref_gff.map{[it[0], it[3]]}
    }

    //
    // MODULE: Quast assembly QC
    //
    if(params.run_illumina_varcall || params.run_nanopore_varcall) {
        QUAST (
            ch_consensus,
            ch_viral_ref,
            ch_viral_gff
        )
        ch_versions      = ch_versions.mix(QUAST.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect{it[1]})
    }

    if(params.run_snpeff && (params.viral_gff || params.annotate_flu_ref || params.use_independant_gff)) {
        //
        // MODULE: Build snpeff db
        //
        SNPEFF_BUILD (
            ch_viral_ref,
            ch_viral_gff
        )
        ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions)

        //
        // CHANNEL: Matchup channels and annotate
        //
        if(!multi_ref) {
            SNPEFF_ANN (
                ch_variants.map{[it[0], it[1]]},
                SNPEFF_BUILD.out.db,
                SNPEFF_BUILD.out.config,
                ch_viral_ref
            )
            ch_versions      = ch_versions.mix(SNPEFF_ANN.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_ANN.out.csv.collect{it[1]})
            ch_vcf_files     = ch_vcf_files.mix(SNPEFF_ANN.out.vcf.map{[it[0], it[1], "snpeff", 4]})
        }
        else {
            ch_var_snpeff_ref = ch_variants
                .map { [it[0].id, it ]}
                .join ( SNPEFF_BUILD.out.db.map { [it[0].id, it[1]] })
                .join ( SNPEFF_BUILD.out.config.map { [it[0].id, it[1]] })
                .join ( ch_viral_ref.map { [it[0].id, it[1]] })
                .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3], it[4]] }

            SNPEFF_ANN (
                ch_var_snpeff_ref.map{[it[0], it[1]]},
                ch_var_snpeff_ref.map{[it[0], it[3]]},
                ch_var_snpeff_ref.map{[it[0], it[4]]},
                ch_var_snpeff_ref.map{[it[0], it[5]]}
            )
            ch_versions      = ch_versions.mix(SNPEFF_ANN.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_ANN.out.csv.collect{it[1]})
            ch_vcf_files     = ch_vcf_files.mix(SNPEFF_ANN.out.vcf.map{[it[0], it[1], "snpeff", 4]})
        }
    }

    ch_count_table = []
    if(params.run_gen_count_table) {
        //
        // MODULE: Calculate the pileip 
        //
        SAMTOOLS_MPILEUP (
            ch_bam_bai_fasta_fai.map{[it[0], it[1], []]},
            ch_bam_bai_fasta_fai.map{[it[3]]}
        )
        ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)

        //
        // MODULE: Generate count table
        //
        GEN_COUNT_TABLE (
            SAMTOOLS_MPILEUP.out.mpileup
        )
        ch_count_table = GEN_COUNT_TABLE.out.csv
    }

    //
    // MODULE: Genome-wide coverage
    //
    MOSDEPTH (
        ch_bam_bai_fasta_fai.map{[it[0], it[1], it[2], []]},
        ch_bam_bai_fasta_fai.map{[it[0], it[3]]},
    )
    ch_versions             = ch_versions.mix(MOSDEPTH.out.versions)
    ch_multiqc_files        = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    ch_report_data_coverage = MOSDEPTH.out.per_base_bed.collect{it[1]}

    if(params.run_nanopore_varcall || params.run_illumina_varcall) {
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
        ch_merged_consensus = MERGE_CONSENSUS.out.file

        //
        // MODULE: MSA
        //
        if(params.run_msa) {
            MUSCLE (
                ch_merged_consensus_ref
            )
            ch_versions = ch_versions.mix(MUSCLE.out.versions)
        }
    }

    //
    // MODULE: Pangolin
    //
    if(params.run_panglolin) {
        PANGOLIN (
            ch_merged_consensus
        )
        ch_versions      = ch_versions.mix(PANGOLIN.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PANGOLIN.out.report.collect{it[1]})
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

    if(params.run_reporting) {
        //
        // CHANNEL: Prepare VCF files for report
        //
        ch_vcf_files_prep = ch_vcf_files
            .groupTuple(by: [0])
            .map { meta, files, callers, order ->
                    def sorted_files_and_callers = [files, callers].transpose().sort { a, b ->
                    order[callers.indexOf(a[1])] <=> order[callers.indexOf(b[1])]
                }.transpose()
                [meta, sorted_files_and_callers[0], sorted_files_and_callers[1]]
            }
            .collect { meta, files, callers ->
                    [callers, files.collect()]
            }
            .map { arr ->
                def callers = arr[0]
                def files = arr.indices.findAll { it % 2 == 1 }.collect { arr[it] }
                [callers, files.flatten()]
            }

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

        if(params.run_report_export) {
            ch_truncation = []
            if(params.calculate_truncation) {
                ch_truncation = ch_seqkit_tsv.collect{it[1]}
            }

            //
            // MODULE: Export report data to pickle file
            //
            EXPORT_REPORT_DATA (
                run_id,
                params.export_report_type,
                json_summary,
                ch_viral_ref_fasta_fai.map{[it[1], it[2]]}.collect(),
                ch_viral_gff.map{it[1]}.collect(),
                ch_samplesheet,
                ch_orig_fastq,
                ch_report_data_host,
                ch_report_data_contam,
                ch_report_data_align,
                ch_report_data_coverage,
                ch_consensus.map{it[1]}.collect(),
                ch_vcf_files_prep,
                ch_compressed_vcf,
                ch_count_table.collect{it[1]},
                ch_truncation
            )
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
