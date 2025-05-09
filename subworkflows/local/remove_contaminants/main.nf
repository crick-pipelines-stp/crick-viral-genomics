//
// Remove contaminants from reads
//

include { LINUX_COMMAND as FORCE_REF_UPPER                 } from '../../../modules/local/linux/command/main'
include { LINUX_COMMAND as EXTRACT_REF_IDS                 } from '../../../modules/local/linux/command/main'
include { LINUX_COMMAND as MERGE_REFS                      } from '../../../modules/local/linux/command/main'
include { BWA_INDEX as BWA_INDEX_HOST                      } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_ALIGN_HOST                        } from '../../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_INDEX                                   } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                                   } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CONTAM          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILT_SEC          } from '../../../modules/nf-core/samtools/view/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_CONTAM_SORT_STATS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_CONTAM            } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_VIRAL             } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_VIRAL           } from '../../../modules/nf-core/samtools/fastq/main'

workflow REMOVE_CONTAMINANTS {
    take:
    ch_fastq          // channel: [ val(meta), path(reads) ]
    viral_fasta       // file: [ path(fasta) ]
    contaminant_fasta // file: [ path(fasta) ]
    mode              // val: ont/illumina
    multi_ref         // val: true/false

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Force reference to upper case
    //
    FORCE_REF_UPPER (
        contaminant_fasta,
        [],
        false,
        "upper"
    )
    ch_contaminant_fasta = FORCE_REF_UPPER.out.file

    //
    // MODULE: Extract contig names from fasta files
    //
    EXTRACT_REF_IDS (
        viral_fasta.map{it[1]}.collect().map{[[id:"viral"], it]},
        [],
        false,
        "contigs"
    )
    ch_viral_ref_ids = EXTRACT_REF_IDS.out.file

    //
    // CHANNEL: Merge refs into one fasta but respectful of different pipeline modes
    //
    ch_merged_refs = viral_fasta
        .combine(ch_contaminant_fasta.map{ it[1] })
        .map{ [it[0], [it[1], it[2] ]] }

    //
    // MODULE: Merge refs into one fasta
    //
    MERGE_REFS (
        ch_merged_refs,
        [],
        true,
        "contam_merged"
    )
    ch_merged_refs = MERGE_REFS.out.file

    //
    // MODULE: Generate index
    //
    ch_contam_index = Channel.empty()
    ch_contam_bam   = Channel.empty()
    if (mode == "illumina") {
        // BWA_INDEX_HOST (ch_merged_refs )
        // ch_versions   = ch_versions.mix(BWA_INDEX_HOST.out.versions)
        // ch_contam_index = BWA_INDEX_HOST.out.index

            //     //
    //     // MODULE: Map raw reads to host
    //     //
    //     BWA_ALIGN_HOST (
    //         ch_fastq,
    //         ch_host_bwa_index.collect(),
    //         ch_host_fasta.collect(),
    //         "view"
    //     )
    //     ch_versions = ch_versions.mix(BWA_ALIGN_HOST.out.versions)
    //     ch_host_bam = BWA_ALIGN_HOST.out.bam
    } else if(mode == "ont") {
        //
        // MODULE: Build minimap index
        //
        MINIMAP2_INDEX ( ch_merged_refs )
        ch_versions     = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_contam_index = MINIMAP2_INDEX.out.index

        if(multi_ref) {
            ch_fastq_mm2_index = ch_fastq
                .map { [it[0].id, it ]}
                .join ( ch_contam_index.map { [it[0].id, it[1]] })
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
            ch_versions   = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
            ch_contam_bam = MINIMAP2_ALIGN.out.bam
        }
        else {
            //
            // MODULE: Minimap align
            //
            MINIMAP2_ALIGN (
                ch_fastq,
                ch_contam_index.collect(),
                true,
                false,
                false,
                false
            )
            ch_versions   = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
            ch_contam_bam = MINIMAP2_ALIGN.out.bam
        }
    }

    //
    // MODULE: Index reads
    //
    SAMTOOLS_INDEX_CONTAM (
        ch_contam_bam,
    )
    ch_versions   = ch_versions.mix(SAMTOOLS_INDEX_CONTAM.out.versions)
    ch_contam_bai = SAMTOOLS_INDEX_CONTAM.out.bai
    ch_contam_bam_bai = ch_contam_bam
    .map { row -> [row[0].id, row ].flatten()}
    .join ( ch_contam_bai.map { row -> [row[0].id, row ].flatten()} )
    .map { row -> [row[1], row[2], row[4]] }

    //
    // MODULE: Filter out secondary alignments
    //
    SAMTOOLS_VIEW_FILT_SEC (
        ch_contam_bam_bai,
        [[],[]],
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FILT_SEC.out.versions)

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_CONTAM_SORT_STATS (
        SAMTOOLS_VIEW_FILT_SEC.out.bam,
        [[],[]]
    )
    ch_versions            = ch_versions.mix(BAM_CONTAM_SORT_STATS.out.versions)
    ch_contam_bam          = BAM_CONTAM_SORT_STATS.out.bam
    ch_contam_bai          = BAM_CONTAM_SORT_STATS.out.bai
    ch_contam_bam_stats    = BAM_CONTAM_SORT_STATS.out.stats
    ch_contam_bam_flagstat = BAM_CONTAM_SORT_STATS.out.flagstat
    ch_contam_bam_idxstats = BAM_CONTAM_SORT_STATS.out.idxstats

    //
    // CHANNEL: Combine bam and bai files
    //
    ch_contam_bam_bai = ch_contam_bam
    .map { row -> [row[0].id, row ].flatten()}
    .join ( ch_contam_bai.map { row -> [row[0].id, row ].flatten()} )
    .map { row -> [row[1], row[2], row[4]] }

    //
    // MODULE: Filter out contamination
    //
    SAMTOOLS_VIEW_CONTAM (
        ch_contam_bam_bai,
        [[],[]],
        [],
        ch_viral_ref_ids.map{it[1]}
    )
    ch_versions  = ch_versions.mix(SAMTOOLS_VIEW_CONTAM.out.versions)
    ch_viral_bam = SAMTOOLS_VIEW_CONTAM.out.bam

    //
    // MODULE: Sort by name
    //
    SAMTOOLS_SORT_VIRAL (
        ch_viral_bam,
        [[],[]]
    )
    ch_versions  = ch_versions.mix(SAMTOOLS_SORT_VIRAL.out.versions)
    ch_viral_bam = SAMTOOLS_SORT_VIRAL.out.bam

    //
    // MODULE: Convert to fastq
    //
    SAMTOOLS_FASTQ_VIRAL (
        ch_viral_bam,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_VIRAL.out.versions)

    // Reads come out of different channels depending if we have paried end or single-end
    if(mode == "illumina") {
        ch_viral_fastq = SAMTOOLS_FASTQ_VIRAL.out.fastq
    } else {
        ch_viral_fastq = SAMTOOLS_FASTQ_VIRAL.out.other
    }

    emit:
    viral_fastq         = ch_viral_fastq
    contam_bam          = ch_contam_bam
    contam_bai          = ch_contam_bai
    contam_bam_stats    = ch_contam_bam_stats
    contam_bam_flagstat = ch_contam_bam_flagstat
    contam_bam_idxstats = ch_contam_bam_idxstats
    versions            = ch_versions
}
