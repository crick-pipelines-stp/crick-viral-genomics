//
// Generate consesnus and call variants from nanopore reads
//

include { MEDAKA_INFERENCE                     } from '../../../modules/local/medaka/inference/main'
include { MEDAKA_VCF                           } from '../../../modules/local/medaka/vcf/main'
include { MEDAKA_ANNOTATE                      } from '../../../modules/local/medaka/annotate/main'
include { ARTIC_VCF_MERGE                      } from '../../../modules/local/artic/vcf_merge/main'
include { BCFTOOLS_PASS_FAIL_SPLIT             } from '../../../modules/local/bcftools/pass_fail_split/main'
include { TABIX_BGZIPTABIX as INDEX_CONSEN_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { ARTIC_MAKE_DEPTH_MASK                } from '../../../modules/local/artic/make_depth_mask/main'
include { ARTIC_MASK                           } from '../../../modules/local/artic/mask/main'
include { BCFTOOLS_CONSENSUS                   } from '../../../modules/nf-core/bcftools/consensus/main'
include { CLAIR3_RUN                           } from '../../../modules/local/clair3/main'
include { BCFTOOLS_VIEW as CLAIR3_FILTER_REF   } from '../../../modules/nf-core/bcftools/view/main'
// include { GUNZIP as GUNZIP_CLAIR3_VCF          } from '../../../modules/nf-core/gunzip/main'
include { LOFREQ_CALL                          } from '../../../modules/local/lofreq/call/main'
include { TABIX_BGZIPTABIX as INDEX_LOFREQ_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { SNIFFLES                             } from '../../../modules/nf-core/sniffles/main'
include { GUNZIP as GUNZIP_SNIFFLES_VCF        } from '../../../modules/nf-core/gunzip/main'
include { LINUX_COMMAND as RENAME_FASTA        } from '../../../modules/local/linux/command/main'

workflow NANOPORE_VARCALL {
    take:
    trimmed_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    primer_trimmed_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    primer_bed             // file
    pool_reads             // val
    reference              // channel: [ val(meta), path(fasta), path(fai) ]
    clair3_model           // val
    clair3_platform        // val
    multi_ref              // val

    main:
    ch_versions = Channel.empty()

    //
    // CHANNEL: join up channels for synced input
    //
    primer_trimmed_bam_bai_ref = Channel.empty()
    if(multi_ref) {
        primer_trimmed_bam_bai_ref = primer_trimmed_bam_bai
            .map { [it[0].id, it ]}
            .join ( reference.map { [it[0].id, it[1], it[2]] })
            .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }
    } else {
        primer_trimmed_bam_bai_ref = primer_trimmed_bam_bai
            .combine(reference.map{[it[1], it[2]]})
            .map{ [it[0], it[1], it[2], it[3], it[4]] }
    }

    //
    // CHANNEL: Get primer pool ids and combine with bam
    //
    ch_pool_trimmed_bam_bai = Channel.empty()
    if(pool_reads) {
        ch_pool_trimmed_bam_bai = primer_bed
        .splitCsv(sep: '\t')
        .map{[it[0][4]]}
        .unique()
        .flatten()
        .combine(trimmed_bam_bai)
        .map{
            def newMeta = it[1].clone()
            newMeta.pool = it[0]
            [it[0], newMeta, it[2], it[3]]
        }
    } else {
        ch_pool_trimmed_bam_bai = trimmed_bam_bai.map{[it[0].id, it[0], it[1], it[2]]}
    }

    //
    // MODULE: Generate inference
    //
    MEDAKA_INFERENCE (
        ch_pool_trimmed_bam_bai.map{[it[1], it[2], it[3]]},
        ch_pool_trimmed_bam_bai.map{[it[0]]}
    )
    // ch_versions   = ch_versions.mix(MEDAKA_INFERENCE.out.versions)
    ch_medaka_hdf = MEDAKA_INFERENCE.out.hdf

    //
    // CHANNEL: join up channels for synced input
    //
    if(multi_ref) {
        ch_hdf_ref = ch_medaka_hdf
            .map { [it[0].id, it ]}
            .join ( reference.map { [it[0].id, it[1], it[2]] })
            .map{ [it[1][0], it[1][1], it[2], it[3]] }
    } else {
        ch_hdf_ref = ch_medaka_hdf
            .combine(reference.map{[it[1], it[2]]})
    }

    //
    // MODULE: Generate vcf
    //
    MEDAKA_VCF (
        ch_hdf_ref.map{[it[0], it[1]]},
        ch_hdf_ref.map{[it[0], it[2], it[3]]}
    )
    // ch_versions   = ch_versions.mix(MEDAKA_VCF.out.versions)
    ch_medaka_vcf = MEDAKA_VCF.out.vcf

    //
    // MODULE: Join pooled VCF with reads
    //
    ch_pool_trimmed_bam_bai_vcf = ch_pool_trimmed_bam_bai
    .map { [it[1].id + "_" + it[1].pool, it ]}
    .join ( ch_medaka_vcf.map { [it[0].id + "_" + it[0].pool, it[1]] })
    .map{ [it[1][0], it[1][1], it[1][2], it[1][3], it[2]] }

    if(multi_ref) {
        ch_pool_trimmed_bam_bai_vcf_ref = ch_pool_trimmed_bam_bai_vcf
            .map { [it[1].id, it ]}
            .join ( reference.map { [it[0].id, it[1], it[2]] })
            .map{ [it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[2], it[3]] }
    } else {
        ch_pool_trimmed_bam_bai_vcf_ref = ch_pool_trimmed_bam_bai_vcf
            .combine(reference.map{[it[1], it[2]]})
    }

    //
    // MODULE: Generate vcf
    //
    MEDAKA_ANNOTATE (
        ch_pool_trimmed_bam_bai_vcf_ref.map{[it[1], it[4]]},
        ch_pool_trimmed_bam_bai_vcf_ref.map{[it[1], it[5], it[6]]},
        ch_pool_trimmed_bam_bai_vcf_ref.map{[it[1], it[2], it[3]]},
        ch_pool_trimmed_bam_bai_vcf_ref.map{[it[0]]}
    )
    // ch_versions   = ch_versions.mix(MEDAKA_ANNOTATE.out.versions)
    ch_medaka_vcf = MEDAKA_ANNOTATE.out.vcf

    if(pool_reads) {
        //
        // MODULE: Merge medaka vcfs
        //
        ch_vcf_pools = ch_medaka_vcf
            .map{[it[0].id, it[0], it[1], it[0].pool]}
            .groupTuple()
            .map{
                it[1][0].remove('pool')
                [it[1][0], it[2], it[3]]
            }
        ARTIC_VCF_MERGE (
            ch_vcf_pools,
            primer_bed
        )
        ch_medaka_vcf = ARTIC_VCF_MERGE.out.vcf
    }

    //
    // MODULE: Split VCF into pass/fail files passed on filter
    //
    BCFTOOLS_PASS_FAIL_SPLIT (
        ch_medaka_vcf
    )
    ch_versions   = ch_versions.mix(BCFTOOLS_PASS_FAIL_SPLIT.out.versions)
    ch_medaka_vcf = BCFTOOLS_PASS_FAIL_SPLIT.out.pass_vcf

    //
    // MODULE: Gzip and index the consensus VCF
    //
    INDEX_CONSEN_VCF (
        ch_medaka_vcf
    )
    ch_versions          = ch_versions.mix(INDEX_CONSEN_VCF.out.versions)
    ch_medaka_vcf_gz_tbi = INDEX_CONSEN_VCF.out.gz_tbi

    //
    // MODULE: Make depth mask
    //
    ARTIC_MAKE_DEPTH_MASK (
        primer_trimmed_bam_bai_ref.map{[it[0], it[1], it[2]]},
        primer_trimmed_bam_bai_ref.map{[it[0], it[3], it[4]]}
    )
    ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.mask

    //
    // MODULE: Build a pre-conensus mask of the coverage mask and N's where the failed variants are
    //
    ch_artic_mask = BCFTOOLS_PASS_FAIL_SPLIT.out.fail_vcf
        .map { [it[0].id, it ]}
        .join ( ch_depth_mask.map { [it[0].id, it[1]] })
        .map{ [it[1][0], it[1][1], it[2]] }

    if(multi_ref) {
        ch_artic_mask_ref = ch_artic_mask
            .map { [it[0].id, it ]}
            .join ( reference.map { [it[0].id, it[1], it[2]] })
            .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }
    } else {
        ch_artic_mask_ref = ch_artic_mask
            .combine(reference.map{[it[1], it[2]]})
    }
    ARTIC_MASK (
        ch_artic_mask_ref.map{[it[0], it[1]]},
        ch_artic_mask_ref.map{[it[0], it[2]]},
        ch_artic_mask_ref.map{[it[0], it[3], it[4]]},
    )
    ch_preconsensus_mask = ARTIC_MASK.out.fasta

    //
    // MODULE: Call the consensus sequence
    //
    ch_vcf_tbi_fasta_mask = ch_medaka_vcf_gz_tbi
        .map { [it[0].id, it ]}
        .join ( ch_preconsensus_mask.map { [it[0].id, it[1]] })
        .join ( ch_depth_mask.map { [it[0].id, it[1]] })
        .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }
    BCFTOOLS_CONSENSUS (
        ch_vcf_tbi_fasta_mask
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
    ch_consensus = BCFTOOLS_CONSENSUS.out.fasta

    // TODO Adjust header with proper fasta header

    //
    // MODULE: Run clair3 variant caller for more accurate variant calling
    //
    CLAIR3_RUN (
        primer_trimmed_bam_bai_ref.map{[it[0], it[1], it[2]]},
        primer_trimmed_bam_bai_ref.map{[it[0], it[3], it[4]]},
        clair3_model,
        clair3_platform
    )
    ch_versions          = ch_versions.mix(CLAIR3_RUN.out.versions)
    ch_clair3_vcf_gz_tbi = CLAIR3_RUN.out.merge_output_gz_tbi
    ch_clair3_vcf_unzip  = CLAIR3_RUN.out.merge_output_gz_tbi

    //
    // MODULE: filter ref calls from clair3
    //
    CLAIR3_FILTER_REF (
        ch_clair3_vcf_gz_tbi,
        [],
        [],
        []
    )
    ch_clair3_vcf = CLAIR3_FILTER_REF.out.vcf

    //
    // MODULE: Unzip clair3 VCF files
    //
    // GUNZIP_CLAIR3_VCF (
    //     ch_clair3_vcf_unzip.map{[it[0], it[1]]}
    // )
    // ch_versions   = ch_versions.mix(GUNZIP_CLAIR3_VCF.out.versions)
    // ch_clair3_vcf = GUNZIP_CLAIR3_VCF.out.gunzip

    //
    // MODULE: Call low frequency variants
    //
    LOFREQ_CALL (
        primer_trimmed_bam_bai_ref.map{[it[0], it[1], it[2]]},
        primer_trimmed_bam_bai_ref.map{[it[0], it[3], it[4]]},
    )
    ch_versions   = ch_versions.mix(LOFREQ_CALL.out.versions)
    ch_lofreq_vcf = LOFREQ_CALL.out.vcf

    //
    // MODULE: Gzip and index the lofreq VCF
    //
    INDEX_LOFREQ_VCF (
        ch_lofreq_vcf
    )
    ch_versions       = ch_versions.mix(INDEX_LOFREQ_VCF.out.versions)
    ch_lofreq_vcf_tbi = INDEX_LOFREQ_VCF.out.gz_tbi

    //
    // MODULE: Call structural variants
    //
    SNIFFLES (
        primer_trimmed_bam_bai_ref.map{[it[0], it[1], it[2]]},
        primer_trimmed_bam_bai_ref.map{[it[0], it[3]]},
        [[],[]],
        true,
        false
    )
    ch_versions         = ch_versions.mix(SNIFFLES.out.versions)
    ch_sniffles_vcf_gz  = SNIFFLES.out.vcf
    ch_sniffles_tbi     = SNIFFLES.out.tbi
    ch_sniffles_vcf_tbi = ch_sniffles_vcf_gz
        .map { [it[0].id, it ]}
        .join ( ch_sniffles_tbi.map { [it[0].id, it[1]] })
        .map{ [it[1][0], it[1][1], it[2]] }

    //
    // MODULE: Unzip sniffles VCF files
    //
    GUNZIP_SNIFFLES_VCF (
        ch_sniffles_vcf_gz
    )
    ch_versions = ch_versions.mix(GUNZIP_SNIFFLES_VCF.out.versions)

    //
    // MODULE: Rename consensus fasta
    //
    RENAME_FASTA (
        ch_consensus,
        [],
        true,
        "consensus"
    )
    ch_consensus = RENAME_FASTA.out.file

    //
    // CHANNEL: Generate merged vcf report channels
    //
    ch_vcf_files = ch_medaka_vcf.map{[it[0], it[1], "medaka", 1]}
        .mix(ch_clair3_vcf.map{[it[0], it[1], "clair3", 2]})
        .mix(ch_lofreq_vcf.map{[it[0], it[1], "lofreq", 3]})

    emit:
    versions         = ch_versions.ifEmpty(null)
    consensus        = ch_consensus
    medaka_vcf_tbi   = ch_medaka_vcf_gz_tbi
    clair3_vcf_tbi   = ch_clair3_vcf_gz_tbi
    lofreq_vcf_tbi   = ch_lofreq_vcf_tbi
    sniffles_vcf_tbi = ch_sniffles_vcf_tbi
    vcf_files        = ch_vcf_files
}
