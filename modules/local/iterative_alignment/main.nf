process ITERATIVE_ALIGNMENT {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/thecrick/pipetech_iterative_alignment:0.2.22"

    input:
    tuple val(meta), path(reads), path(ref)

    output:
    tuple val(meta), path("${meta.id}/*.bam")              , emit: bam
    tuple val(meta), path("${meta.id}/*.bai")              , emit: bai
    tuple val(meta), path("${meta.id}/*consensus.fasta")   , emit: consensus_wref
    tuple val(meta), path("${meta.id}/*consensus.wn.fasta"), emit: consensus_wn
    tuple val(meta), path("${meta.id}/*final_ref.fasta")   , emit: final_ref
    tuple val(meta), path("${meta.id}/*.csv")  , emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    max_iterations = task.ext.max_iterations ?: 1
    min_cov        = task.ext.min_cov ?: 10
    var_thresh     = task.ext.var_thresh ?: 0.001
    """
    #!/opt/conda/bin/python

    import logging
    from crick_genome_tools.workflows.iterative_alignment import Aligner, IterationMode, IterativeAlignment

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

    iter_align = IterativeAlignment(
        num_cores=${task.cpus},
        output_path=".",
        min_iterations=1,
        max_iterations=${max_iterations},
        aligner=Aligner.BWA,
        iteration_mode=IterationMode.COUNT,
        var_thresh=${var_thresh},
        min_coverage=${min_cov},
        bwa_args=["-T10"],
    )

    iter_align.run_sample(
        "${meta.id}",
        "${reads[0]}",
        "${reads[1]}",
        "${ref}",
    )
    """
}
