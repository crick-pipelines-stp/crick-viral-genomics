process EXPORT_REPORT_DATA {
    tag "$run_id"
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:0.5.13"

    input:
    val(run_id)
    val(report_type)
    val(summary)
    path("data/ref/*")
    path("data/annotation/*")
    path("data/samplesheet/*")
    path("data/toulligqc/*")
    path("data/samtools_host/*")
    path("data/samtools_contaminent/*")
    path("data/samtools_alignment/*")
    path("data/coverage/*")
    path("data/consensus/*")
    tuple val(vcf_tools), path("data/variants/*")
    path("data/variants_compressed/*")
    path("data/count_table/*")
    path("data/truncation/*")

    output:
    path("*.pkl"), emit: pkl
    path("*.sh") , emit: sh

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/local/bin/python

    import json
    import logging
    import os
    from crick_genome_tools.reporting.report_data_parser import ReportDataParser

    logging.basicConfig(level=logging.INFO)

    vcf_tools = ${vcf_tools.collect{"\"${it}\""}}

    summary_dict = json.loads('${summary}')
    with open("./data/summary.json", "w") as f:
        json.dump(summary_dict, f)

    parser = ReportDataParser("data")
    parser.get_data(vcf_tools)
    parser.save_data("${run_id}.pkl")

    run_script = '#!/bin/bash\\n\\nSCRIPT_DIR=\$(dirname "\$(realpath "\$0")")\\ndocker run -p 8501:8501 -v "\$SCRIPT_DIR:/data" crick/st_report_server:latest streamlit run app/app.py --server.port 8501 --server.headless true -- --data_path /data/${run_id}.pkl --report_type ${report_type}'

    rep_path = "run_interative_report.sh"
    with open(rep_path, "w") as f:
        f.write(run_script)
    os.chmod(rep_path, 0o755)
    """
}
