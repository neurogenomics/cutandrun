process FIND_UNIQUE_READS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(input)
    path mqc_header

    output:
    tuple val(meta), path('*.linear_dedup.sam'), emit: sam
    tuple val(meta), path('*metrics.txt'), emit: metrics
    tuple val(meta), path('*mqc.tsv'), emit: linear_metrics_mqc
    path  "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find_unique_reads.py \\
        --sam_path $input \\
        --output_path "${prefix}.linear_dedup.sam" \\
        --metrics_path "${prefix}_metrics.txt" \\
        --header_path $mqc_header \\
        --mqc_path "${prefix}_mqc.tsv" \\
        --tmpdir . \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """
}
