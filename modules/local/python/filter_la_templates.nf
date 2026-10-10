process FILTER_LA_TEMPLATES {
    tag "$meta.id"
    label 'process_single'
    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(original), path(accepted)

    output:
    tuple val(meta), path('*.complete_templates.sam'), emit: sam
    tuple val(meta), path('*.la_filter.json'), emit: metrics
    path 'versions.yml', emit: versions

    script:
    """
    filter_la_templates.py \\
        --original "$original" --accepted "$accepted" \\
        --output "${meta.id}.complete_templates.sam" \\
        --report "${meta.id}.la_filter.json" --tmpdir .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
