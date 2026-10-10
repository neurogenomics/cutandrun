process LA_PAIR_INPUTS {
    tag "$meta.id"
    label 'process_low'
    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(original), path(accepted)

    output:
    tuple val(meta), path('original.sam'), path('accepted.sam'), emit: sam
    path 'versions.yml', emit: versions

    script:
    """
    samtools view -h --no-PG --threads ${task.cpus-1} -o original.sam "$original"
    samtools view -h --no-PG --threads ${task.cpus-1} -o accepted.sam "$accepted"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
