nextflow.enable.dsl=2

include { SAMTOOLS_VIEW_SORT_STATS as FILTER_READS } from '../../../subworkflows/local/samtools_view_sort_stats'
include { COMPLETE_LA_TEMPLATES } from '../../../subworkflows/local/filter_la_templates'
include { DEDUPLICATE_LINEAR } from '../../../subworkflows/local/deduplicate_linear'
include { SAMTOOLS_VIEW as SAMTOOLS_MAKE_BAM } from '../../../modules/nf-core/samtools/view/main'

process CHECK_LINEAR_OUTPUT {
    container 'biocontainers/samtools:1.17--h00cdaf9_0'
    input:
    tuple val(meta), path(bam)
    path expected
    output:
    path 'verified.txt'
    script:
    """
    samtools quickcheck $bam
    samtools view $bam > actual.sam
    diff -u $expected actual.sam
    test \$(samtools view -c -f 64 $bam) -eq 2
    test \$(samtools view -c -f 128 $bam) -eq 2
    echo 'PASS: exact expected four SAM records; two complete primary templates' > verified.txt
    """
}

workflow CUTANDRUN {
    main:
    meta = [id:'fixture', single_end:false, is_control:false]
    ch_sam = Channel.of([meta, file(params.test_pair_filter ? "$projectDir/fixtures/input-with-filter-failure.sam" : "$projectDir/fixtures/input.sam"), []])
    ch_fasta = Channel.value([[id:'reference'], file("$projectDir/fixtures/reference.fa")])
    SAMTOOLS_MAKE_BAM(ch_sam, ch_fasta, [])
    ch_input = SAMTOOLS_MAKE_BAM.out.bam
    if (params.test_pair_filter) {
        FILTER_READS(ch_input, [], ch_fasta)
        COMPLETE_LA_TEMPLATES(ch_input, FILTER_READS.out.bam, ch_fasta)
        ch_input = COMPLETE_LA_TEMPLATES.out.bam
    }
    DEDUPLICATE_LINEAR(ch_input, Channel.empty(), ch_fasta, Channel.empty(), true,
                      file("$projectDir/../../../assets/multiqc/linear_duplication_header.txt"))
    CHECK_LINEAR_OUTPUT(DEDUPLICATE_LINEAR.out.bam, file("$projectDir/fixtures/expected.sam"))
}

workflow NFCORE_CUTANDRUN {
    main:
    CUTANDRUN()
}

workflow {
    NFCORE_CUTANDRUN()
}
