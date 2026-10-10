/*
 * Select complete primary templates by the R1 aligned 5' anchor. Retain lossless
 * SAM identity/flags/tags throughout: BED names plus a QNAME allowlist cannot
 * distinguish primary templates that share a name across read groups/cells.
 */
include { FIND_UNIQUE_READS } from '../../modules/local/python/find_unique_reads'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_LINEAR_SAM } from '../../modules/nf-core/samtools/view/main.nf'

workflow DEDUPLICATE_LINEAR {
    take:
    bam
    bai            // Retained interface; sequential SAM conversion needs no index.
    fasta
    fai            // Retained interface for callers.
    process_target
    mqc_header

    main:
    ch_versions = Channel.empty()
    ch_targets = Channel.empty()
    if (process_target) {
        ch_process = bam
    } else {
        bam.branch {
            target: it[0].is_control == false
            control: it[0].is_control == true
        }.set { ch_split }
        ch_process = ch_split.control
        ch_targets = ch_split.target
    }

    // Preserve the previous coordinate-sorted R1 tie order.
    SAMTOOLS_SORT(ch_process)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    SAMTOOLS_VIEW_LINEAR_SAM(
        SAMTOOLS_SORT.out.bam.map { meta, file -> [meta, file, []] },
        fasta,
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_LINEAR_SAM.out.versions)
    FIND_UNIQUE_READS(SAMTOOLS_VIEW_LINEAR_SAM.out.sam, mqc_header)
    ch_versions = ch_versions.mix(FIND_UNIQUE_READS.out.versions)
    ch_selected = FIND_UNIQUE_READS.out.sam
    if (!process_target) {
        // Preserve stable target/control order on resume as before.
        ch_targets = ch_targets.toSortedList { row -> row[0].id }.flatMap()
        ch_selected = ch_targets.concat(ch_selected.toSortedList { row -> row[0].id }.flatMap())
    }
    BAM_SORT_STATS_SAMTOOLS(ch_selected, fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai = BAM_SORT_STATS_SAMTOOLS.out.bai
    stats = BAM_SORT_STATS_SAMTOOLS.out.stats
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    metrics = FIND_UNIQUE_READS.out.metrics
    linear_metrics_mqc = FIND_UNIQUE_READS.out.linear_metrics_mqc
    versions = ch_versions
}
