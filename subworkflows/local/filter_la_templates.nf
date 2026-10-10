/* Explicit candidate method: keep both original primary mates only if both
 * survive the existing MAPQ/region/mitochondrial filters. */
include { LA_PAIR_INPUTS } from '../../modules/local/samtools/la_pair_inputs/main'
include { FILTER_LA_TEMPLATES } from '../../modules/local/python/filter_la_templates'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow COMPLETE_LA_TEMPLATES {
    take:
    original
    accepted
    fasta

    main:
    ch_pairs = original.join(accepted, by: [0], failOnDuplicate: true, failOnMismatch: true)
    LA_PAIR_INPUTS(ch_pairs)
    FILTER_LA_TEMPLATES(LA_PAIR_INPUTS.out.sam)
    BAM_SORT_STATS_SAMTOOLS(FILTER_LA_TEMPLATES.out.sam, fasta)

    emit:
    bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai = BAM_SORT_STATS_SAMTOOLS.out.bai
    stats = BAM_SORT_STATS_SAMTOOLS.out.stats
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    metrics = FILTER_LA_TEMPLATES.out.metrics
    versions = LA_PAIR_INPUTS.out.versions.mix(FILTER_LA_TEMPLATES.out.versions).mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
}
