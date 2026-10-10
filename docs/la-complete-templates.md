# Candidate complete-template filtering for linear amplification

`--la_complete_templates true` is an explicit candidate method change for paired-end LA workflows. It defaults to false and requires active read filtering and linear duplicate removal. It applies to both target and control samples. Workflows without LA duplicate removal reject this flag; their existing filtering path is unchanged.

Historical filtering applies MAPQ, allowed-region and mitochondrial predicates to each record. A pair with MAPQ 42 and 9 at a threshold of 10 therefore leaves one record. The candidate evaluates those same existing predicates, but retains both original primary records only when both survive. It discards secondary/supplementary records and incomplete original templates. It never recovers a missing mate or guesses pair identity from an alignment coordinate. Identity is the exact tuple of QNAME, RG (or absence), and CB (or absence).

The original pre-filter BAM and the existing filtered BAM are converted to SAM. The report binds both SAM inputs to SHA256 values computed during their reads and before/after file identities; source changes during filtering fail the task. Reference `@SQ` dictionaries must match, while harmless sort/program header differences are allowed. The helper compares complete primary records, including flags and tags, and rejects foreign or altered accepted records and repeated primary identity/mate combinations. A bounded pending-mate cache spills to SQLite. The LA deduplicator subsequently validates reciprocal mate geometry; this filter alone does not certify structural integrity. Readwise MAPQ or blacklist predicates should be described as alignment-record criteria; requiring both records to survive does not establish that the entire physical fragment avoids a blacklist.

`02_alignment/<aligner>/la_template_filter/<sample>.la_filter.json` reports original/accepted record counts, complete and incomplete original templates, retained templates, templates rejected because one or both mates failed, nonprimary/unpaired exclusions, SQLite bytes and pending-cache spills. The original readwise filter statistics remain separately named. Downstream alignment statistics describe the complete-template BAM. Biological comparisons must include these filter fates and distinguish this option from the unchanged historical readwise baseline.

This implementation prioritises an auditable boundary over speed: it writes two uncompressed SAMs and a temporary SQLite database containing accepted records, then sorts the selected SAM. Disk demand and full-depth runtime have not been benchmarked. Reserve sufficient writable task scratch and retain trace/resource metrics before adopting it for large libraries. No production launchers enable this option automatically.

Run the focused regression suite with:

```bash
python3 -m unittest discover -s tests/unit -p 'test_filter_la_templates.py' -v
```

Tests cover asymmetric MAPQ filtering through the real SAMtools CLI when available, RG/CB collisions, whole-template rejection, nonprimary records, missing mates, disk spill, duplicate/foreign records, and refusing to overwrite existing products.
