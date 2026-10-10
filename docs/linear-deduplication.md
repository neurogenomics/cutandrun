# Linear duplicate removal: count and input contract

Linear duplicate removal chooses complete primary templates using the R1 alignment boundary. It does not establish physical original molecules, IVT families or biological cells. The helper uses no UMIs.

For ordinary bulk input without `CB`, the duplicate group remains `(reference, R1 strand, R1 boundary)`. A forward R1 uses its zero-based alignment start; a reverse R1 uses its exclusive reference end. This preserves the former BED boundary rule, including terminal deletions/skips. It is not a Tn5-shifted insertion coordinate or the first/last query-aligned base. Soft clips do not extend the boundary. Numeric R1 MAPQ chooses the representative; ties keep the first R1 in the coordinate-sorted input, regardless of when its mate arrives. MAPQ 255 remains the literal numeric value for compatibility and must not be interpreted as calibrated confidence.

When `CB` is present, the exact cell barcode also namespaces duplicate groups. This conditional change prevents collapsing records across different cells; it is not a cell call. `RG` distinguishes template identities but does not namespace duplicate groups, because sequencing lanes of the same bulk library may contain duplicates. Distinct physical libraries must not be pooled into one deduplication invocation without an explicit library-level rule.

Every primary template has one full `(QNAME, RG-or-absent, CB-or-absent)` identity and exactly one R1 and R2. Both records must be marked paired, mapped, and have mapped mates; they must agree on reciprocal mate coordinates, strand flags, equal-and-opposite TLEN and proper-pair status. Nonproper complete pairs remain supported; zero TLEN is allowed. Mapped coordinates and reference-consuming CIGAR spans must be in header bounds, there must be a query-aligned base, and present SEQ/QUAL lengths must agree with the CIGAR. The rule does not infer mates, strip names, check aligner-specific proper-pair geometry, or compare TLEN magnitude with alignment span.

Repeated primary roles, identity reuse, missing mates, conflicting identity tags or inconsistent mate metadata fail the task. In particular, readwise upstream filtering can leave an orphan: the deduplicator cannot reconstruct its lost mate. The separately documented opt-in `--la_complete_templates` candidate filter retains an original pair only when both exact primary records survived the existing readwise filters. That filter changes the retained population and must be distinguished from historical processing.

Secondary and supplementary records cannot select a winner. Such records are retained only if their exact identity belongs to a selected primary pair. Records without such a primary decision are dropped and counted; no ancestry is guessed from a shared name with different or missing tags. Associated records do not contribute to primary-template counts.

## Implementation and metrics

The subworkflow converts coordinate-sorted BAM to SAM using the existing pinned SAMtools module. `find_unique_reads.py` validates pairs and selects winners, then writes unchanged SAM records for selected exact identities. SAMtools sorts and indexes the result. The old BED/name-list interface has been removed because it cannot retain read-group/cell identity or distinguish primary from secondary records.

The legacy metric labels remain compatible, but `Reads before filtering`, `LA duplicates removed (n)` and `Unique reads after LA duplicate removal` now explicitly count complete primary templates, not all SAM records. Additional fields state the count unit, identity rule, source SAM SHA-256, primary/nonprimary counts, mate-cache peak, spill events and SQLite bytes. Empty input yields zero templates and `NA` duplicate percentage.

The in-memory pending-mate cache defaults to 100,000 identities. Older pending mates and every completed identity are stored in SQLite; the page cache is 8 MiB. This avoids an unbounded Python dictionary, but temporary disk grows with the number of templates. Sequential SAM intermediates are uncompressed and may substantially increase disk and I/O. The local synthetic tests do not establish full-depth memory, speed or scratch-disk requirements; measure all three before a production rollout. Set helper `--max_pending` and `--tmpdir` through `FIND_UNIQUE_READS` `ext.args` when needed.

Input SAM bytes and file identity must agree across validation and emission. Output files are staged until both scans succeed, so validation failures preserve existing files. Final file renames are sequential, not a filesystem transaction; a rename failure can publish only part of the output set, and Nextflow must treat the failed task as failed.

## Developer validation

```bash
python -m unittest discover -s tests/unit -p test_find_unique_reads.py -v
nextflow -C tests/integration/la_dedup/nextflow.config \
  run tests/integration/la_dedup/main.nf \
  -work-dir /new/private/tiny-work --outdir /new/private/tiny-results
```

The unit suite covers numeric MAPQ, same-name identities across read groups/cells, secondary-record resurrection, R1-order ties, reverse boundary equivalence, invalid/incomplete primary pairs, forced disk spill and changed-input rejection. Where SAMtools and BEDTools are installed, it also compares full retained SAM records with frozen legacy revision `8e3bc3a753767282247b1ad5c5de2fe92f5231ce` on valid ordinary bulk fixtures. The tiny Nextflow harness checks the actual subworkflow against four expected records representing two complete templates; enable its optional pair-filter mode to exercise asymmetric readwise filtering too. It uses synthetic data only.
