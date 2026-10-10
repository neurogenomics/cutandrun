import importlib.util
import io
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('linear', ROOT / 'bin/find_unique_reads.py')
linear = importlib.util.module_from_spec(spec)
spec.loader.exec_module(linear)
HEADER = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n@RG\tID:rg1\tSM:bulk\n@RG\tID:rg2\tSM:bulk\n'


def pair(name, start=101, mapq=42, rg='rg1', cb=None, reverse=False, secondary=False):
    tags = ([] if rg is None else ['RG:Z:' + rg]) + ([] if cb is None else ['CB:Z:' + cb])
    rows = []
    for first in (True, False):
        flag = (83 if first else 163) if reverse else (99 if first else 147)
        if secondary: flag |= 256
        values = [name, flag, 'chr1', start if first else start+40, mapq, '20M', '=',
                  start+40 if first else start, 60 if first else -60, 'A'*20, 'I'*20, *tags]
        rows.append('\t'.join(str(v) for v in values) + '\n')
    return rows


def field(row, index, value):
    values = row.rstrip('\n').split('\t')
    values[index] = str(value)
    return '\t'.join(values) + '\n'


class LinearTemplateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source, self.output, self.metrics, self.header, self.mqc = [self.root / name for name in
            ('input.sam', 'output.sam', 'metrics.txt', 'header.txt', 'mqc.tsv')]
        self.header.write_text('# fixture\n')

    def execute(self, rows, **kwargs):
        self.source.write_text(HEADER + ''.join(rows))
        counts = linear.select_templates(self.source, self.output, self.metrics, self.header, self.mqc, **kwargs)
        records = [line.split('\t') for line in self.output.read_text().splitlines() if not line.startswith('@')]
        return counts, records

    def reject(self, rows, pattern=None):
        for path in (self.output, self.metrics, self.mqc): path.write_text('previous valid output\n')
        with self.assertRaisesRegex(ValueError, pattern or '.'):
            self.execute(rows)
        for path in (self.output, self.metrics, self.mqc): self.assertEqual(path.read_text(), 'previous valid output\n')

    def test_numeric_mapq_and_original_records(self):
        for ordered in (pair('low', mapq=9)+pair('high'), pair('high')+pair('low', mapq=9)):
            counts, records = self.execute(ordered)
            self.assertEqual([r[0] for r in records], ['high', 'high'])
            self.assertEqual(counts['retained_templates'], 1)
            self.assertEqual(self.output.read_text(), HEADER + ''.join(pair('high')))

    def test_same_qname_readgroup_does_not_resurrect_loser(self):
        counts, records = self.execute(pair('A') + pair('B', mapq=30) + pair('B', start=301, rg='rg2'))
        self.assertEqual(counts['retained_templates'], 2)
        self.assertEqual([(r[0], r[3], r[11]) for r in records],
                         [('A','101','RG:Z:rg1'),('A','141','RG:Z:rg1'),('B','301','RG:Z:rg2'),('B','341','RG:Z:rg2')])

    def test_nonprimary_cannot_choose_or_resurrect_templates(self):
        rows = pair('A') + pair('B', mapq=30) + pair('B', start=301, secondary=True)
        rows += pair('A', start=501, secondary=True) + pair('unknown', start=701, secondary=True)
        counts, records = self.execute(rows)
        self.assertEqual([r[0] for r in records], ['A']*4)
        self.assertEqual(counts['retained_templates'], 1)
        self.assertEqual(counts['retained_nonprimary_records'], 2)
        self.assertEqual(counts['dropped_nonprimary_records'], 4)

    def test_cb_partitions_anchors_but_rg_does_not(self):
        counts, _ = self.execute(pair('same', cb='CELL1') + pair('same', cb='CELL2'))
        self.assertEqual(counts['retained_templates'], 2)
        counts, records = self.execute(pair('same', rg='rg1', mapq=20) + pair('same', rg='rg2', mapq=42))
        self.assertEqual(counts['retained_templates'], 1)
        self.assertTrue(all('RG:Z:rg2' in r for r in records))

    def test_tie_is_first_r1_not_first_complete_pair(self):
        a, b = pair('A'), pair('B')
        _, records = self.execute([a[0], b[0], b[1], a[1]])
        self.assertEqual([r[0] for r in records], ['A', 'A'])

    def test_reverse_anchor_and_opposite_strands_remain_distinct(self):
        a, b = pair('A', start=101, reverse=True), pair('B', start=91, reverse=True)
        b[0] = field(b[0], 5, '30M')
        b[0] = field(field(b[0], 9, 'A'*30), 10, 'I'*30)
        counts, records = self.execute(a+b+pair('forward', start=101))
        self.assertEqual(counts['retained_templates'], 2)
        self.assertEqual([r[0] for r in records], ['A','A','forward','forward'])

    def test_orphans_duplicate_roles_and_late_identity_reuse_fail(self):
        a = pair('A')
        for rows in ([a[0]], [a[1]], [a[0], a[0]], a+pair('B', start=301)+a):
            with self.subTest(rows=rows): self.reject(rows)

    def test_mate_identity_and_reciprocity_fail_closed(self):
        for kind in ('rg', 'cb', 'coordinate', 'strand', 'tlen', 'proper'):
            a = pair('A')
            if kind == 'rg': a[1] = a[1].replace('RG:Z:rg1', 'RG:Z:rg2')
            if kind == 'cb': a[1] = a[1].rstrip('\n') + '\tCB:Z:CELL\n'
            if kind == 'coordinate': a[1] = field(a[1], 7, 999)
            if kind == 'strand': a[1] = field(a[1], 1, 147|32)
            if kind == 'tlen': a[1] = field(a[1], 8, -61)
            if kind == 'proper': a[1] = field(a[1], 1, 147 & ~2)
            with self.subTest(kind=kind): self.reject(a)

    def test_nonproper_zero_tlen_and_absent_readgroup_are_supported(self):
        a = pair('A', rg=None)
        a = [field(field(row, 1, int(row.split('\t')[1]) & ~2), 8, 0) for row in a]
        counts, records = self.execute(a)
        self.assertEqual(counts['retained_templates'], 1)
        self.assertEqual(len(records), 2)

    def test_invalid_mapped_records_or_tags_fail(self):
        for kind in ('cigar', 'span', 'deletion_only', 'skip_only', 'query_length', 'quality_length',
                     'position', 'flag', 'mapq', 'duplicate_tag', 'tag_type'):
            a = pair('A')
            if kind == 'cigar': a[0] = field(a[0], 5, '*')
            if kind == 'span': a[0] = field(a[0], 5, '20S')
            if kind == 'deletion_only': a[0] = field(a[0], 5, '20D20S')
            if kind == 'skip_only': a[0] = field(a[0], 5, '20N20S')
            if kind == 'query_length': a[0] = field(a[0], 5, '10M')
            if kind == 'quality_length': a[0] = field(a[0], 10, 'I'*19)
            if kind == 'position': a[0] = field(a[0], 3, 10001)
            if kind == 'flag': a[0] = field(a[0], 1, -1)
            if kind == 'mapq': a[0] = field(a[0], 4, 300)
            if kind == 'duplicate_tag': a[0] = a[0].rstrip('\n')+'\tRG:Z:rg1\n'
            if kind == 'tag_type': a[0] = a[0].replace('RG:Z:rg1', 'RG:i:1')
            with self.subTest(kind=kind): self.reject(a)

    def test_empty_and_nonprimary_only_do_not_divide_by_zero(self):
        for rows in ([], pair('A', secondary=True)):
            counts, records = self.execute(rows)
            self.assertEqual(records, [])
            self.assertEqual(counts['retained_templates'], 0)
            self.assertIn('LA duplicates removed (%)\tNA', self.metrics.read_text())

    def test_pending_spill_is_lossless(self):
        pairs = [pair('read%d' % i, start=101+i*100) for i in range(5)]
        rows = [p[0] for p in pairs] + [p[1] for p in pairs]
        counts, _ = self.execute(rows, max_pending=1)
        self.assertEqual(counts['complete_templates'], 5)
        self.assertEqual(counts['pending_spill_events'], 4)
        self.assertEqual(counts['peak_pending_in_memory'], 1)
        self.assertEqual(self.output.read_text(), HEADER+''.join(rows))

    def test_second_scan_changes_preserve_old_outputs(self):
        rows = pair('A')
        self.source.write_text(HEADER+''.join(rows))
        for path in (self.output, self.metrics, self.mqc): path.write_text('previous\n')
        original_open = Path.open
        reads = 0
        def opened(path, *args, **kwargs):
            nonlocal reads
            if path == self.source and args == ('rb',):
                reads += 1
                if reads == 2: return io.BytesIO((HEADER+''.join(pair('B'))).encode())
            return original_open(path, *args, **kwargs)
        with patch.object(Path, 'open', opened), self.assertRaisesRegex(ValueError, 'changed'):
            linear.select_templates(self.source, self.output, self.metrics, self.header, self.mqc)
        for path in (self.output, self.metrics, self.mqc): self.assertEqual(path.read_text(), 'previous\n')
        self.assertFalse(list(self.root.glob('.output.sam*')))

    @unittest.skipUnless(shutil.which('bedtools') and shutil.which('samtools'), 'bedtools/samtools integration tools unavailable')
    def test_valid_bulk_output_matches_frozen_legacy_full_records(self):
        rows = pair('low', mapq=20)+pair('high')+pair('distinct', start=501)+pair('reverse',start=701,reverse=True)
        self.execute(rows)
        old = self.root/'legacy.py'
        old.write_bytes(subprocess.check_output(['git','-C',str(ROOT),'show','8e3bc3a753767282247b1ad5c5de2fe92f5231ce:bin/find_unique_reads.py']))
        bam, bed, names = (self.root/name for name in ('input.bam','input.bed','names.txt'))
        subprocess.run(['samtools','view','-b','-o',str(bam),str(self.source)],check=True,capture_output=True)
        bed.write_bytes(subprocess.check_output(['bedtools','bamtobed','-i',str(bam)]))
        subprocess.run([sys.executable,str(old),'--bed_path',str(bed),'--output_path',str(names),
                        '--metrics_path',str(self.root/'old.metrics'),'--header_path',str(self.header),
                        '--mqc_path',str(self.root/'old.mqc')],check=True,capture_output=True)
        legacy = subprocess.check_output(['samtools','view','-N',str(names),str(bam)]).decode()
        candidate = ''.join(line for line in self.output.read_text().splitlines(True) if not line.startswith('@'))
        self.assertEqual(candidate, legacy)


if __name__ == '__main__': unittest.main()
