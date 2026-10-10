#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2023 @chris-cheshire
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


"""Retain complete primary templates by R1 alignment-boundary anchor (Python 3.8+).

Optional CB namespaces anchors. RG disambiguates identity, not duplicate groups.
The coordinate rule is an assay heuristic, not proof of original molecules.
"""
import argparse
from collections import Counter, OrderedDict
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile


def file_identity(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def template_identity(fields):
    if fields[0] == '*' or not re.fullmatch(r'[!-?A-~]{1,254}', fields[0]):
        raise ValueError('Missing QNAME')
    tags = {}
    for item in fields[11:]:
        if item[:3] in ('RG:', 'CB:'):
            parts = item.split(':', 2)
            if len(parts) != 3 or parts[1] != 'Z' or not parts[2] or parts[0] in tags:
                raise ValueError('RG/CB must be unique nonempty Z tags')
            tags[parts[0]] = parts[2]
    return json.dumps([fields[0], tags.get('RG'), tags.get('CB')], separators=(',', ':'))


def parse_primary(fields, references, ordinal):
    flag, pos, mapq, mate_pos, tlen = (int(fields[i]) for i in (1, 3, 4, 7, 8))
    if not flag & 1 or flag & (4 | 8) or bool(flag & 64) == bool(flag & 128):
        raise ValueError('Primary must be mapped paired R1 or R2 with mapped mate')
    if not 0 <= mapq <= 255:
        raise ValueError('MAPQ outside SAM range')
    chrom, cigar = fields[2], fields[5]
    mate_chrom = chrom if fields[6] == '=' else fields[6]
    if (chrom not in references or mate_chrom not in references
            or not 1 <= pos <= references[chrom] or not 1 <= mate_pos <= references[mate_chrom]):
        raise ValueError('Mapped coordinate outside reference bounds')
    if not re.fullmatch(r'(?:[1-9][0-9]*[MIDNSHP=X])+', cigar):
        raise ValueError('Invalid mapped CIGAR')
    operations = [(int(n), op) for n, op in re.findall(r'([0-9]+)([MIDNSHP=X])', cigar)]
    span = sum(n for n, op in operations if op in 'MDN=X')
    if not any(op in 'M=X' for _, op in operations) or pos - 1 + span > references[chrom]:
        raise ValueError('Invalid mapped reference span')
    query_length = sum(n for n, op in operations if op in 'MIS=X')
    seq, qual = fields[9:11]
    if ((seq == '*' and qual != '*') or (seq != '*' and len(seq) != query_length)
            or (qual != '*' and len(qual) != len(seq))):
        raise ValueError('CIGAR, SEQ and QUAL lengths disagree')
    return dict(flag=flag, role=1 if flag & 64 else 2, chrom=chrom, pos=pos,
                end=pos - 1 + span, mate_chrom=mate_chrom, mate_pos=mate_pos,
                mapq=mapq, tlen=tlen, ordinal=ordinal)


def check_pair(a, b):
    if a['role'] == b['role']:
        raise ValueError('Duplicate primary mate role')
    for left, right in ((a, b), (b, a)):
        if left['mate_chrom'] != right['chrom'] or left['mate_pos'] != right['pos']:
            raise ValueError('Nonreciprocal primary mate coordinates')
        if bool(left['flag'] & 32) != bool(right['flag'] & 16):
            raise ValueError('Nonreciprocal primary mate orientation')
    if a['tlen'] != -b['tlen'] or bool(a['flag'] & 2) != bool(b['flag'] & 2):
        raise ValueError('Primary mate TLEN/proper-pair flags disagree')
    return a if a['role'] == 1 else b


def select_templates(source, output, metrics, header, mqc, max_pending=100000, tmpdir=None):
    """Fail closed on invalid pairs; bounded mate cache plus disk identity index.

    Secondary/supplementary records cannot choose winners. They survive only
    when their exact QNAME/RG/CB identity belongs to a selected primary pair.
    All output files are staged until validation and emission agree on input
    bytes. Their final renames are sequential, not a multi-file transaction.
    """
    if max_pending < 1:
        raise ValueError('max_pending must be positive')
    source, output, metrics, header, mqc = [Path(p) for p in (source, output, metrics, header, mqc)]
    if len({p.resolve() for p in (source, output, metrics, header, mqc)}) != 5:
        raise ValueError('Input/output paths must be distinct')
    before, first_hash = file_identity(source), hashlib.sha256()
    counts, references, pending = Counter(), {}, OrderedDict()
    with tempfile.TemporaryDirectory(prefix='linear-dedup-', dir=tmpdir) as scratch:
        database = Path(scratch) / 'templates.sqlite'
        db = sqlite3.connect(str(database))
        try:
            db.executescript("""
                PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
                PRAGMA cache_size=-8192; PRAGMA mmap_size=0;
                CREATE TABLE completed (identity TEXT PRIMARY KEY) WITHOUT ROWID;
                CREATE TABLE pending (identity TEXT PRIMARY KEY, record TEXT) WITHOUT ROWID;
                CREATE TABLE winners (anchor TEXT PRIMARY KEY, identity TEXT, mapq INTEGER, ordinal INTEGER) WITHOUT ROWID;
            """)
            seen_record = False
            with source.open('rb') as handle:
                for raw in handle:
                    first_hash.update(raw)
                    line = raw.decode('utf-8').rstrip('\r\n')
                    if line.startswith('@'):
                        if seen_record:
                            raise ValueError('SAM header encountered after alignments')
                        if line.startswith('@SQ\t'):
                            tags = dict(item.split(':', 1) for item in line.split('\t')[1:])
                            name, length = tags['SN'], int(tags['LN'])
                            if name in references or length <= 0:
                                raise ValueError('Invalid or duplicate @SQ reference')
                            references[name] = length
                        continue
                    seen_record = True
                    fields = line.split('\t')
                    if len(fields) < 11:
                        raise ValueError('SAM alignment has fewer than 11 fields')
                    key = template_identity(fields)
                    counts['input_alignment_records'] += 1
                    flag = int(fields[1])
                    if not 0 <= flag <= 65535:
                        raise ValueError('FLAG outside SAM range')
                    if flag & (256 | 2048):
                        counts['input_nonprimary_records'] += 1
                        continue
                    record = parse_primary(fields, references, counts['input_alignment_records'])
                    counts['primary_records'] += 1
                    counts['primary_r%d' % record['role']] += 1
                    if db.execute('SELECT 1 FROM completed WHERE identity=?', (key,)).fetchone():
                        raise ValueError('Repeated completed primary template identity')
                    mate = pending.pop(key, None)
                    if mate is None:
                        spilled = db.execute('SELECT record FROM pending WHERE identity=?', (key,)).fetchone()
                        if spilled:
                            mate = json.loads(spilled[0])
                            db.execute('DELETE FROM pending WHERE identity=?', (key,))
                    if mate is None:
                        pending[key] = record
                        if len(pending) > max_pending:
                            old_key, old_record = pending.popitem(last=False)
                            db.execute('INSERT INTO pending VALUES (?,?)', (old_key, json.dumps(old_record)))
                            counts['pending_spill_events'] += 1
                        counts['peak_pending_in_memory'] = max(counts['peak_pending_in_memory'], len(pending))
                        continue
                    r1 = check_pair(mate, record)
                    db.execute('INSERT INTO completed VALUES (?)', (key,))
                    reverse = bool(r1['flag'] & 16)
                    anchor = json.dumps([json.loads(key)[2], r1['chrom'], reverse, r1['end'] if reverse else r1['pos'] - 1])
                    winner = db.execute('SELECT mapq,ordinal FROM winners WHERE anchor=?', (anchor,)).fetchone()
                    if winner is None:
                        db.execute('INSERT INTO winners VALUES (?,?,?,?)', (anchor, key, r1['mapq'], r1['ordinal']))
                    elif (r1['mapq'], -r1['ordinal']) > (winner[0], -winner[1]):
                        db.execute('UPDATE winners SET identity=?,mapq=?,ordinal=? WHERE anchor=?', (key, r1['mapq'], r1['ordinal'], anchor))
                    counts['complete_templates'] += 1
                    if counts['complete_templates'] % 10000 == 0:
                        db.commit()
            if pending or db.execute('SELECT 1 FROM pending LIMIT 1').fetchone():
                raise ValueError('Orphan primary: both mates must survive filtering with identical QNAME/RG/CB')
            if file_identity(source) != before:
                raise ValueError('Input SAM changed during validation')
            db.execute('CREATE UNIQUE INDEX winner_identity ON winners(identity)')
            counts['retained_templates'] = db.execute('SELECT COUNT(*) FROM winners').fetchone()[0]
            db.commit()
            counts['sqlite_bytes'] = database.stat().st_size
            staged = []
            try:
                with tempfile.NamedTemporaryFile(dir=output.parent, prefix='.' + output.name, delete=False) as sam:
                    staged.append((Path(sam.name), output))
                    second_hash = hashlib.sha256()
                    with source.open('rb') as handle:
                        for raw in handle:
                            second_hash.update(raw)
                            if raw.startswith(b'@'):
                                sam.write(raw)
                                continue
                            fields = raw.decode('utf-8').rstrip('\r\n').split('\t')
                            key = template_identity(fields)
                            if db.execute('SELECT 1 FROM winners WHERE identity=?', (key,)).fetchone():
                                sam.write(raw)
                                counts['retained_nonprimary_records' if int(fields[1]) & (256 | 2048) else 'retained_primary_records'] += 1
                    if file_identity(source) != before or first_hash.digest() != second_hash.digest():
                        raise ValueError('Input SAM changed between validation and emission')
                if counts['retained_primary_records'] != 2 * counts['retained_templates']:
                    raise ValueError('Output primary-template counts do not reconcile')
                total, retained = counts['complete_templates'], counts['retained_templates']
                mqc_line = 'LA duplicates removed (%)\t' + str(round((total - retained) / total * 100, 2) if total else 'NA')
                report = 'LINEAR AMPLIFICATION DUPLICATION METRICS\nReads before filtering\t%d\nLA duplicates removed (n)\t%d\n%s\nUnique reads after LA duplicate removal\t%d' % (total, total - retained, mqc_line, retained)
                report += '\nCount unit\tComplete primary templates; R1 alignment anchors are not proven original molecules'
                report += '\nTemplate identity\tQNAME, RG-or-absent, CB-or-absent'
                report += '\nDuplicate key\tCB-or-absent, reference, R1 strand, R1 alignment-boundary 5-prime anchor'
                report += '\nInput SAM SHA256\t' + first_hash.hexdigest()
                counts['dropped_nonprimary_records'] = counts['input_nonprimary_records'] - counts['retained_nonprimary_records']
                for label, value in sorted(counts.items()):
                    report += '\n%s\t%s' % (label, value)
                for target, text in ((metrics, report + '\n'), (mqc, header.read_text() + mqc_line + '\n')):
                    with tempfile.NamedTemporaryFile(mode='w', dir=target.parent, prefix='.' + target.name, delete=False) as handle:
                        staged.append((Path(handle.name), target))
                        handle.write(text)
                for temporary, target in staged:
                    os.replace(str(temporary), str(target))
            finally:
                for temporary, _ in staged:
                    if temporary.exists():
                        temporary.unlink()
            return dict(counts)
        finally:
            db.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('sam_path', 'output_path', 'metrics_path', 'header_path', 'mqc_path'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--max_pending', type=int, default=100000)
    parser.add_argument('--tmpdir')
    args = parser.parse_args()
    select_templates(args.sam_path, args.output_path, args.metrics_path, args.header_path,
                     args.mqc_path, args.max_pending, args.tmpdir)


if __name__ == '__main__':
    main()
