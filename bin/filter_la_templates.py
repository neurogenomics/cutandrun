#!/usr/bin/env python3
"""Opt-in LA filter: retain original primary pairs only when both records pass.

The accepted SAM must come from the existing readwise SAMtools filters on the
same original BAM. This helper does not reproduce MAPQ/region logic or recover
missing mates. It intersects exact original records, then keeps whole templates.
"""
import argparse
import hashlib
from collections import Counter, OrderedDict
import json
from pathlib import Path
import sqlite3
import tempfile


def record(line):
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 11:
        raise ValueError('SAM record has fewer than 11 fields')
    flag = int(fields[1])
    if flag < 0 or flag > 65535:
        raise ValueError('Invalid SAM flag')
    if flag & 0x900:
        return None, None, 'nonprimary'
    if not flag & 1:
        return None, None, 'unpaired'
    first, second = bool(flag & 64), bool(flag & 128)
    if first == second:
        raise ValueError('Primary paired record lacks exactly one mate flag')
    tags = {}
    for field in fields[11:]:
        if field[:3] in ('RG:', 'CB:'):
            bits = field.split(':', 2)
            if len(bits) != 3 or bits[1] != 'Z' or not bits[2] or bits[0] in tags:
                raise ValueError('RG/CB must be unique nonempty Z tags')
            tags[bits[0]] = bits[2]
    if not fields[0] or fields[0] == '*':
        raise ValueError('Template identity requires QNAME')
    key = json.dumps([fields[0], tags.get('RG'), tags.get('CB')], separators=(',', ':'))
    return key, 1 if first else 2, 'primary'


def stat_identity(path):
    stat = Path(path).stat()
    return dict(device=stat.st_dev, inode=stat.st_ino, size=stat.st_size,
                mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns)


def source_lines(path, digest, dictionary):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        for line in handle:
            digest.update(line.encode('utf-8'))
            if line.startswith('@SQ\t'):
                fields = line.rstrip('\r\n').split('\t')[1:]
                tags = dict(field.split(':', 1) for field in fields)
                name = tags.get('SN')
                if not name or name in dictionary or 'LN' not in tags:
                    raise ValueError('Invalid or duplicate @SQ header entry')
                dictionary[name] = tags
            yield line


def select_pairs(original, accepted, output, report, tmpdir=None, max_pending=100000):
    if max_pending < 1:
        raise ValueError('max_pending must be positive')
    if Path(output).exists() or Path(report).exists():
        raise FileExistsError('Refusing to replace an existing output or report')
    if Path(output).resolve() == Path(report).resolve():
        raise ValueError('Output and report must be different paths')
    output_created = False
    report_created = False
    sources = {name: dict(path=str(Path(path).resolve()), before=stat_identity(path))
               for name, path in [('original', original), ('accepted', accepted)]}
    script_before = stat_identity(__file__)
    script_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    digests = {name: hashlib.sha256() for name in sources}
    dictionaries = {name: {} for name in sources}
    counts = Counter()
    pending = OrderedDict()
    with tempfile.TemporaryDirectory(prefix='la-pair-filter-', dir=tmpdir) as scratch:
        database = Path(scratch) / 'identity.sqlite'
        db = sqlite3.connect(str(database))
        try:
            db.execute('PRAGMA cache_size=-8192')
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('PRAGMA synchronous=OFF')
            db.execute('CREATE TABLE accepted (identity TEXT, mate INTEGER, record TEXT, seen INTEGER DEFAULT 0, PRIMARY KEY(identity,mate))')
            db.execute('CREATE TABLE pending (identity TEXT PRIMARY KEY, mate INTEGER, record TEXT, passed INTEGER)')
            db.execute('CREATE TABLE completed (identity TEXT PRIMARY KEY)')
            for line in source_lines(accepted, digests['accepted'], dictionaries['accepted']):
                if line.startswith('@'):
                    continue
                key, mate, kind = record(line)
                counts['accepted_' + kind + '_records'] += 1
                if kind == 'nonprimary':
                    continue
                if kind == 'unpaired':
                    raise ValueError('Accepted SAM contains unpaired primary record')
                try:
                    db.execute('INSERT INTO accepted(identity,mate,record) VALUES(?,?,?)', (key, mate, line))
                except sqlite3.IntegrityError as error:
                    raise ValueError('Duplicate accepted primary identity/mate') from error
            with Path(output).open('x') as out:
                output_created = True
                handle = source_lines(original, digests['original'], dictionaries['original'])
                _select_original(handle, out, db, counts, pending, max_pending)
            counts['incomplete_original_templates_rejected'] = len(pending) + db.execute('SELECT count(*) FROM pending').fetchone()[0]
            if db.execute('SELECT count(*) FROM accepted WHERE seen=0').fetchone()[0]:
                raise ValueError('Accepted SAM has primary records absent from original SAM')
            if dictionaries['original'] != dictionaries['accepted']:
                raise ValueError('Original and accepted @SQ dictionaries differ')
            for name, path in [('original', original), ('accepted', accepted)]:
                sources[name]['after'] = stat_identity(path)
                sources[name]['sha256'] = digests[name].hexdigest()
                if sources[name]['before'] != sources[name]['after']:
                    raise ValueError('Input changed during filtering: ' + name)
            if script_before != stat_identity(__file__):
                raise ValueError('Filter script changed during execution')
            db.commit()
            counts['sqlite_bytes'] = database.stat().st_size
            for name in ('retained_templates','complete_original_templates','rejected_one_mate_failed_read_filters','rejected_both_mates_failed_read_filters','original_primary_records','accepted_primary_records','original_nonprimary_records','accepted_nonprimary_records','original_unpaired_records','pending_spill_events'):
                counts.setdefault(name, 0)
            result = dict(schema='la-complete-template-filter/v1',status='passed',policy='Both exact original primary records must survive existing readwise filters; incomplete original templates rejected; nonprimary records excluded.',identity=['QNAME','RG_or_absent','CB_or_absent'],sources=sources,script_sha256=script_sha256,parameters=dict(max_pending=max_pending,sqlite_cache_mib=8),counts=dict(counts),limits=['Explicit candidate method change; not equivalent to historical readwise filtering.','Reciprocal mate geometry is checked by the downstream LA deduplicator.','Temporary SQLite size and full-depth runtime require benchmarking.'])
            with Path(report).open('x') as handle:
                report_created = True
                json.dump(result, handle, indent=2)
                handle.write('\n')
            return result
        except Exception:
            # Remove only files created by this call; never remove a pre-existing output.
            for path, created in ((output, output_created), (report, report_created)):
                if created:
                    try:
                        Path(path).unlink()
                    except FileNotFoundError:
                        pass
            raise
        finally:
            db.close()


def _select_original(handle, out, db, counts, pending, max_pending):
    for line in handle:
        if line.startswith('@'):
            if line.startswith('@HD\t'):
                line = '\t'.join(f for f in line.rstrip('\n').split('\t') if not f.startswith('SO:')) + '\tSO:unsorted\n'
            out.write(line)
            continue
        key, mate, kind = record(line)
        counts['original_' + kind + '_records'] += 1
        if kind != 'primary':
            continue
        accepted_row = db.execute('SELECT record FROM accepted WHERE identity=? AND mate=?', (key, mate)).fetchone()
        passed = accepted_row is not None
        if passed:
            if accepted_row[0] != line:
                raise ValueError('Accepted record differs from its original record')
            db.execute('UPDATE accepted SET seen=1 WHERE identity=? AND mate=?', (key, mate))
        if db.execute('SELECT 1 FROM completed WHERE identity=?', (key,)).fetchone():
            raise ValueError('Repeated completed primary template identity')
        prior = pending.pop(key, None)
        if prior is None:
            prior = db.execute('SELECT mate,record,passed FROM pending WHERE identity=?', (key,)).fetchone()
            if prior is not None:
                db.execute('DELETE FROM pending WHERE identity=?', (key,))
        if prior is None:
            pending[key] = (mate, line, passed)
            if len(pending) > max_pending:
                oldkey, values = pending.popitem(last=False)
                db.execute('INSERT INTO pending VALUES(?,?,?,?)', (oldkey,) + values)
                counts['pending_spill_events'] += 1
            counts['peak_pending_in_memory'] = max(counts['peak_pending_in_memory'], len(pending))
            continue
        if prior[0] == mate:
            raise ValueError('Duplicate original primary identity/mate')
        db.execute('INSERT INTO completed VALUES(?)', (key,))
        counts['complete_original_templates'] += 1
        if prior[2] and passed:
            for value in sorted([prior, (mate, line, passed)]):
                out.write(value[1])
            counts['retained_templates'] += 1
        elif prior[2] or passed:
            counts['rejected_one_mate_failed_read_filters'] += 1
        else:
            counts['rejected_both_mates_failed_read_filters'] += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('original','accepted','output','report'):
        parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--tmpdir',type=Path)
    parser.add_argument('--max-pending',type=int,default=100000)
    args = parser.parse_args()
    select_pairs(**vars(args))


if __name__ == '__main__':
    main()
