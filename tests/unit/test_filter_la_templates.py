import importlib.util
import hashlib
from unittest.mock import patch
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('la_filter', ROOT / 'bin/filter_la_templates.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
HEADER = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n'


def pair(name='a', rg=None, cb=None, mapqs=(42, 42)):
    tags = ''.join('\t{}:Z:{}'.format(k,v) for k,v in [('RG',rg),('CB',cb)] if v is not None)
    return ['{}\t{}\tchr1\t{}\t{}\t5M\t=\t{}\t{}\tAAAAA\tIIIII{}\n'.format(name,flag,pos,q,mpos,tlen,tags) for flag,pos,q,mpos,tlen in [(99,10,mapqs[0],30,25),(147,30,mapqs[1],10,-25)]]


class FilterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def run_filter(self, original, accepted, **kwargs):
        self.original = self.root/'original.sam'; self.original.write_text(HEADER+''.join(original))
        self.accepted = self.root/'accepted.sam'; self.accepted.write_text(HEADER+''.join(accepted))
        self.output = self.root/'output.sam'; self.report = self.root/'report.json'
        result = module.select_pairs(self.original,self.accepted,self.output,self.report,tmpdir=self.root,**kwargs)
        return result['counts'], [l+'\n' for l in self.output.read_text().splitlines() if not l.startswith('@')]

    def test_report_hashes_exact_inputs(self):
        self.run_filter(pair(),pair())
        import json
        report=json.loads(self.report.read_text())
        for name in ('original','accepted'):
            self.assertEqual(report['sources'][name]['sha256'],hashlib.sha256(getattr(self,name).read_bytes()).hexdigest())
            self.assertEqual(report['sources'][name]['before'],report['sources'][name]['after'])

    def test_mutation_during_window_fails_and_cleans(self):
        real=module._select_original
        def mutate(handle,out,db,counts,pending,max_pending):
            real(handle,out,db,counts,pending,max_pending)
            with self.accepted.open('a') as stream:stream.write('@CO\tconcurrent replacement\n')
        with patch.object(module,'_select_original',mutate):
            with self.assertRaisesRegex(ValueError,'Input changed'):
                self.run_filter(pair(),pair())
        self.assertFalse(self.output.exists());self.assertFalse(self.report.exists())

    def test_reference_dictionary_mismatch_fails(self):
        # Change the parsed original dictionary through the actual input bytes,
        # before either source is opened, rather than mocking validation itself.
        original=self.root/'original.sam';accepted=self.root/'accepted.sam'
        original.write_text(HEADER+''.join(pair()))
        accepted.write_text(HEADER.replace('LN:10000','LN:9999')+''.join(pair()))
        with self.assertRaisesRegex(ValueError,'@SQ dictionaries differ'):
            module.select_pairs(original,accepted,self.root/'out.sam',self.root/'report.json')
        self.assertFalse((self.root/'out.sam').exists())

    def test_complete_template_fates_and_conservation(self):
        a,b,c,d = pair('a'),pair('b'),pair('c'),pair('d')
        counts, output = self.run_filter(a+b+c+d[:1], a+b[:1]+d[:1])
        self.assertEqual(output,a)
        for name in ['retained_templates','rejected_one_mate_failed_read_filters','rejected_both_mates_failed_read_filters','incomplete_original_templates_rejected']:
            self.assertEqual(counts[name],1)
        self.assertEqual(counts['original_primary_records'],2*counts['complete_original_templates']+counts['incomplete_original_templates_rejected'])
        self.assertIn('SO:unsorted',self.output.read_text())

    def test_identity_separates_rg_and_cb(self):
        a,b,c = pair('same','r1','c1'),pair('same','r2','c1'),pair('same','r1','c2')
        counts, output = self.run_filter(a+b+c,a+b[:1]+c)
        self.assertEqual(output,a+c); self.assertEqual(counts['retained_templates'],2)

    def test_nonprimary_and_unpaired_excluded(self):
        a=pair(); secondary=a[0].replace('\t99\t','\t355\t'); unpaired=a[0].replace('\t99\t','\t0\t')
        counts, output=self.run_filter(a+[secondary,unpaired],a+[secondary])
        self.assertEqual(output,a); self.assertEqual(counts['original_nonprimary_records'],1)
        self.assertEqual(counts['original_unpaired_records'],1)

    def test_spill_keeps_exact_records(self):
        a,b,c=pair('a'),pair('b'),pair('c')
        records=[a[0],b[0],c[0],c[1],b[1],a[1]]
        counts,output=self.run_filter(records,records,max_pending=1)
        self.assertEqual(sorted(output),sorted(records));self.assertGreater(counts['pending_spill_events'],0)
        self.assertEqual(list(self.root.glob('la-pair-filter-*')),[])

    def test_empty_is_valid(self):
        counts,output=self.run_filter([],[])
        self.assertEqual(output,[]);self.assertEqual(counts['retained_templates'],0)

    def test_foreign_accepted_record_cleans_output(self):
        with self.assertRaisesRegex(ValueError,'absent from original'):
            self.run_filter(pair(),pair()+pair('foreign'))
        self.assertFalse(self.output.exists());self.assertFalse(self.report.exists())

    def test_mutated_accepted_record_fails(self):
        a=pair()
        with self.assertRaisesRegex(ValueError,'differs'):
            self.run_filter(a,[a[0].replace('\t42\t','\t43\t'),a[1]])
        self.assertFalse(self.output.exists())

    def test_duplicate_original_mate_fails(self):
        a=pair()
        with self.assertRaisesRegex(ValueError,'Duplicate original'):
            self.run_filter([a[0],a[0]],[])

    def test_repeated_original_template_fails(self):
        a=pair()
        with self.assertRaisesRegex(ValueError,'Repeated completed'):
            self.run_filter(a+a,[])

    def test_duplicate_accepted_fails(self):
        a=pair()
        with self.assertRaisesRegex(ValueError,'Duplicate accepted'):
            self.run_filter(a,a+a[:1])

    def test_existing_products_are_preserved(self):
        self.run_filter(pair(),pair())
        before=self.output.read_bytes()
        with self.assertRaises(FileExistsError):
            module.select_pairs(self.original,self.accepted,self.output,self.report)
        self.assertEqual(self.output.read_bytes(),before)

    def test_invalid_identity_tag_fails(self):
        a=pair(rg='x')
        with self.assertRaisesRegex(ValueError,'unique nonempty'):
            self.run_filter([a[0].rstrip('\n')+'\tRG:Z:y\n',a[1]],[])

    @unittest.skipUnless(shutil.which('samtools'),'samtools CLI unavailable')
    def test_real_mapq_filter_drops_whole_asymmetric_template(self):
        a,b=pair('asymmetric',mapqs=(42,9)),pair('pass')
        source=self.root/'source.sam';source.write_text(HEADER+''.join(a+b))
        accepted=subprocess.check_output(['samtools','view','-h','-q','10','-f','1','-F','12',str(source)],text=True)
        counts,output=self.run_filter(a+b,[l+'\n' for l in accepted.splitlines() if not l.startswith('@')])
        self.assertEqual(output,b);self.assertEqual(counts['rejected_one_mate_failed_read_filters'],1)


if __name__=='__main__':unittest.main()
