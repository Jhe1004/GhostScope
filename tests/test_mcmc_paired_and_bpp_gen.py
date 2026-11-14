import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from ghostscope.mcmc_engine import *
from ghostscope.bpp_generator import *

PACKAGE=Path(__file__).resolve().parents[1]

class MCMCTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=PACKAGE/'tests');self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def traces(self,n=2000,offset=0):
        rng=np.random.default_rng(21);paths=[]
        for i in range(2):
            a=np.c_[np.arange(n),rng.normal(-100+i*offset,1,n),rng.normal(.1,.01,n),rng.normal(.08,.01,n),rng.normal(.2,.02,n)]
            p=self.p/('chain'+str(i)+'.txt');np.savetxt(p,a,header='Gen lnL theta:R tau:R phi:2<-1:H<-R',comments='');paths.append(p)
        return paths
    def test_dual_windows_and_large_mass_ratios(self):
        r=audit_chains(self.traces());self.assertTrue(r['stable_operational']);self.assertTrue(r['stable_strict'])
        self.assertTrue(r['cumulative']['rare_tail'][0]['epsilon_audits'][1]['decision_agreement'])
    def test_reference_numerical_replay(self):
        path=PACKAGE.parent/('08_Higher_Order_4to6_and_ArbitraryN_R8' if (PACKAGE.parent/'08_Higher_Order_4to6_and_ArbitraryN_R8').exists() else 'round8_456_boundaries_20261108')/'tools/dbpp_convergence.py'
        spec=importlib.util.spec_from_file_location('read_only_reference',path);ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
        paths=self.traces();x=ref.diagnostics(paths);y=audit_chains(paths)['cumulative']
        self.assertEqual(x['max_rhat'],y['max_rhat']);self.assertEqual(x['min_ess_bulk'],y['min_ess_bulk'])
        traces=[ref.read_trace(p) for p in paths];mapping={str(p):(h,a[len(a)//2:]) for p,(h,a) in zip(paths,traces)}
        ref.read_trace=lambda p:mapping[str(p)]
        x=ref.diagnostics(paths);y=audit_chains(paths)['recent_half']
        self.assertEqual(x['max_rhat'],y['max_rhat']);self.assertEqual(x['min_ess_bulk'],y['min_ess_bulk'])
    def test_different_lengths_rejected(self):
        paths=self.traces();lines=paths[1].read_text().splitlines();paths[1].write_text('\n'.join(lines[:-1])+'\n')
        with self.assertRaises(ValueError):audit_chains(paths)
    def test_shifted_chains_fail(self):self.assertFalse(audit_chains(self.traces(offset=5))['stable_operational'])
    def test_constant_parameter_fails(self):
        paths=self.traces()
        for p in paths:
            h,a=read_trace(p);a[:,2]=1.;np.savetxt(p,a,header=' '.join(h),comments='')
        self.assertFalse(audit_chains(paths)['stable_operational'])
    def test_zero_hits_bound_is_labelled_iid_only(self):
        r=rare_tail([0,76],[500000,500000]);self.assertEqual(r['critical_counts'],[50,50]);self.assertFalse(r['decision_agreement'])
        self.assertEqual(r['interval_mass_ratio'][0],float('inf'));self.assertIn('independent',r['cp_scope'])
    def test_paired_graph_trap(self):
        a=audit_chains(self.traces());r=paired_audit({'ghost':a,'inflow':a})
        self.assertIn('SINGLE_GRAPH_BF_ILLUSION_DETECTED',r['flags'])
    def test_msci_roundtrip_both_parent_semantics(self):
        t=load_tree('((A,B),C);')
        for m in ('ghost','inflow','outflow'):
            s=msci_network(t,'A','B','C',m);self.assertTrue(validate_msci(s,['A','B','C'],m))
    def test_arbitrary_backbone_keeps_all_species(self):
        t=load_tree('((((A,B),C),D),E);')
        for m in ('ghost','inflow','outflow'):
            s=msci_network(t,'A','B','D',m);self.assertTrue(validate_msci(s,list('ABCDE'),m))
    def test_generate_all_18_and_refuse_overwrite(self):
        data=self.p/'d.bpp';data.write_text('3 4\na AAAA\nb AAAC\nc AACC\n'*2);imap=self.p/'imap';imap.write_text('a A\nb B\nc C\n')
        r=make_bpp('((A,B),C);','A','B','C',data,imap,2,self.p/'suite')
        self.assertEqual(r['n_controls'],18)
        manifest=json.loads((self.p/'suite/MANIFEST.json').read_text());self.assertEqual(len({t['seed'] for t in manifest['controls']}),18)
        self.assertTrue(manifest['inputs_normalized_for_BPP']);self.assertTrue((self.p/'suite/normalized.bpp').exists())
        with self.assertRaises(ValueError):make_bpp('((A,B),C);','A','B','C',data,imap,2,self.p/'suite')
        with self.assertRaises(ValueError):make_bpp('((A,B),C);','A','B','C',data,imap,1,self.p/'invalid_single_locus_suite')

if __name__=='__main__':unittest.main()
