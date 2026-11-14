"""Read-only regression fixtures from the local research project.

These tests test numerical/behavioral consistency, not a promised power of 100%.
"""
import json
from pathlib import Path
import unittest
from ghostscope.io_parsers import read_trees,parse_newick,unroot_key,load_tree
from ghostscope.analysis import audit_trees,scan_trees
from ghostscope.mcmc_engine import audit_summary
from ghostscope.bpp_generator import msci_network
from ghostscope.topology_engine import node,six_table,_freeze

ROOT=Path(__file__).resolve().parents[2]
R8=ROOT/('08_Higher_Order_4to6_and_ArbitraryN_R8' if (ROOT/'08_Higher_Order_4to6_and_ArbitraryN_R8').exists() else 'round8_456_boundaries_20261108')

class PaperRegressionTests(unittest.TestCase):
    def test_direct_outer_donor_must_not_be_called_ghost(self):
        p=R8/'ingroup_6/6I-Exp1-E_r010/genealogies.jsonl'
        r=audit_trees(list(read_trees(p)),taxa=['A1_1','A2_1','B_1','C_1','D_1','E_1'],outgroups=['Outgroup_1'],reps=200,distances=False)
        # This real counterexample can satisfy a four-leaf projected signature.
        self.assertTrue(r['four_leaf']['passed']);self.assertEqual(r['evidence']['tier'],'D')
        direct=next(c for c in r['six_leaf']['contrasts'] if c['cell']==[0,'a','e'])
        self.assertLess(abs(direct['z']),1.96)
    def test_true_ghost_500_loci_are_not_forged_into_seventy_rejections(self):
        p=R8/'ingroup_6/6I-Exp0_r001/genealogies.jsonl'
        r=audit_trees(list(read_trees(p)),taxa=['A1_1','A2_1','B_1','C_1','D_1','E_1'],reps=200,distances=False)
        self.assertEqual(r['n_input_trees'],500);self.assertLess(r['six_leaf']['n_significant_excluded'],70)
        self.assertEqual(r['evidence']['tier'],'D')
    def test_exactly_five_saved_tail_reversals(self):
        p=R8/'empirical_extension68_20261108/results/AMD_FINAL/05_results/EXTENSION68_ALL_RESULTS.json'
        r=audit_summary(p);self.assertEqual(r['rare_tail_split_pairs'],5)
        self.assertEqual(r['source_kind'],'SAVED_AUDIT_SUMMARY_NOT_RAW_DRAWS')
    def test_real_audit_tsv_is_imported_with_stage_provenance(self):
        r=audit_summary(R8/'empirical_new_boundaries_20261108/results/E4A_PLANT_BPP_DUALCHAIN_AUDIT.tsv')
        self.assertEqual(r['source_kind'],'IMPORTED_R8_TSV_NOT_RAW_DRAWS')
        self.assertTrue(r['rows']);self.assertIn('dual-window',r['rows'][0]['gate_scope'])
    def test_scan_handles_every_recipient_even_without_terminal_spine(self):
        p=R8/'ingroup_4/4I-Exp0_r003/genealogies.jsonl'
        r=scan_trees(list(read_trees(p)),'(((A1,A2),B),C);',outgroups=['Outgroup'],delimiter='_',reps=100,distances=False)
        self.assertEqual({x['recipient'] for x in r['rows']},{'A1','A2','B','C'})
    def test_original_msci_graphs_semantically_match(self):
        base=ROOT/('06_Empirical_Data_Audits_R5' if (ROOT/'06_Empirical_Data_Audits_R5').exists() else 'round5_empirical_scouting')/'round5B/extracted/data/emprical-result/Thuja/BPP'
        def signature(t):
            return (t.name if not t.children and not t.annotations else ('HYBRID' if t.annotations else ''),
                    tuple(t.annotations),tuple(sorted([signature(c) for c in t.children],key=repr)))
        # Recover species names from the historical graph instead of analytical
        # code special-casing any taxon.
        ghostpath=next((base/'ghost').glob('*.ctl'));text=ghostpath.read_text()
        original=next(line.strip() for line in text.splitlines() if 'H[&phi' in line)
        graph=parse_newick(original);hybrid=next(c for c in graph.children if c.annotations)
        recipient=hybrid.children[0].name
        species=set(graph.tips())-{'H'}
        from ghostscope.io_parsers import canonical,induced
        sibling_branch=next(c for c in graph.children if not c.annotations)
        # The hybrid leaf and its ordinary sister define the second anchor.
        def sister_of_h(n):
            if any(c.name=='H' and not c.children for c in n.children):return next(c.name for c in n.children if c.name!='H')
            for c in n.children:
                r=sister_of_h(c)
                if r:return r
        sister=sister_of_h(sibling_branch);outer=next(iter(species-{recipient,sister}))
        backbone=load_tree(f'(({recipient},{sister}),{outer});')
        for model,folder in [('ghost','ghost'),('inflow','inflow'),('outflow','ouflow')]:
            ctl=next((base/folder).glob('*.ctl'));old=next(l.strip() for l in ctl.read_text().splitlines() if 'H[&phi' in l)
            new=msci_network(backbone,recipient,sister,outer,model)
            self.assertEqual(signature(parse_newick(new)),signature(parse_newick(old)))

class CompleteCoordinateTests(unittest.TestCase):
    def test_all_945_rooted_topologies_project_to_105_certificate_coordinates(self):
        from functools import lru_cache
        @lru_cache(None)
        def trees(tips):
            if len(tips)==1:return [tips[0]]
            answer=[];first=tips[0]
            for bits in range(1<<(len(tips)-1)):
                left=(first,)+tuple(tips[i+1] for i in range(len(tips)-1) if bits&(1<<i))
                right=tuple(t for t in tips if t not in left)
                if not right:continue
                answer.extend(node(a,b) for a in trees(left) for b in trees(right))
            return answer
        def nw(t):return t if isinstance(t,str) else '('+','.join(nw(x) for x in t)+')'
        topologies=trees(tuple('ahbcde'));self.assertEqual(len(topologies),945)
        keys={unroot_key(parse_newick(nw(t)+';')) for t in topologies};self.assertEqual(len(keys),105)
        self.assertTrue(all(_freeze(t['coordinate']) in keys for c in six_table()['contrasts'] for t in c['terms']))

if __name__=='__main__':unittest.main()
