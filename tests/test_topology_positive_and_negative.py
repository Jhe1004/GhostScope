from collections import Counter
from importlib.resources import files
import json
import unittest
from ghostscope.io_parsers import parse_newick,unroot_key
from ghostscope.topology_engine import *
from ghostscope.grading import grade
from ghostscope.cli import demo

class TopologyTests(unittest.TestCase):
    def test_wedge_routes_and_does_not_reject(self):
        trees=[parse_newick('((A:1,B:1):1,C:2);')]*60+[parse_newick('((B:1,C:1):1,A:2);')]*30+[parse_newick('((A:1,C:1):1,B:2);')]*10
        r=triplet(trees,['A','B','C']);self.assertTrue(r['open_wedge']);self.assertEqual(r['delta_log_mpl'],0.)
        e=grade({'triplet':r});self.assertEqual(e['tier'],'D')
    def test_six_certificate_really_has_eight_term_cell(self):
        rows=six_table()['contrasts'];self.assertEqual(len(rows),70)
        special=next(r for r in rows if r['cell']==[0,'a','e'])
        self.assertEqual(len(special['terms']),8)
    def test_unroot_key_rotation(self):
        a=parse_newick('(((a:1,h:1):1,b:1):1,(c:1,(d:1,e:1):1):1);')
        b=parse_newick('(c:1,((d:1,e:1):1,(b:1,(h:1,a:1):1):2):0);')
        self.assertEqual(unroot_key(a),unroot_key(b))
    def test_sparse_loci_no_spurious_seventy_exclusions(self):
        trees=[parse_newick('(((a:1,h:1):1,b:1):1,(c:1,(d:1,e:1):1):1);')]*100
        r=six_leaf(trees,list('ahbcde'),reps=100)
        self.assertFalse(r['passed']);self.assertLess(r['n_significant_excluded'],70)
    def test_parallel_is_conditional_warning(self):
        r=spine(parse_newick('((((A,H),B),C),(D,E));'),'A')
        self.assertTrue(r['type_III_possible']);self.assertEqual(grade({'spine':r})['tier'],'D')
    def test_rooted_left_right_canonical(self):
        self.assertEqual(canonical(parse_newick('(((d,b),c),a);')),canonical(parse_newick('(a,(c,(b,d)));')))
    def test_mpl_closure_is_not_claimed_as_finite_fit(self):
        r=profile([0,1,9]);self.assertEqual(r['status'],'CLOSURE_SUPREMUM_ONLY')
        self.assertEqual(r['finite_maximizer_histories'],[])
    def test_six_failure_cannot_be_overridden_by_four_projection(self):
        r={'n_ingroup':6,'four_leaf':{'passed':True},'six_leaf':{'passed':False}}
        self.assertEqual(grade(r)['tier'],'D')
    def test_demo_real_calculation_both_directions(self):
        r=demo();pos=r['cases']['Positive topology control'];neg=r['cases']['Ambiguous depth control']
        self.assertEqual(pos['evidence']['tier'],'A');self.assertEqual(pos['six_leaf']['n_significant_excluded'],70)
        self.assertEqual(neg['evidence']['tier'],'E');self.assertTrue(neg['distance_calibration']['sign_flip'])
        self.assertIn('SINGLE_GRAPH_BF_ILLUSION_DETECTED',r['cases']['Paired-network trap']['flags'])

if __name__=='__main__':unittest.main()
