import gzip
from pathlib import Path
import tempfile
import unittest
import numpy as np
from ghostscope.io_parsers import *

PACKAGE=Path(__file__).resolve().parents[1]

class IOTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=PACKAGE/'tests');self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def test_quoted_comments_scientific_and_roundtrip(self):
        t=parse_newick("(('sample:a':1e-3[comment],B:2E-3)99:0.2,'O''name':.3);")
        self.assertEqual(set(t.tips()),{'sample:a','B',"O'name"})
        self.assertTrue(np.allclose(patristic(t,t.tips()),patristic(parse_newick(newick(t)),t.tips())))
    def test_lca_regression_not_pendant_twigs(self):
        t=parse_newick('(((A:1,X:2):3,B:4):5,(C:6,O:7):8);')
        d=patristic(t,['A','B','C']);sub=induced(t,['A','B','C'])
        self.assertEqual(d[0,1],8);self.assertEqual(d[0,2],23)
        np.testing.assert_allclose(d,patristic(sub,['A','B','C']))
    def test_root_path_invariance(self):
        t=parse_newick('((A:1,B:2):3,(C:4,O:5):6);');r=reroot_outgroup(t,['O'])
        np.testing.assert_allclose(patristic(t,t.tips()),patristic(r,t.tips()))
        self.assertEqual(canonical(induced(r,['A','B','C'])),canonical(parse_newick('((A,B),C);')))
    def test_multitip_outgroup(self):
        t=parse_newick('(((A:1,B:1):2,C:3):1,(O:1,P:1):3);')
        r=reroot_outgroup(t,['O','P']);np.testing.assert_allclose(patristic(t,t.tips()),patristic(r,t.tips()))
        with self.assertRaises(ValueError):reroot_outgroup(t,['O','A'])
    def test_polytomy_and_zero_internal_edges(self):
        t=parse_newick('((A:1,B:1):0,C:1);');self.assertFalse(binary(canonical(t)))
        q=parse_newick('((A:1,B:1):0,(C:1,D:1):0);');self.assertEqual(unroot_key(q),())
    def test_negative_length_not_distance(self):
        with self.assertRaises(ValueError):patristic(parse_newick('(A:-1,B:1);'),['A','B'])
    def test_imap_and_delimiter(self):
        (self.p/'imap').write_text('a A\nb B\n')
        m=TaxonMapper(self.p/'imap','_',{'Species_one'})
        self.assertEqual(m('prefix^a'),'A');self.assertEqual(m('C_gene1'),'C');self.assertEqual(m('Species_one'),'Species_one')
        with self.assertRaises(ValueError):TaxonMapper(delimiter='_').tree(parse_newick('(A_1,A_2,B_1);'))
    def test_fasta_gzip_and_bpp(self):
        with gzip.open(self.p/'x.fa.gz','wt') as f:f.write('>A\nACGT\n>B\nACGA\n')
        (self.p/'x.bpp').write_text('2 4\nA ACGT\nB ACGA\n\n2 3\nA AAA\nB AAC\n')
        self.assertEqual(len(list(read_alignments(self.p/'x.bpp'))),2)
        self.assertEqual(next(read_alignments(self.p/'x.fa.gz'))[1]['A'],'ACGT')
    def test_interleaved_phylip(self):
        (self.p/'x.phy').write_text('3 8\nA ACGT\nB ACGT\nC ACGA\n\nACGT\nACGA\nACGA\n')
        seq=next(read_alignments(self.p/'x.phy'))[1]
        self.assertEqual(seq['A'],'ACGTACGT')
    def test_jc_mask_and_saturation(self):
        names,d,n=distance_matrix({'A':'AANN','B':'ACGT'})
        self.assertEqual(n[0,1],2);self.assertAlmostEqual(d[0,1],-.75*np.log(1-4*.5/3))
        self.assertFalse(np.isfinite(distance_matrix({'A':'AAAA','B':'CCCC'})[1]).all())
    def test_nj_additive_quartet(self):
        t=parse_newick('((A:1,B:2):3,(C:4,D:5):2);');names=['A','B','C','D'];d=patristic(t,names)
        q=neighbor_joining(names,d)
        np.testing.assert_allclose(patristic(q,names),d)
        self.assertEqual(unroot_key(q),unroot_key(t))
    def test_nj_star_no_artificial_topology(self):
        d=np.ones((6,6))-np.eye(6);q=neighbor_joining(list('ABCDEF'),d)
        self.assertEqual(unroot_key(q),())
    def test_tree_stream_gzip(self):
        with gzip.open(self.p/'x.tre.gz','wt') as f:f.write('(A:1,B:1);\n(A:2,B:2);')
        self.assertEqual(len(list(read_trees(self.p/'x.tre.gz'))),2)
    def test_malformed_inputs(self):
        for text in ['(A,B','(A,A);','(A:NaN,B:1);','(A:,B:1);']:
            with self.subTest(text=text),self.assertRaises(ValueError):parse_newick(text)

if __name__=='__main__':unittest.main()
