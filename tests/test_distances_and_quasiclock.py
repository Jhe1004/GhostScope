import unittest
import numpy as np
from ghostscope.distance_engine import *
from ghostscope.grading import grade

class DistanceTests(unittest.TestCase):
    def test_positive_physical_six_moment_reconstruction(self):
        paths=np.array([[.012,.006,.025],[.02,.015,.03]]);w=np.array([.3,.7]);k=.002
        pair=paths[:,[0,1,0]]+paths[:,[1,2,2]];mean=w@pair
        result=six_moment_fit(means=mean+k,variances=w@((pair-mean)**2)+k*k,k=k)
        self.assertTrue(result['compatible']);self.assertLess(result['max_abs_residual'],1e-15)
        self.assertGreater(min(result['C_paths']),0)
    def test_failed_fit_not_universal_rejection(self):
        r=six_moment_fit(means=[1,10,1],variances=[.001]*3,refine=False)
        self.assertFalse(r['compatible']);self.assertFalse(r['rejection_proved'])
    def test_depth_identity_and_sign_flip(self):
        # AB star A=.02 B=.005 C=.035; BC star A=.04 B=.006 C=.01.
        x=np.tile([.025,.04,.055],(60,1));y=np.tile([.046,.016,.05],(40,1))
        r=calibrate(np.r_[x,y],np.r_[np.zeros(60),np.ones(40)],reps=200)
        self.assertAlmostEqual(r['identity_residual'],0);self.assertTrue(r['sign_flip'])
    def test_missing_outgroup_cannot_pass_all_four_checks(self):
        rng=np.random.default_rng(1);d=np.r_[rng.normal(.1,.001,(50,3)),rng.normal(.09,.001,(50,3))]
        r=calibrate(d,np.r_[np.zeros(50),np.ones(50)],reps=100)
        self.assertIsNone(r['quasi_clock_checks']['outgroup_contrast_small']);self.assertFalse(r['quasi_clock_pass'])
    def test_sequence_multi_hap_proxy_not_tier_b(self):
        seq={'X':'AAAAACCCCC','Y':'CCCCCAAAAA','B':'AAAACCCCCC','C':'AAAAACCCCC'}
        row=sequence_multi_rows(seq,['X','Y'],'B','C');self.assertEqual(len(row),6)
        self.assertEqual(grade({'multi_haplotype':{'tier_b_eligible':False,'code':'SUFFICIENT_ZONE_CANDIDATE'}})['tier'],'D')
    def test_latent_multi_hap_requires_assumptions(self):
        names=['X','Y','B','C'];d=np.zeros((4,4));d[0,3]=d[3,0]=1
        r=multi_haplotype_distances([d]*100,names,['X','Y'],'B','C',reps=100,latent=True,exchangeable=False)
        self.assertFalse(r['tier_b_eligible'])
        r=multi_haplotype_distances([d]*100,names,['X','Y'],'B','C',reps=100,latent=True,exchangeable=True)
        self.assertTrue(r['tier_b_eligible'])
    def test_survival_is_not_automatic_theorem(self):
        r=survival_curvature(np.linspace(.01,1,100));self.assertEqual(r['status'],'EXPLORATORY_ONLY')
        self.assertEqual(len(r['survival']),6)
    def test_site_mask_and_pattern_orientation(self):
        r=site_patterns({'A':'ACAN','B':'CAAN','C':'CCAN','O':'AAAN'},['A','B','C','O'])
        self.assertEqual(r[5],3);s=summarize_patterns([r,r]);self.assertIn('gprime_pattern',s)

if __name__=='__main__':unittest.main()
