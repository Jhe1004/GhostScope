"""Distance decomposition and constructive moment compatibility, never a history verdict."""
import itertools
import numpy as np
from scipy import stats
from scipy.optimize import least_squares
from .topology_engine import interval_summary, DEFAULT_SEED

def bootstrap_mean(x,reps=2000,seed=DEFAULT_SEED):
    x=np.asarray(x,dtype=float)
    if x.ndim==1:x=x[:,None]
    if len(x)<2 or not np.isfinite(x).all():raise ValueError('Need at least two finite locus observations')
    rng=np.random.default_rng(seed);boots=[]
    for _ in range(reps):boots.append(x[rng.integers(len(x),size=len(x))].mean(0))
    return interval_summary(x.mean(0),np.asarray(boots))

def moment_vector(d):
    d=np.asarray(d,dtype=float)
    if d.ndim!=2 or d.shape[1]!=3 or len(d)<2 or not np.isfinite(d).all():raise ValueError('Six moments require at least two finite AB,BC,AC observations')
    return np.r_[d.mean(0),d.var(0,ddof=1)]

def six_moment_fit(d=None,means=None,variances=None,k=None,refine=True):
    """Lemma S1.1 positive construction with a common finite-epoch distance hazard.

    Failure of this sufficient construction is NOT proof of infeasibility of every
    pure-tree rate model. No unconstrained independent exponential epochs are used.
    """
    if d is not None:
        target=moment_vector(d);m=target[:3];v=target[3:]
    else:
        m=np.asarray(means,dtype=float);v=np.asarray(variances,dtype=float);target=np.r_[m,v]
    if m.shape!=(3,) or v.shape!=(3,) or not np.isfinite(target).all() or min(target)<=0:
        return {'status':'DEGENERATE_MOMENTS','compatible':False,'rejection_proved':False}
    cap=min(*m,*np.sqrt(v))
    ks=[float(k)] if k is not None else list(cap*np.r_[np.geomspace(1e-5,.95,48),.1,.5])
    def record(paths,w,common):
        pair=paths[:,[0,1,0]]+paths[:,[1,2,2]]
        pm=w@pair;pred=np.r_[pm+common,w@((pair-pm)**2)+common**2]
        residual=float(np.max(abs(pred-target)));epoch=float(paths[:,2].min()/2)
        return {'status':'FEASIBLE_POSITIVE_MSC','compatible':True,'rejection_proved':False,
                'weights':w.tolist(),'class_paths_A_B_Ctilde':paths.tolist(),'k_T':float(common),'k_Q':float(common),
                'epoch_u':epoch,'C_paths':(paths[:,2]-epoch).tolist(),'predicted_moments':pred.tolist(),
                'target_moments':target.tolist(),'max_abs_residual':residual,
                'physical_bridge':'Common positive AB/root distance hazard and positive finite epoch u; residual C_i=ctilde_i-u>0',
                'scope':'Six marginal moments only; no claim of full joint or sequence-law equality'}
    for common in ks:
        if not 0<common<cap:continue
        e=m-common;g=np.sqrt(v-common*common)
        a=np.array([(e[0]+e[2]-e[1])/2,(e[0]+e[1]-e[2])/2,(e[1]+e[2]-e[0])/2])
        if min(a)<=0:continue
        for signs in itertools.product([1,-1],repeat=3):
            z=g*np.array(signs);s=np.array([(z[0]+z[2]-z[1])/2,(z[0]+z[1]-z[2])/2,(z[1]+z[2]-z[0])/2])
            lo=max([0.]+[-s[i]/a[i] for i in range(3) if s[i]<0]);hi=min([a[i]/s[i] for i in range(3) if s[i]>0]+[np.inf])
            if not lo<hi:continue
            h=(lo+hi)/2 if np.isfinite(hi) else max(1.,lo*2)
            paths=a[None,:]+np.array([-h,1/h])[:,None]*s
            if paths.min()<=0:continue
            w=np.array([1,h*h])/(1+h*h);r=record(paths,w,common)
            if r['max_abs_residual']<1e-12:return r
    if refine:
        # Bounded nonlinear least squares restricted to the SAME realizable MSC family.
        scale=max(m.max(),np.sqrt(v).max());normalizer=np.r_[np.maximum(m,scale*1e-5),np.maximum(v,scale**2*1e-5)]
        def decode(x):return x[:6].reshape(2,3),np.array([x[6],1-x[6]]),x[7]
        def residual(x):
            paths,w,common=decode(x);pair=paths[:,[0,1,0]]+paths[:,[1,2,2]];pm=w@pair
            return (np.r_[pm+common,w@((pair-pm)**2)+common*common]-target)/normalizer
        for fraction in (.2,.5,.8):
            x=np.r_[np.tile(m.mean()/4,6),fraction,cap/10]
            lo=np.r_[np.full(6,scale*1e-12),1e-10,scale*1e-12];hi=np.r_[np.full(6,np.inf),1-1e-10,cap*(1-1e-10)]
            fit=least_squares(residual,x,bounds=(lo,hi),max_nfev=1200,ftol=1e-13,xtol=1e-13,gtol=1e-13)
            paths,w,common=decode(fit.x);r=record(paths,w,common)
            if r['max_abs_residual']<1e-12:return r
    return {'status':'NO_CONSTRUCTIVE_SOLUTION','compatible':False,'rejection_proved':False,
            'target_moments':target.tolist(),'scope':'Search failed in a sufficient subfamily; broader pure-tree rate models remain unexcluded'}

def calibrate(d,labels,reps=2000,seed=DEFAULT_SEED,outgroup_distances=None,resolved_fraction=1.):
    d=np.asarray(d,dtype=float);labels=np.asarray(labels)
    star=np.c_[(d[:,0]+d[:,2]-d[:,1])/2,(d[:,0]+d[:,1]-d[:,2])/2,(d[:,1]+d[:,2]-d[:,0])/2]
    if np.min(star)<-1e-10:raise ValueError('Distances violate the triangle inequality')
    x=star[labels==0];y=star[labels==1]
    if min(len(x),len(y))<2:return {'status':'INSUFFICIENT_TOPOLOGIES','quasi_clock_pass':False}
    def calc(a,b):
        aa=a.mean(0);bb=b.mean(0);db=aa[1]-bb[1];ac=aa[0]-bb[2]
        return np.array([db+ac,db,ac,2*db,aa[0]/aa[1] if aa[1]>0 else np.nan,bb[2]/bb[1] if bb[1]>0 else np.nan])
    point=calc(x,y);rng=np.random.default_rng(seed)
    boots=np.array([calc(x[rng.integers(len(x),size=len(x))],y[rng.integers(len(y),size=len(y))]) for _ in range(reps)])
    names=['Delta_raw','Delta_B','Delta_AC','Delta_B_clock','rA_over_rB','rC_over_rB']
    result={k:float(v) for k,v in zip(names,point)}
    result['bootstrap']={k:dict(zip(('estimate','se','ci95','p_one_sided','p_two_sided'),
                                         (s[0] for s in [q['estimate'],q['se'],q['ci95'],q['p_one_sided'],q['p_two_sided']])))
                         for k,i in zip(names,range(6)) for q in [interval_summary(point[i:i+1],boots[:,i:i+1])]}
    fraction=abs(point[2])/(abs(point[1])+abs(point[2])) if abs(point[1])+abs(point[2]) else None
    of=None
    if outgroup_distances is not None:
        od=np.asarray(outgroup_distances);den=od[:,1].mean()
        of=float(abs((od[:,0]-od[:,1]).mean())/den) if den>0 else None
    checks={'path_ratio_symmetry':bool(np.isfinite(point[4:]).all() and np.max(abs(point[4:]-1))<=.1),
            'small_conditional_path_offset':fraction is not None and fraction<=.2,
            'outgroup_contrast_small':None if of is None else of<=.05,
            'resolved_fraction_healthy':resolved_fraction>=.5}
    quasi=all(v is True for v in checks.values())
    cb=result['bootstrap']['Delta_B_clock'];positive=bool(cb['ci95'][0]>0 and cb['p_one_sided']<.01)
    result.update(status='COMPUTED',counts_12_23=[len(x),len(y)],Delta_AC_fraction=float(abs(point[2]/point[0])) if point[0] else None,
                  conditional_path_offset_fraction=fraction,outgroup_contrast=of,quasi_clock_checks=checks,
                  quasi_clock_pass=quasi,quasi_clock_support=bool(quasi and positive),
                  sign_flip=bool(point[0]>0 and point[3]<0),
                  significance_lost=bool(result['bootstrap']['Delta_raw']['p_two_sided']<.05 and cb['p_two_sided']>=.05),
                  identity_residual=float(abs(point[0]-point[1]-point[2])),
                  scope='Conditional path ratios are not identified mutation rates. Quasi-clock checks are heuristic consistency screens; they cannot rule out ancestral structure, rate mixtures, or outgroup violations.')
    return result

def multi_haplotype_distances(matrices,names,haplotypes,b,c,reps=2000,seed=DEFAULT_SEED,latent=False,exchangeable=False):
    indices=[names.index(h) for h in haplotypes];ib=names.index(b);ic=names.index(c);values=[]
    if len(indices)<2:raise ValueError('Need at least two recipient haplotypes')
    for d in matrices:
        pairs=[((d[i,ib]-d[j,ib])**2,(d[i,ic]-d[j,ic])**2) for i,j in itertools.combinations(indices,2)]
        z=np.mean(pairs,axis=0);values.append([z[0],z[1],z[0]-z[1]])
    result=bootstrap_mean(values,reps,seed);negative=result['ci95'][2][1]<0
    result.update(status='SUFFICIENT_ZONE_CANDIDATE' if negative else 'SCOPE_ROUTING',
                  code='MULTI_HAPLOTYPE_SUFFICIENT_ZONE' if negative else 'MULTI_HAPLOTYPE_IN_CONFOUNDED_CONE',
                  tier_b_eligible=bool(negative and latent and exchangeable),n=len(values),
                  scope='Tier B requires same-species exchangeable recipient copies, common per-locus rate, finite second moments, latent genealogical distances and the specified canonical candidate models. Finite-sequence estimates are proxies.')
    return result

def sequence_multi_rows(seq,haplotypes,b,c,outgroup=None):
    names=list(haplotypes)+[b,c]+([outgroup] if outgroup else [])
    if not set(names)<=seq.keys():return None
    x=np.stack([np.frombuffer(seq[n].encode(),dtype=np.uint8) for n in names]);x=x[:,np.all(np.isin(x,list(b'ACGT')),axis=0)];n=x.shape[1]
    if n<2:return None
    bseq=x[len(haplotypes)];cseq=x[len(haplotypes)+1];pairs=[]
    for i,j in itertools.combinations(range(len(haplotypes)),2):
        zb=(x[i]!=bseq).astype(float)-(x[j]!=bseq).astype(float);zc=(x[i]!=cseq).astype(float)-(x[j]!=cseq).astype(float)
        ib=zb.mean()**2;ic=zc.mean()**2
        ub=(zb.sum()**2-(zb*zb).sum())/(n*(n-1));uc=(zc.sum()**2-(zc*zc).sum())/(n*(n-1))
        pairs.append([ib,ic,ib-ic,ub,uc,ub-uc])
    return np.mean(pairs,axis=0)

def survival_curvature(distances,window=None,reps=2000,seed=DEFAULT_SEED):
    d=np.asarray(distances);d=d[np.isfinite(d)]
    if len(d)<10:return {'status':'INSUFFICIENT_OBSERVATIONS'}
    lo,hi=np.quantile(d,[.2,.8]) if window is None else window
    if not hi>lo:return {'status':'DEGENERATE_WINDOW'}
    points=np.linspace(lo,hi,6);survival=np.array([(d>p).mean() for p in points])
    if min(survival)<=0:return {'status':'ZERO_SURVIVAL','window':[lo,hi]}
    hankel=np.array([[survival[i+j] for j in range(3)] for i in range(3)])
    # Six equally spaced survival observations cannot reject an arbitrary <=3-mode
    # family by a 3x3 determinant alone; expose the moments without overclaiming.
    return {'status':'EXPLORATORY_ONLY','points':points.tolist(),'survival':survival.tolist(),
            'hankel_singular_values':np.linalg.svd(hankel,compute_uv=False).tolist(),
            'scope':'Population separation requires a certified common breakpoint-free window and mode constraints. This finite-sample screen has no calibrated test or automatic evidence tier.'}

def site_patterns(seq,names):
    if not set(names)<=seq.keys():return None
    x=np.stack([np.frombuffer(seq[n].encode(),dtype=np.uint8) for n in names]);x=x[:,np.all(np.isin(x,list(b'ACGT')),axis=0)]
    if x.shape[1]<2:return None
    a,b,c,o=x
    return np.array([np.sum((a==b)&(c==o)&(a!=c)),np.sum((a==o)&(b==c)&(a!=b)),
                     np.sum((b==o)&(a==c)&(a!=b)),np.sum((b==o)&(a!=o)),np.sum((a==o)&(b!=o)),x.shape[1]])

def summarize_patterns(rows):
    x=np.asarray(rows);s=x.sum(0);aa,bb,cc,ua,ub,n=s;den=bb+cc
    if not den:return {'status':'NO_ABBA_BABA','called_sites':int(n)}
    d=(bb-cc)/den;leave=s-x;valid=leave[:,1]+leave[:,2]>0
    if len(x)>1 and valid.all():
        ds=(leave[:,1]-leave[:,2])/(leave[:,1]+leave[:,2]);se=np.sqrt((len(x)-1)/len(x)*np.sum((ds-ds.mean())**2))
    else:se=np.nan
    chi=(ua-ub)**2/(ua+ub) if ua+ub else 0.
    return {'status':'COMPUTED','BBAA_ABBA_BABA':[int(aa),int(bb),int(cc)],'called_sites':int(n),'D':float(d),
            'locus_jackknife_se':float(se),'p_D':float(2*stats.norm.sf(abs(d/se))) if se>0 else None,
            'gprime_pattern':float((bb-cc)/(aa+bb-2*cc)) if aa+bb!=2*cc else None,
            'tajima_unique_A_B':[int(ua),int(ub)],'tajima_chi2':float(chi),'tajima_p':float(stats.chi2.sf(chi,1)),
            'scope':'Site-pattern screen. Gamma is a pattern ratio, not a reimplementation of the official HyDe variance test. D uses deposited loci as jackknife blocks; physical linkage not inferred. Tajima detects cumulative path differences.'}
