"""Topology contrasts, with candidate-specific nulls and explicit finite-data scope."""
from collections import Counter, defaultdict
from importlib.resources import files
import csv
import json
import numpy as np
from scipy.stats import norm, chi2
from .io_parsers import induced, canonical, binary, unroot_key
from ._mpl import profile

DEFAULT_SEED = 20261111

def node(*children):
    return tuple(sorted(children,key=repr))

def bootstrap_categories(keys, values, reps=2000, seed=DEFAULT_SEED):
    if reps < 100: raise ValueError('At least 100 bootstrap replicates are required')
    counts = Counter(keys); n = sum(counts.values())
    if n < 2: raise ValueError('At least two loci required')
    categories = list(counts)
    probabilities = np.array([counts[k]/n for k in categories])
    weights = np.array([values(k) for k in categories],dtype=float)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n, probabilities,size=reps).astype(float)
    boots = np.einsum('ij,jk->ik',draws,weights,optimize=False)/n
    return np.einsum('i,ij->j',probabilities,weights,optimize=False), boots

def interval_summary(point, boots):
    point = np.asarray(point); se = np.std(boots,axis=0,ddof=1)
    z = np.divide(point,se,out=np.zeros_like(point,dtype=float),where=se>0)
    nonzero=(se==0)&(point!=0)
    z[nonzero]=np.sign(point[nonzero])*np.inf
    return {'estimate':point.tolist(),'se':se.tolist(),'ci95':np.quantile(boots,[.025,.975],axis=0).T.tolist(),
            'p_one_sided':norm.sf(z).tolist(),'p_two_sided':(2*norm.sf(abs(z))).tolist(),
            'bootstrap_nonpositive_fraction':np.mean(boots<=0,axis=0).tolist()}

def triplet(trees, names):
    a,b,c = names
    index = {node(node(a,b),c):0,node(node(b,c),a):1,node(node(a,c),b):2}
    counts = np.zeros(3,dtype=int); labels=[]; missing=0
    for t in trees:
        if not set(names)<=set(t.tips()): missing+=1; labels.append(-2); continue
        label = index.get(canonical(induced(t,names)),-1)
        labels.append(label)
        if label>=0: counts[label]+=1
    n=int(counts.sum()); p=counts/n if n else np.zeros(3)
    aa=float(p[0]-p[2]); bb=float(p[1]-p[2]); wedge=bool(n and aa>0 and bb>0)
    return {'taxa':names,'counts_AB_BC_AC':counts.tolist(),'resolved':n,'missing':missing,
            'unresolved':len(trees)-n-missing,'probabilities_AB_BC_AC':p.tolist(),'a':aa,'b':bb,
            'open_wedge':wedge,'status':'SCOPE_ROUTING' if wedge else 'OUTSIDE_EMPIRICAL_OPEN_WEDGE',
            'gamma_open_interval':[bb,1-aa] if wedge else None,
            'delta_log_mpl':0. if n else None,'MPL':profile(counts),
            'gprime_ISM':bb/(aa+bb) if wedge else None,
            'hyde_orientation':'gprime_ISM corresponds to 1-Gamma for the paper ordering; verify software order',
            'mpl_scope':'Exact equality of profiled topology likelihoods on the specified open wedge; not a finite-DNA likelihood result',
            'labels':labels}

def four_leaf(trees,names,reps=2000,seed=DEFAULT_SEED,gtee_floor=.015):
    if gtee_floor<0: raise ValueError('GTEE floor cannot be negative')
    a,d,b,c=names
    delta={node(node(node(d,b),c),a):1,node(node(node(a,b),c),d):-1}
    i1={node(a,node(d,node(b,c))):1,node(a,node(b,node(c,d))):-1}
    terms=[(-1,node(b,node(d,node(a,c)))),(-1,node(c,node(a,node(b,d)))),
           (1,node(c,node(d,node(a,b)))),(-1,node(d,node(c,node(a,b)))),(1,node(node(a,c),node(b,d)))]
    i5=defaultdict(int)
    for s,t in terms: i5[t]+=s
    keys=[]; missing=0; unresolved=0
    for t in trees:
        if not set(names)<=set(t.tips()):missing+=1;continue
        key=canonical(induced(t,names))
        if not binary(key): unresolved+=1;continue
        keys.append(key)
    if len(keys)<2:return {'status':'INSUFFICIENT_RESOLVED_LOCI','n':len(keys),'passed':False,'missing':missing,'unresolved':unresolved}
    point,boots=bootstrap_categories(keys,lambda k:[delta.get(k,0),i1.get(k,0),i5.get(k,0)],reps,seed)
    result=interval_summary(point,boots); cov=np.cov(boots[:,1:].T,ddof=1); rank=int(np.linalg.matrix_rank(cov))
    inv=np.linalg.pinv(cov); off=point[1:]-cov@inv@point[1:]
    statistic=float(point[1:]@inv@point[1:])
    exact_p=0. if np.linalg.norm(off)>1e-10 else (float(chi2.sf(statistic,rank)) if rank else 1.)
    ci=np.asarray(result['ci95']); dominant=abs(point[0])/max(abs(point[1]),abs(point[2]),1e-9)
    overlaps=(ci[1:,0]<=gtee_floor)&(ci[1:,1]>=-gtee_floor)
    # The directive's point-estimate and CI-overlap rules disagree: require BOTH.
    tolerance=bool(np.all(abs(point[1:])<=gtee_floor) and np.all(overlaps))
    positive=bool(ci[0,0]>0 and result['p_one_sided'][0]<.01)
    passed=bool(positive and tolerance and dominant>=3.)
    return dict(result,status='GTEE_HEURISTIC_SUPPORT' if passed else 'NO_FOUR_LEAF_SIGNATURE',
                passed=passed,exact_zero_p=exact_p,exact_zero_df=rank,exact_zero_statistic=statistic,
                exact_zero_nonsignificant=exact_p>=.01,gtee_floor=gtee_floor,tolerance_pass=tolerance,
                dominance_ratio=float(dominant),names=names,n=len(keys),missing=missing,unresolved=unresolved,
                conditioning='Resolved induced rooted trees only; unresolved/missing loci reported separately',
                validity='GTEE floor is a user-specified tolerance, not an estimated error bound or calibrated false-positive rate',
                candidate_scope='Fixed rooted ladder; specified single-pulse direct/private-route competitors; not every introgression network')

def six_table():
    return json.loads(files('ghostscope').joinpath('data/six_leaf_70_contrasts.json').read_text())

def _freeze(x):
    return tuple(_freeze(v) for v in x) if isinstance(x,list) else x

def six_leaf(trees,names,reps=2000,seed=DEFAULT_SEED):
    standard=six_table(); mapping=dict(zip(names,standard['standard_tips'])); keys=[];missing=0;unresolved=0
    for tree in trees:
        if not set(names)<=set(tree.tips()):missing+=1;continue
        sub=induced(tree,names); key=unroot_key(sub)
        if len(key)!=3:unresolved+=1;continue
        key=tuple(sorted(tuple(sorted(mapping[x] for x in split)) for split in key))
        keys.append(key)
    return six_leaf_counts(Counter(keys),names,reps,seed,missing,unresolved)

def six_leaf_counts(counts,names,reps=2000,seed=DEFAULT_SEED,missing=0,unresolved=0):
    """Same statistic from compact raw topology multiplicities (used by the demo)."""
    standard=six_table();n=sum(counts.values())
    if n<2:return {'status':'INSUFFICIENT_RESOLVED_LOCI','passed':False,'n':n,'missing':missing,'unresolved':unresolved}
    contrasts=standard['contrasts']
    weights=[{_freeze(t['coordinate']):t['coefficient'] for t in r['terms']} for r in contrasts]
    point,boots=bootstrap_categories(counts,lambda k:[w.get(k,0) for w in weights],reps,seed)
    summary=interval_summary(point,boots);se=np.asarray(summary['se']);signs=np.array([r['target_sign'] for r in contrasts])
    z=np.divide(point,se,out=np.zeros(70),where=se>0)
    nonzero=(se==0)&(point!=0)
    z[nonzero]=np.sign(point[nonzero])*np.inf
    signed_z=z*signs; ps=norm.sf(signed_z); consistent=point*signs>0; excluded=(signed_z>=1.96)&consistent
    rows=[dict(cell=r['cell'],contrast=float(point[i]),se=float(se[i]),z=float(z[i]),target_sign=int(signs[i]),
               p_directional=float(ps[i]),ci95=summary['ci95'][i],sign_consistent=bool(consistent[i]),excluded=bool(excluded[i])) for i,r in enumerate(contrasts)]
    weakest=int(np.argmax(ps)); passed=bool(excluded.all())
    return {'status':'ALL_DIRECT_CELLS_EXCLUDED' if passed else 'DIRECT_CANDIDATES_REMAIN',
            'passed':passed,'names':names,'n':n,'missing':missing,'unresolved':unresolved,
            'n_sign_consistent':int(consistent.sum()),'n_significant_excluded':int(excluded.sum()),'n_candidates':70,
            'intersection_union_p':float(max(ps)),'weakest_link':rows[weakest],'contrasts':rows,
            'candidate_scope':standard['scope'],
            'finite_data_scope':'Exact nulls apply to latent genealogies; inferred gene trees can introduce systematic error. No GTEE correction for these 70 contrasts.',
            'test_logic':'Tier A requires every candidate-specific null rejected in the certified direction; high sign count alone is insufficient'}

def spine(root,recipient,outgroups=()):
    keep=set(root.tips())-set(outgroups)
    if recipient not in keep:raise ValueError('Recipient absent from ingroup species tree')
    tree=induced(root,keep); path=[]
    def find(n):
        if not n.children:return n.name==recipient
        for c in n.children:
            if find(c):
                path.append([x for x in n.children if x is not c]);return True
        return False
    find(tree)
    clades=[sorted(x for n in siblings for x in n.tips()) for siblings in path]
    binary_path=all(len(siblings)==1 for siblings in path)
    # A multi-tip sister clade is a sampling opportunity, not an identified donor.
    parallel=any(len(c)>=2 for c in clades[1:])
    return {'recipient':recipient,'sibling_clades':clades,'clade_sizes':[len(c) for c in clades],
            'm_eff':len(clades)+1,'binary_spine':binary_path,'terminal_ladder':binary_path and all(len(c)==1 for c in clades),
            'order':[recipient]+[c[0] for c in clades if len(c)==1],
            'type_III_possible':parallel,
            'warning':'TYPE_III_PARALLEL_DONOR_WARNING' if parallel else None,
            'fiber_lower_bound_conditional':2*len(keep)+4,
            'fiber_scope':'Independent positive branch rates, time gauge, terminal single recipient/pulse and private singleton segments; topology alone does not establish these assumptions',
            'type_III_scope':'Potential ambiguity only: exact equivalence additionally requires rate-path coupling, positive chronology and the specified stem/crown construction'}

def cross_c_table(path,reference_taxa=()):
    with open(path,encoding='utf-8') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows:raise ValueError('Empty cross-C table')
    required={'Taxon_A','Taxon_B','Taxon_C','AB_count','BC_count','AC_count','Test_conclusion'}
    if not required<=rows[0].keys():raise ValueError('Unsupported cross-C header; expected ghostparser columns')
    groups=defaultdict(list);unordered=defaultdict(list); interior=0; ghostrows=0;ghostinside=0
    refnames=sorted(reference_taxa,key=len,reverse=True)
    for r in rows:
        ab,bc,ac=(float(r[k]) for k in ('AB_count','BC_count','AC_count'))
        inside=ab>ac and bc>ac;interior+=int(inside)
        ghost='unsampled' in r['Test_conclusion'].lower();ghostrows+=int(ghost);ghostinside+=int(ghost and inside)
        groups[(r['Taxon_A'],r['Taxon_B'])].append((r['Taxon_C'],ghost,r['Test_conclusion']))
        unordered[tuple(sorted((r['Taxon_A'],r['Taxon_B'])))].append((r['Taxon_C'],ghost,r['Test_conclusion']))
    def summarize(gs):
        selected={k:v for k,v in gs.items() if any(x[1] for x in v)}
        changing={k:v for k,v in selected.items() if any(not x[1] and any(x[0]!=g[0] for g in v if g[1]) for x in v)}
        return {'tested_pairs':len(gs),'ghost_ever_pairs':len(selected),'binary_label_switch_pairs':len(changing),
                'switching_rate':len(changing)/len(selected) if selected else None}
    def level(mapper):
        ordered=defaultdict(list);unordered=defaultdict(list);unmapped=set()
        for row in rows:
            vals=[mapper(row['Taxon_'+k]) for k in 'ABC']
            if not all(vals):
                unmapped.update(row['Taxon_'+k] for k,v in zip('ABC',vals) if v is None);continue
            if len(set(vals))<3:continue
            ghost='unsampled' in row['Test_conclusion'].lower();value=(vals[2],ghost,row['Test_conclusion'])
            ordered[tuple(vals[:2])].append(value);unordered[tuple(sorted(vals[:2]))].append(value)
        return {'ordered':summarize(ordered),'unordered':summarize(unordered),'unmapped':sorted(unmapped)}
    result={'rows':len(rows),'open_wedge_rows':interior,'ghost_label_rows':ghostrows,'ghost_label_open_wedge_rows':ghostinside,
            'ordered':summarize(groups),'unordered':summarize(unordered),
            'scope':'Author-oriented table labels; changing C and taxon orientation is sensitivity, not independent biological validation'}
    result['species_binomial']=level(lambda name:'_'.join(name.split('_')[:2]))
    if refnames:result['reference_taxon']=level(lambda name:next((n for n in refnames if name==n or name.startswith(n+'_')),None))
    return result
