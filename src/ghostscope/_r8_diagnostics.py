# Adapted from project dbpp_convergence.py; formulas retained for numerical replay.
"""Outcome-blind, versioned strict and operational diagnostics from real draws.

Operational thresholds are the supplied experiment policy, not a claim that all
posterior quantities are converged. Every column and failed rule is retained.
"""
import math,re
from pathlib import Path
import numpy as np
from scipy.stats import rankdata,norm
POLICY={'version':'R8_TWO_TIER_20261108_v1','eps':.01,'bf_cutoff':100.,'active_mean_phi':.02,
 'strict_rhat':1.05,'strict_ess':200,'decision_rhat':1.05,'active_phi_ess':150,'lnL_ess':150,
 'active_or_backbone_rhat':1.08,'active_or_backbone_ess':100,'all_rhat':1.10,
 'dormant_phi_ess_floor':None,'constant_free_parameter':'FAIL_BOTH_GATES',
 'rate_hyperparameters':'Apply active/backbone thresholds; never treat as dormant'}
def read_trace(p):
 import gzip
 p = Path(p)
 if not p.exists() and p.with_name(p.name + '.gz').exists(): p = p.with_name(p.name + '.gz')
 opener = (lambda x: gzip.open(x, 'rt')) if p.suffix == '.gz' else (lambda x: Path(x).open())
 with opener(p) as f: header=f.readline().split(); a=np.loadtxt(f,ndmin=2)
 if a.shape[1]!=len(header) or not np.isfinite(a).all():raise ValueError('Malformed/nonfinite trace')
 return header,a
def ess_bulk(x):
 m,n=x.shape;centered=x-x.mean(axis=1,keepdims=True)
 ac=np.fft.irfft(abs(np.fft.rfft(centered,n=2*n,axis=1))**2,axis=1)[:,:n]/n
 W=ac[:,0].mean();B=n*np.var(x.mean(axis=1),ddof=1);vp=(n-1)/n*W+B/n
 if not vp>0:return 0.0
 rho=1-(W-ac.mean(axis=0))/vp;total=0.;previous=float('inf')
 for j in range(0,n-1,2):
  pair=float(rho[j]+rho[j+1])
  if pair<0:break
  pair=min(pair,previous);total+=pair;previous=pair
 return min(m*n,m*n/max(1,-1+2*total))
def diagnostics(paths=None, traces=None):
 if traces is None:traces=[read_trace(p) for p in paths]
 if len(traces)<2:raise ValueError('At least two independent chains required')
 header=traces[0][0]
 if any(t[0]!=header for t in traces):raise ValueError('Chain column mismatch')
 lengths=[len(a) for _,a in traces]
 if len(set(lengths))!=1:raise ValueError('Different chain lengths: refuse silent truncation')
 n=lengths[0];h=n//2
 if h<20:return {'accepted':False,'strict_all_columns_converged':False,'operational_decision_converged':False,'reason':'Too few samples','samples_per_chain':n,'parameters':[],'policy':POLICY}
 edge_rows=[];active_nodes=set();known_nodes=set()
 for j,name in enumerate(header):
  if not name.startswith(('phi:','phi_')):continue
  pair=[a[:,j] for _,a in traces]
  if any(np.any((z<0)|(z>1)) for z in pair):raise ValueError('phi outside [0,1]')
  means=[float(z.mean()) for z in pair];counts=[int(np.sum(z<.01)) for z in pair];supported=[.01/(c/n)>=100 if c else True for c in counts]
  active=any(supported) or any(m>=.02 for m in means);label=name.rsplit(':',1)[-1]
  nodes=[x.strip() for x in label.split('<-')];mapped=len(nodes)==2
  if mapped:known_nodes.update(nodes)
  if active and mapped:active_nodes.update(nodes)
  edge_rows.append({'parameter':name,'scenario':label,'nodes':nodes if mapped else [],'active':active,'means_by_chain':means,'posterior_mean_phi':float(np.mean(means)),
   'count_below_eps_by_chain':counts,'samples_by_chain':lengths,'b10_eps0.01_by_chain':['Inf' if not c else .01*n/c for c in counts],
   'supported_by_chain':supported,'decision_agreement':len(set(supported))==1,'mapping_known':mapped})
 lookup={r['parameter']:r for r in edge_rows};records=[]
 for j,name in enumerate(header[1:],1):
  raw=np.stack([a[:2*h,j] for _,a in traces]);x=raw.reshape(2*len(traces),h)
  constant=np.ptp(x)==0
  if constant:rh=float('inf');ess=0.
  else:
   def rhat(z):
    W=np.var(z,axis=1,ddof=1).mean();B=h*np.var(z.mean(axis=1),ddof=1)
    return float(np.sqrt(((h-1)*W/h+B/h)/W)) if W>0 else float('inf')
   def rank(z):return norm.ppf((rankdata(z.ravel()).reshape(z.shape)-.375)/(z.size+.25))
   z=rank(x);rh=max(rhat(z),rhat(rank(abs(x-np.median(x)))));ess=ess_bulk(z)
  strict=bool(not constant and rh<=1.05 and ess>=200)
  if name in lookup:
   active=lookup[name]['active'];role='active_phi' if active else 'dormant_phi';limit=1.05;floor=150 if active else None
  elif name=='lnL':role='log_likelihood';limit=1.05;floor=150
  elif name.startswith(('tau:','theta:')):
   node=name.rsplit(':',1)[-1]
   if re.fullmatch(r'Z\d+',node) and node in known_nodes and node not in active_nodes:role='dormant_auxiliary';limit=1.10;floor=None
   else:role='active_auxiliary' if node in active_nodes else ('unmapped_auxiliary_active' if re.fullmatch(r'Z\d+',node) else 'backbone');limit=1.08;floor=100
  else:role='other_active_hyperparameter';limit=1.08;floor=100
  passed=bool(not constant and rh<=limit and (floor is None or ess>=floor));failed=[]
  if constant:failed.append('CONSTANT_FREE_TRACE')
  if rh>limit:failed.append('RHAT_LIMIT')
  if floor is not None and ess<floor:failed.append('ESS_LIMIT')
  records.append({'parameter':name,'constant':bool(constant),'rhat':rh,'ess_bulk':ess,'accepted':strict,'strict_pass':strict,'operational_role':role,'operational_rhat_limit':limit,'operational_ess_floor':floor,'operational_pass':passed,'operational_failures':failed})
 agreement=all(x['decision_agreement'] for x in edge_rows);mapping=all(x['mapping_known'] for x in edge_rows)
 strict=all(r['strict_pass'] for r in records) and agreement;oper=all(r['operational_pass'] for r in records) and agreement and mapping and all(r['rhat']<=1.10 for r in records)
 return {'accepted':bool(oper),'legacy_all_columns_converged':all(r['strict_pass'] for r in records),'strict_all_columns_converged':bool(strict),
 'operational_decision_converged':bool(oper),'edge_decision_agreement':agreement,'auxiliary_mapping_complete':mapping,'samples_per_chain':n,'parameters':records,'candidate_edges':edge_rows,
 'max_rhat':max((r['rhat'] for r in records),default=0),'min_ess_bulk':min((r['ess_bulk'] for r in records),default=0),'policy':POLICY,
 'criterion':'Both gates retained; accepted means operational decision gate, not whole-posterior convergence; no truth labels used'}
