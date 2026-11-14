"""Draw-level and explicitly labelled summary-level competing-network audits."""
from pathlib import Path
import csv
import json
import math
import re
import numpy as np
from scipy.stats import beta
from ._r8_diagnostics import diagnostics, read_trace, ess_bulk, POLICY

RATIO_SCOPE = ('epsilon / estimated posterior mass(phi<epsilon), assuming Uniform(0,1). '
               'A finite-interval mass ratio, not an exact point-null Savage-Dickey BF or a between-network BF. '
               'Zero hits do not establish zero probability. Chain means of lnL are descriptive, not marginal likelihoods.')

def rare_tail(counts,lengths,eps=.01,cutoff=100.):
    if not 0<eps<1 or cutoff<=0:raise ValueError('Invalid epsilon or cutoff')
    if len(counts)!=len(lengths) or not counts:raise ValueError('Missing chain counts')
    if any(n<=0 or c<0 or c>n for c,n in zip(counts,lengths)):raise ValueError('Invalid rare-tail counts')
    ratios=[float(eps*n/c) if c else float('inf') for c,n in zip(counts,lengths)]
    limits=[math.floor(eps*n/cutoff) for n in lengths]
    decisions=[c<=k for c,k in zip(counts,limits)]
    upper=[float(beta.ppf(.975,c+1,n-c)) if c<n else 1. for c,n in zip(counts,lengths)]
    return {'epsilon':eps,'counts':list(map(int,counts)),'samples':list(map(int,lengths)),
            'critical_counts':limits,'interval_mass_ratio':ratios,'supported_by_chain':decisions,
            'decision_agreement':len(set(decisions))==1,
            'iid_only_cp95_mass_upper':upper,'iid_only_ratio_lower':[eps/u for u in upper],
            'cp_scope':'Clopper-Pearson coverage assumes independent Bernoulli draws; these bounds are not conservative guarantees for correlated MCMC. Indicator ESS and temporal visits must be examined.',
            'scope':RATIO_SCOPE}

def _window(traces):
    base=diagnostics(traces=traces); header=traces[0][0];n=len(traces[0][1]);extras=[]
    if n<40:return base
    for j,name in enumerate(header):
        if not name.startswith(('phi:','phi_')):continue
        pair=[a[:,j] for _,a in traces];epsrows=[]
        for eps in (.005,.01,.02):
            counts=[int((x<eps).sum()) for x in pair];r=rare_tail(counts,[n]*len(pair),eps)
            indicators=np.array([x<eps for x in pair],dtype=float)
            r['indicator_ess']=float(ess_bulk(indicators)) if not all(c==0 for c in counts) else None
            r['counts_in_10_time_blocks']=[[int(v.sum()) for v in np.array_split(x<eps,10)] for x in pair]
            epsrows.append(r)
        extras.append({'parameter':name,'epsilon_audits':epsrows,
                       'quantile_ci95_by_chain':[np.quantile(x,[.025,.975]).tolist() for x in pair]})
    continuous=bool(base.get('parameters') and all(r['operational_pass'] for r in base['parameters']))
    split=any(not r['epsilon_audits'][1]['decision_agreement'] for r in extras)
    base.update(continuous_parameter_gate=continuous,rare_tail=extras,
                flags=['BULK_CONVERGED_RARE_TAIL_SPLIT'] if continuous and split else [],
                odd_final_sample_excluded_from_rank_split=n%2,scope=RATIO_SCOPE)
    return base

def audit_chains(paths):
    if len(paths)<2:raise ValueError('At least two independent chain files required')
    traces=[read_trace(p) for p in paths]
    if len({len(a) for _,a in traces})!=1:raise ValueError('Different chain lengths: refuse silent truncation')
    if len({tuple(h) for h,_ in traces})!=1:raise ValueError('Chain column mismatch')
    n=len(traces[0][1]);full=_window(traces);recent=_window([(h,a[n//2:]) for h,a in traces])
    header=traces[0][0];lnl=[float(a[:,header.index('lnL')].mean()) for _,a in traces] if 'lnL' in header else None
    return {'source_kind':'RAW_TRACES','chains':[str(p) for p in paths],'saved_samples_per_chain':n,
            'cumulative':full,'recent_half':recent,
            'stable_operational':bool(full.get('accepted') and recent.get('accepted')),
            'stable_strict':bool(full.get('strict_all_columns_converged') and recent.get('strict_all_columns_converged')),
            'lnL_means':lnl,'flags':sorted(set(full.get('flags',[])+recent.get('flags',[]))),
            'policy':POLICY,'scope':RATIO_SCOPE}

def _network_summary(raw):
    if 'cumulative' in raw:
        edges=raw['cumulative'].get('candidate_edges',[])
        return {'accepted':raw['stable_operational'],'edges':edges,'lnL_means':raw.get('lnL_means')}
    return {'accepted':raw.get('stable_operational',raw.get('reported_cumulative_operational',False)),
            'dual_window':raw.get('stable_operational',False),'edges':raw.get('phi',[]),'lnL_means':raw.get('lnL_means')}

def paired_audit(networks):
    summaries={k:_network_summary(v) for k,v in networks.items()};comparisons=[];flags=[]
    ghost=summaries.get('ghost')
    def supports(s):return bool(s and s['accepted'] and s['edges'] and all(all(e.get('supported_by_chain',[])) and bool(e.get('supported_by_chain')) for e in s['edges']))
    if ghost:
        for name,s in summaries.items():
            if name=='ghost':continue
            both=supports(ghost) and supports(s)
            if both:flags.append('SINGLE_GRAPH_BF_ILLUSION_DETECTED')
            valid_ll=bool(ghost['lnL_means'] and s['lnL_means'] and all(v is not None for v in ghost['lnL_means']+s['lnL_means']))
            delta=float(np.mean(ghost['lnL_means'])-np.mean(s['lnL_means'])) if valid_ll else None
            comparisons.append({'competitor':name,'both_graphs_pass_with_large_interval_ratios':both,
                                'both_dual_window_gates_verified':bool(ghost.get('dual_window',ghost['accepted']) and s.get('dual_window',s['accepted'])),
                                'delta_mean_lnL_ghost_minus_competitor':delta})
    return {'networks':networks,'comparisons':comparisons,'flags':sorted(set(flags)),
            'same_data_requirement':'Comparisons presume identical loci, sites, likelihood scaling and observations; verify input provenance. No automatic between-network evidence from lnL.',
            'scope':RATIO_SCOPE}

def audit_summary(path):
    """Recompute flags from recorded numbers; never pretend to re-read missing draws."""
    p=Path(path)
    if p.suffix=='.tsv':
        with p.open() as f:table=list(csv.DictReader(f,delimiter='\t'))
        if not table or not {'id','BF10_eps01_chain1','BF10_eps01_chain2','actual_saved_samples_per_chain'}<=table[0].keys():
            raise ValueError('Unsupported audit TSV; supply the R8 chain-summary schema or raw traces')
        rows=[]
        def number(x):
            try:return float(x)
            except (ValueError,TypeError):return None
        for row in table:
            n=int(number(row['actual_saved_samples_per_chain']) or 0);ratios=[number(row['BF10_eps01_chain'+str(i)]) for i in (1,2)]
            r={k:row.get(k) for k in ('id','system','model','mode')}
            r.update(loci=int(number(row.get('loci')) or 0),partition=row.get('partition'),samples_config=row.get('samples'),
                     reported_cumulative_operational=row.get('Tier2')=='1',phi=[],
                     lnL_means=[number(row.get('lnL_chain'+str(i))) for i in (1,2)],
                     gate_scope='Imported initial TSV; dual-window diagnostics not available; hit counts reconstructed only if reported ratios determine integers')
            if n and all(v is not None and v>0 for v in ratios):
                inferred=[0. if math.isinf(v) else .01*n/v for v in ratios]
                if any(abs(v-round(v))>1e-5 for v in inferred):raise ValueError('Rounded TSV ratios do not identify integer tail hits; use exact JSON or raw traces')
                counts=[int(round(v)) for v in inferred];a=rare_tail(counts,[n,n]);means=[number(row.get('phi_chain'+str(i))) for i in (1,2)]
                r['phi']=[{'parameter':row.get('parameter'),'means_by_chain':means,'posterior_mean_phi':number(row.get('phi_pooled')),
                           'count_below_eps_by_chain':counts,'samples_by_chain':[n,n],'supported_by_chain':a['supported_by_chain'],
                           'count_provenance':'Reconstructed from recorded interval ratios and sample count; not newly observed trace visits'}]
            rows.append(r)
        result=_audit_summary_rows(rows);result['source_kind']='IMPORTED_R8_TSV_NOT_RAW_DRAWS'
        return result
    obj=json.loads(p.read_text());rows=obj if isinstance(obj,list) else obj.get('tasks',obj.get('summary',obj.get('rows',[]))) if isinstance(obj,dict) else []
    if not isinstance(rows,list) or not rows:raise ValueError('Expected audit JSON with a summary list')
    return _audit_summary_rows(rows)

def _audit_summary_rows(rows):
    outputs=[];groups={}
    for r in rows:
        if 'tier2' in r and 'stable_operational' not in r:
            r=dict(r,reported_cumulative_operational=bool(r['tier2']),
                   gate_scope='Initial cumulative-window audit only; recent-half diagnostics absent from this saved product')
        audits=[];flags=[]
        for edge in r.get('phi',[]):
            counts=edge.get('count_below_eps_by_chain');lengths=edge.get('samples_by_chain')
            if not counts or not lengths:continue
            a=rare_tail(counts,lengths);audits.append(dict(parameter=edge.get('parameter'),**a))
        continuous=bool(r.get('cumulative_all_parameter_columns_strict') and r.get('recent_all_parameter_columns_strict'))
        if continuous and any(not a['decision_agreement'] for a in audits):flags.append('BULK_CONVERGED_RARE_TAIL_SPLIT')
        out=dict(r);out.update(recomputed_rare_tail=audits,flags=flags);outputs.append(out)
        key=(r.get('system'),r.get('loci'),r.get('mode'),r.get('partition'),r.get('samples_config'))
        groups.setdefault(key,{})[r.get('model',r.get('id'))]=out
    paired=[dict(group=list(k),**paired_audit(v)) for k,v in groups.items() if 'ghost' in v and len(v)>1]
    return {'source_kind':'SAVED_AUDIT_SUMMARY_NOT_RAW_DRAWS','rows':outputs,'paired':paired,'rate_model_shifts':rate_shifts(outputs),
            'flags':sorted(set(f for r in outputs+paired for f in r['flags'])),
            'rare_tail_split_pairs':sum('BULK_CONVERGED_RARE_TAIL_SPLIT' in r['flags'] for r in outputs),
            'scope':'Recomputed counts/decisions from saved summaries. Original convergence diagnostics and means are reported, not independently recomputed from unavailable 500k trajectories. '+RATIO_SCOPE}

def rate_shifts(rows):
    groups={};out=[]
    for r in rows:
        key=(r.get('system'),r.get('loci'),r.get('model'),r.get('partition'),r.get('samples_config','single'))
        groups.setdefault(key,{})[r.get('mode')]=r
    for key,modes in groups.items():
        if 'default' not in modes:continue
        ref={r['parameter']:r for r in modes['default'].get('phi',[])}
        for mode,r in modes.items():
            if mode=='default':continue
            for edge in r.get('phi',[]):
                name=edge['parameter']
                if name not in ref:continue
                a=ref[name].get('posterior_mean_phi');b=edge.get('posterior_mean_phi')
                if a is None or b is None:continue
                out.append({'group':list(key),'parameter':name,'mode':mode,'reference_mode':'default',
                            'mean_phi_shift':b-a,
                            'reference_gate':_network_summary(modes['default'])['accepted'],'mode_gate':_network_summary(r)['accepted'],
                            'scope':'Descriptive sensitivity of fixed-graph posteriors; retain stage-specific gate qualifications.'})
    return out

def audit_summaries(paths):
    """Merge initial then continuation records by ID, never pool different stages."""
    if len(paths)==1:return audit_summary(paths[0])
    records={};sources=[]
    for path in paths:
        r=audit_summary(path);sources.append({'path':str(path),'source_kind':r['source_kind']})
        for row in r.get('rows',[]):
            key=row.get('id')
            if not key:raise ValueError('Summary merge needs task IDs')
            records[key]=row
    groups={}
    for r in records.values():
        key=(r.get('system'),r.get('loci'),r.get('mode'),r.get('partition'),r.get('samples_config','single'))
        groups.setdefault(key,{})[r.get('model',r['id'])]=r
    pairs=[dict(group=list(k),**paired_audit(v)) for k,v in groups.items() if 'ghost' in v and len(v)>1]
    rows=list(records.values())
    return {'source_kind':'MERGED_STAGE_SUMMARIES_NOT_RAW_DRAWS','sources':sources,'rows':rows,'paired':pairs,'rate_model_shifts':rate_shifts(rows),
            'rare_tail_split_pairs':sum('BULK_CONVERGED_RARE_TAIL_SPLIT' in r['flags'] for r in rows),
            'flags':sorted(set(f for r in rows+pairs for f in r['flags'])),
            'scope':'Later files replace earlier task records by ID. Initial cumulative eligibility is not dual-window convergence. '+RATIO_SCOPE}

def ti_audit(directories):
    if len(directories)!=2:raise ValueError('TI requires two independently run node directories')
    sets=[]
    for folder in directories:
        found={}
        for p in Path(folder).glob('*.mcmc.txt'):
            match=re.search(r'\.b(\d+)\.mcmc\.txt$',p.name)
            if match:found[int(match.group(1))]=p
        sets.append(found)
    if not sets[0] or sets[0].keys()!=sets[1].keys():raise ValueError('TI node sets are missing or differ')
    rows=[]
    for key in sorted(sets[0]):
        r=audit_chains([s[key] for s in sets]);rows.append({'node':key,**r})
    accepted=all(r['stable_strict'] for r in rows)
    return {'source_kind':'TI_NODE_RAW_TRACES','nodes':rows,'all_nodes_strict':accepted,
            'scope':'Every supplied node is checked in both windows. Complete quadrature weights and endpoint coverage are additionally required; this command does not invent a marginal likelihood from unweighted node means.'}
