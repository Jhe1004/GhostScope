"""Compact terminal summaries and standards-compliant, lossless JSON exports."""
from pathlib import Path
import csv
import json
import math
import numpy as np

def clean(value):
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not math.isfinite(value):return 'Inf' if value>0 else ('-Inf' if value<0 else None)
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    return value

def flatten(value,prefix=''):
    for key,v in value.items():
        name=prefix+'.'+str(key) if prefix else str(key)
        if isinstance(v,dict):yield from flatten(v,name)
        else:yield name,json.dumps(clean(v),ensure_ascii=False) if isinstance(v,(list,tuple)) else str(clean(v))

def card(r):
    title='GhostScope 1.0.0'; lines=[title,'='*len(title)]
    if 'evidence' in r:
        e=r['evidence'];lines += [f"Evidence    Tier {e['tier']} — {e['title']}",f"Loci        {r.get('n_trees',r.get('n_loci','NA'))}"]
        tri=r.get('triplet',{})
        if tri:lines.append(f"Triplet     a={tri['a']:.5g}  b={tri['b']:.5g}  wedge={tri['open_wedge']}")
        q=r.get('four_leaf',{})
        if 'estimate' in q:lines.append(f"Rooted-4    Delta={q['estimate'][0]:.5g}  I1={q['estimate'][1]:.5g}  I5={q['estimate'][2]:.5g}  Rdom={q['dominance_ratio']:.3g}")
        q=r.get('six_leaf',{})
        if 'n_candidates' in q:lines.append(f"Unrooted-6  signs={q['n_sign_consistent']}/70  excluded={q['n_significant_excluded']}/70  p_union={q['intersection_union_p']:.3g}")
        q=r.get('distance_calibration',{})
        if q.get('status')=='COMPUTED':lines.append(f"Depth       raw={q['Delta_raw']:.6g}  B-clock={q['Delta_B_clock']:.6g}")
        q=r.get('cross_c',{})
        if q:
            groups=q['species_binomial']['ordered'];lines.append(f"Cross-C     {q['open_wedge_rows']}/{q['rows']} wedge rows; species pairs switched {groups['binary_label_switch_pairs']}/{groups['ghost_ever_pairs']}")
        lines.append('Compatible  '+'; '.join(e['compatible_candidate_histories']))
        lines.extend('Flag        '+f for f in e['flags'])
        lines.extend(['Next        '+e['next_step'],'Scope       '+e['scope']])
        if r.get('fixture_scope'):lines.append('Fixture     '+r['fixture_scope'])
    elif 'rows' in r and r.get('analysis')=='exploratory_scan':
        lines.append(f"{'Recipient':<24} {'Tier':<6} Interpretation")
        for row in r['rows']:lines.append(f"{row['recipient']:<24} {row['tier']:<6} {row['result']['evidence']['title']}")
        lines.append('Scope       '+r['scope'])
    elif 'n_controls' in r:
        lines += [f"BPP suite   {r['n_controls']} controls",'Directory   '+r['suite'],'Next        '+r['next_command'],'Scope       '+r['scope']]
    elif r.get('analysis')=='demo':
        for name,item in r['cases'].items():lines += ['\n'+name,card(item)]
    elif r.get('analysis')=='paper_replay':
        lines += [f"Checks      {r['passed']}/{r['total']} PASS",'Status      '+r['status'],'Details     '+r['details']]
    else:
        lines.append('Source      '+r.get('source_kind','PAIRED_NETWORKS'))
        if 'stable_operational' in r:
            lines.append('Evidence    Tier '+('D — Conditional model audit; histories remain compatible' if r['stable_operational'] else 'E — Operational convergence not established'))
        elif r.get('flags'):
            lines.append('Evidence    Tier D — Multiple histories / stage-specific diagnostics')
        if 'stable_operational' in r:lines.append(f"Convergence operational={r['stable_operational']}  strict={r['stable_strict']}  samples/chain={r['saved_samples_per_chain']}")
        for window in ('cumulative','recent_half'):
            if r.get(window,{}).get('reason'):lines.append('Window      '+window+': '+r[window]['reason'])
        if 'rare_tail_split_pairs' in r:lines.append('Tail splits '+str(r['rare_tail_split_pairs']))
        for model,entry in r.get('networks',{}).items():lines.append(f"{model:<12} stable_operational={entry.get('stable_operational')}")
        for p in r.get('comparisons',[]):lines.append(f"lnL delta   ghost - {p['competitor']}: {p['delta_mean_lnL_ghost_minus_competitor']}")
        for pair in r.get('paired',[]):
            if pair['flags']:
                lines.append('Paired      '+str(pair['group'])+' '+','.join(pair['flags']))
                for p in pair['comparisons']:
                    lines.append(f"lnL delta   ghost - {p['competitor']}: {p['delta_mean_lnL_ghost_minus_competitor']}")
        lines.extend('Flag        '+f for f in r.get('flags',[]))
        lines.extend(['Compatible  Multiple candidate histories remain unless independently compared.',
                      'Next        Inspect both windows, rare-tail visits and all TI nodes before model selection.'])
        if r.get('scope'):lines.append('Scope       '+r['scope'])
    return '\n'.join(lines)

def export(r,outdir):
    p=Path(outdir);p.mkdir(parents=True,exist_ok=True)
    (p/'report.json').write_text(json.dumps(clean(r),ensure_ascii=False,allow_nan=False,indent=2)+'\n')
    with (p/'report.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['field','value']);w.writerows(flatten(r))
    (p/'report.md').write_text('# GhostScope audit\n\n```text\n'+card(r)+'\n```\n\nFull statistics and assumptions: `report.json`.\n')
