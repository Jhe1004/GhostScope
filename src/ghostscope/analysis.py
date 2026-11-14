"""High-level reusable Python analysis without filesystem side effects."""
import numpy as np
from pathlib import Path
from .io_parsers import TaxonMapper, reroot_outgroup, induced, patristic, load_tree
from .topology_engine import triplet, four_leaf, six_leaf, spine
from .distance_engine import calibrate,six_moment_fit
from .grading import grade

def audit_trees(trees,taxa=None,species_tree=None,recipient=None,outgroups=(),imap=None,delimiter=None,
                reps=2000,seed=20261111,gtee_floor=.015,rooted=False,distances=True):
    if isinstance(species_tree,(str,Path)):species_tree=load_tree(species_tree)
    known=set(species_tree.tips()) if species_tree else set(taxa or [])|set(outgroups)
    mapper=TaxonMapper(imap,delimiter,known);mapped=[];ambiguous=0;omitted=set();total_input=0
    selected=set(taxa or [])|set(outgroups) if taxa else None
    for t in trees:
        total_input+=1
        if selected:
            retained=[tip for tip in t.tips() if mapper(tip) in selected]
            omitted.update(mapper(tip) for tip in t.tips() if mapper(tip) not in selected)
            if not retained:continue
            t=induced(t,retained)
        try:mapped.append(mapper.tree(t))
        except ValueError:
            ambiguous+=1
    if ambiguous and not mapped:raise ValueError('All loci have multiple samples per mapped species; choose representatives explicitly')
    if not mapped:raise ValueError('No gene trees supplied')
    root_failures=0;processed=[]
    for t in mapped:
        if outgroups:
            if not set(outgroups)<=set(t.tips()):root_failures+=1;continue
            try:t=reroot_outgroup(t,outgroups)
            except ValueError:root_failures+=1;continue
        processed.append(t)
    if not processed:raise ValueError('No trees can be rooted with the specified monophyletic outgroup')
    nroot=bool(outgroups or rooted)
    result={'n_trees':len(processed),'n_input_trees':total_input,'outgroup_rooting_failures':root_failures,
            'rooting':'outgroup' if outgroups else ('declared_rooted' if rooted else 'unrooted'),
            'flags':[],'bootstrap_reps':reps,'seed':seed}
    result['ambiguous_multiple_sample_loci_excluded']=ambiguous
    result['omitted_sampled_taxa']=sorted(omitted)
    if species_tree:
        if not recipient:recipient=next(t for t in species_tree.tips() if t not in outgroups)
        info=spine(species_tree,recipient,outgroups);result['spine']=info
        if taxa is None:
            if info['terminal_ladder']:taxa=info['order']
            else:
                # Deterministic projection, flagged as a projection, never a full-ladder proof.
                taxa=[recipient]+[c[0] for c in info['sibling_clades']]
                result['flags'].append('NONTERMINAL_SPINE_PROJECTION')
            if len(taxa)<3:
                result.update(taxa=taxa,n_ingroup=len(set(species_tree.tips())-set(outgroups)))
                result['flags'].append('INSUFFICIENT_TERMINAL_SPINE')
                result['evidence']=grade(result)
                return result
    elif taxa is None:
        taxa=[n for n in processed[0].tips() if n not in outgroups]
    taxa=list(taxa)
    if len(taxa)<3 or len(set(taxa))!=len(taxa):raise ValueError('Supply at least three distinct ingroup taxa in backbone order')
    result.update(taxa=taxa,n_ingroup=len(taxa))
    # For a ladder, the familiar three-species projection drops the new recipient sister.
    three=[taxa[0],taxa[2],taxa[3]] if len(taxa)>=4 else taxa[:3]
    if nroot:
        tri=triplet(processed,three);labels=np.array(tri.pop('labels'));result['triplet']=tri
        if distances:
            ds=[];ls=[];ods=[];missing_lengths=0
            for t,label in zip(processed,labels):
                if label==-2:continue
                try:
                    d=patristic(t,three);ds.append([d[0,1],d[1,2],d[0,2]]);ls.append(label)
                    if outgroups:
                        od=patristic(t,three+[outgroups[0]]);ods.append(od[:3,3])
                except ValueError:missing_lengths+=1
            if missing_lengths:result['flags'].append('MISSING_OR_INVALID_BRANCH_LENGTHS')
            result['n_distance_loci']=len(ds)
            result['distance_missing_or_invalid_lengths']=missing_lengths
            if len(ds)>=2:
                ds=np.array(ds);ls=np.array(ls)
                result['six_moment']=six_moment_fit(ds)
                result['distance_calibration']=calibrate(ds,ls,reps,seed,np.array(ods) if len(ods)==len(ds) else None,
                                                         tri['resolved']/max(1,len(processed)))
        if len(taxa)>=4:result['four_leaf']=four_leaf(processed,taxa[:4],reps,seed,gtee_floor)
        if len(taxa)>4:
            grid=[]
            for reference in taxa[3:]:
                q=triplet(processed,[taxa[0],taxa[2],reference]);q.pop('labels')
                grid.append({'reference_C':reference,'open_wedge':q['open_wedge'],'a':q['a'],'b':q['b'],
                             'counts_AB_BC_AC':q['counts_AB_BC_AC'],'resolved':q['resolved']})
            result['cross_c_grid']={'rows':grid,'wedge_status_changes':len({r['open_wedge'] for r in grid})>1,
                                    'scope':'Fixed recipient and B, varying sampled outer C; projection sensitivity, not a ghost-history classification.'}
    if len(taxa)>=6 and (not species_tree or result.get('spine',{}).get('terminal_ladder')):
        # For n>6 the certified first five are held fixed; every outer representative
        # must pass its own full 70-candidate reduction. This is conservative.
        marginals=[six_leaf(processed,taxa[:5]+[outer],reps,seed) for outer in taxa[5:]]
        result['six_leaf']=dict(marginals[0]);result['six_leaf']['all_outer_marginals']=marginals
        result['six_leaf']['passed']=all(m['passed'] for m in marginals)
        result['six_leaf']['n_outer_marginals']=len(marginals)
    if not nroot and len(taxa)<6:result['flags'].append('ROOT_REQUIRED_FOR_TRIPLET_OR_FOUR_LEAF_TEST')
    if result.get('distance_calibration',{}).get('quasi_clock_checks',{}).get('resolved_fraction_healthy') is False:
        result['flags'].append('LOW_DIVERGENCE_STAR_TREE_WARNING')
    result['evidence']=grade(result)
    return result

def scan_trees(trees,species_tree,outgroups=(),**kwargs):
    tree=load_tree(species_tree) if isinstance(species_tree,(str,Path)) else species_tree
    trees=list(trees);rows=[]
    for recipient in tree.tips():
        if recipient in outgroups:continue
        r=audit_trees(trees,species_tree=tree,recipient=recipient,outgroups=outgroups,**kwargs)
        rows.append({'recipient':recipient,'tier':r['evidence']['tier'],'result':r})
    rows.sort(key=lambda r:(r['tier'],r['recipient']))
    return {'analysis':'exploratory_scan','rows':rows,
            'scope':'Exploratory ranking across recipients and chosen projections; no family-wise calibration. Confirm candidates on held-out loci before inference.'}
