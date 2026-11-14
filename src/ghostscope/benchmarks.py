"""Read-only project replay. Expected outcomes are checks, never analysis inputs.

Numerical baselines and finite-locus power outcomes are verified read-only
against frozen project artifacts without threshold tuning or hard-coded answers.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import types
import json
import subprocess
import sys
import time
import numpy as np
from .analysis import audit_trees,scan_trees
from .cli import sequence_audit
from .io_parsers import read_trees,load_tree,TaxonMapper,induced,patristic,canonical,binary
from .topology_engine import cross_c_table,triplet
from .distance_engine import calibrate,six_moment_fit,bootstrap_mean,sequence_multi_rows
from .mcmc_engine import audit_chains,audit_summary,audit_summaries,ti_audit
from .report import clean,export

PACKAGE=Path(__file__).resolve().parents[2]
R8='round8_456_boundaries_20261108'
R5='round5_empirical_scouting'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def write_tsv(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
        for r in rows:w.writerow({k:json.dumps(clean(v),ensure_ascii=False) if isinstance(v,(dict,list,tuple)) else v for k,v in r.items()})

class Replay:
    def __init__(self,root,out,reps):
        self.root=Path(root).resolve();self.out=Path(out);self.reps=reps;self.checks=[];self.results={};self.inputs={};self.controls=[]
        self.r8='08_Higher_Order_4to6_and_ArbitraryN_R8' if (self.root/'08_Higher_Order_4to6_and_ArbitraryN_R8').exists() else R8
        self.r5='06_Empirical_Data_Audits_R5' if (self.root/'06_Empirical_Data_Audits_R5').exists() else R5
    def source(self,relative):
        p=self.root/relative
        if not p.is_file():raise FileNotFoundError(str(p))
        self.inputs[str(p.relative_to(self.root))]=sha(p)
        return p
    def check(self,name,ok,observed=None,expected=None,kind='REPLAY'):
        self.checks.append({'check':name,'status':'PASS' if ok else 'FAIL','kind':kind,'observed':clean(observed),'expected':clean(expected)})
    def section(self,name,run):
        start=time.monotonic()
        try:run();print(name+': completed',flush=True)
        except Exception as error:
            self.check(name+': execution',False,type(error).__name__+': '+str(error),'Complete available input processing',kind='EXECUTION')
            print(name+': '+str(error),flush=True)
        self.results.setdefault('timings',{})[name]=time.monotonic()-start
    def control(self,folder,sequence=False,regime='four_ghost_500'):
        d=json.loads(self.source(folder+'/design.json').read_text());tree=load_tree(d['trees'][0]);labels=dict(d['labels'])
        out=[sp for sp in labels.values() if sp not in tree.tips()]
        kw=dict(taxa=tree.tips(),species_tree=tree,recipient=tree.tips()[0],outgroups=out,delimiter='_',reps=self.reps,distances=False)
        if sequence:
            args=argparse.Namespace(input=str(self.root/folder/'loci'),format='auto',species_tree=d['trees'][0],taxa=','.join(tree.tips()),recipient=tree.tips()[0],outgroup=','.join(out),imap=None,taxon_delim='_',bootstrap=self.reps,seed=20261111,gtee_floor=.015,rooted=False,haplotypes_a=None,taxon_b=None,taxon_c=None,nj_trees=True,survival_window=None)
            r=sequence_audit(args)
        else:r=audit_trees(list(read_trees(self.source(folder+'/genealogies.jsonl'))),**kw)
        actual=r['evidence']['tier'];four=r.get('four_leaf',{});six=r.get('six_leaf',{})
        if regime=='four_ghost_500':
            delta_sig=bool(four.get('ci95') and four['ci95'][0][0]>0 and four['p_one_sided'][0]<.01)
            wedge_routed='SCOPE_ROUTING: OPEN_TOPOLOGY_WEDGE' in r['evidence']['flags']
            ok=bool(delta_sig and wedge_routed and actual in ('A','D'))
            expected='Delta>0 (p<0.01), open-wedge routed, tier in {A,D} at 500 loci'
        elif regime=='six_ghost_500':
            wedge_routed='SCOPE_ROUTING: OPEN_TOPOLOGY_WEDGE' in r['evidence']['flags']
            power_bounded=bool(six and six.get('n_significant_excluded',70)<70 and actual=='D')
            ok=bool(wedge_routed and power_bounded)
            expected='Finite-500-locus power bound (<70/70 IUT rejections -> Tier D, not forged Tier A)'
        else:
            ok=actual!='A'
            expected='not A'
        row={'dataset':Path(folder).name,'input':'500bp_NJ' if sequence else 'latent_genealogies','loci':r.get('n_loci',r['n_trees']),
             'expected_regime':expected,'observed_tier':actual,'status':'PASS' if ok else 'FAIL',
             'four_estimates':four.get('estimate'),'four_tolerance_pass':four.get('tolerance_pass'),
             'six_excluded':six.get('n_significant_excluded'),'six_signs':six.get('n_sign_consistent'),
             'flags':r['evidence']['flags']}
        self.controls.append(row);self.results.setdefault('controls',[]).append(r)
        self.check('Simulation control '+row['dataset']+' '+row['input'],ok,actual,expected,kind='FINITE_LOCUS_BENCHMARK')
    def positive_controls(self):
        for i in range(1,6):self.control(f'{self.r8}/ingroup_4/4I-Exp0_r{i:03d}',regime='four_ghost_500')
        for i in range(1,4):self.control(f'{self.r8}/ingroup_4/4I-Exp0_r{i:03d}',True,regime='four_ghost_500')
        for exp in (1,2):
            for i in range(1,4):self.control(f'{self.r8}/ingroup_4/4I-Exp{exp}_r{i:03d}',regime='non_ghost')
        for i in range(1,4):self.control(f'{self.r8}/ingroup_6/6I-Exp0_r{i:03d}',regime='six_ghost_500')
        self.control(f'{self.r8}/ingroup_6/6I-Exp0_r001',True,regime='six_ghost_500')
        self.control(f'{self.r8}/ingroup_6/6I-Exp1-E_r010',regime='non_ghost')
        self.control(f'{self.r8}/ingroup_6/6I-Parallel_r010',regime='non_ghost')
        self.check('Parallel donor warning',any('TYPE_III_PARALLEL_DONOR_WARNING' in r['flags'] for r in self.controls if 'Parallel' in r['dataset']),True,True)
    def thuja(self):
        refroot=self.root/self.r5/'round5B/results'
        baseline=list(csv.DictReader(self.source(f'{self.r5}/round5B/results/FIGURE5_DECOMPOSITION_DATA.tsv').open(),delimiter='\t'))
        row=next(r for r in baseline if r['dataset_id']=='pang2024')
        old=json.loads(self.source(f'{self.r5}/round5B/results/THUJA_FOUR_TIP_TREE_DIAGNOSTICS.json').read_text());taxa=list(old['taxa'].values())
        path=self.source(f'{self.r5}/round5B/extracted/data/emprical-result/Thuja/PhyloNet/PhyloNet.input')
        raw=list(read_trees(path));distance_complete=[]
        for t in raw:
            try:patristic(t,taxa)
            except ValueError:continue
            distance_complete.append(t)
        # The historical combined distance/topology audit excluded branchless
        # author members. A topology-only CLI may legitimately count them.
        r=audit_trees(distance_complete,taxa=taxa,rooted=True,reps=self.reps,seed=20261106)
        r['author_input_members_without_complete_distance_paths']=len(raw)-len(distance_complete)
        self.results['thuja_author_trios']=r
        for key in ('Delta_raw','Delta_B','Delta_AC','Delta_B_clock','rA_over_rB','rC_over_rB'):
            a=r['distance_calibration'][key];b=float(row[key]);self.check('Thuja author-trio '+key,abs(a-b)<1e-12,a,b)
        self.check('Thuja author-trio wedge',r['triplet']['open_wedge'],r['triplet']['counts_AB_BC_AC'],[int(row[k]) for k in ('nAB','nBC','nAC')])
        seqdir=self.root/self.r5/'round5B/extracted/data/emprical-result/Thuja/sequence'
        # Directly reparse the real quartet trees, and keep this version separate.
        files=sorted(seqdir.glob('*.treefile.gz'));ts=[next(read_trees(p)) for p in files]
        q=audit_trees(ts,taxa=taxa,outgroups=[old['outgroup']],reps=self.reps,seed=20261106)
        self.results['thuja_four_tip_reparse']=q
        self.check('Thuja four-tip topology counts',q['triplet']['counts_AB_BC_AC']==old['counts'],q['triplet']['counts_AB_BC_AC'],old['counts'])
        contrees=[next(read_trees(p)) for p in sorted(seqdir.glob('*.contree.gz'))]
        consensus=audit_trees(contrees,taxa=taxa,outgroups=[old['outgroup']],reps=self.reps,seed=20261106)
        self.results['thuja_consensus_four_tip_reparse']=consensus
        self.check('Thuja uncollapsed treefiles vs consensus quartet kept distinct',q['triplet']['counts_AB_BC_AC']!=consensus['triplet']['counts_AB_BC_AC'],q['triplet']['counts_AB_BC_AC'],consensus['triplet']['counts_AB_BC_AC'])
        self.check('Thuja consensus quartet matches author-trio counts',consensus['triplet']['counts_AB_BC_AC']==r['triplet']['counts_AB_BC_AC'],consensus['triplet']['counts_AB_BC_AC'],r['triplet']['counts_AB_BC_AC'])
        self.check('Thuja six-moment positive construction',r['six_moment']['compatible'] and r['six_moment']['max_abs_residual']<1e-15,r['six_moment'],True)
    def brassicaceae(self):
        path=self.source(f'{self.r5}/datasets/brassicaceae2020/extracted/Brassicaceae_phylogenomics-master/R_test_dir/cat_trees_nuc170213.txt')
        old=json.loads(self.source(f'{self.r5}/round5B/results/brassicaceae2020.json').read_text())
        table=list(csv.DictReader(self.source(f'{self.r5}/round5B/results/FIGURE5_DECOMPOSITION_DATA.tsv').open(),delimiter='\t'))
        row=next(r for r in table if r['dataset_id']=='brassicaceae2020')
        # Read the role prefixes from the archived original analysis metadata.
        metadata=json.loads(self.source(f'{self.r5}/round5B/history/EMPIRICAL_THEOREM_DIAGNOSTICS.json').read_text())
        base=next(r for r in metadata if r['dataset_id']=='brassicaceae2020')
        taxa=[x.rstrip('_') for x in base['taxa'].values()];og=base['outgroup'].rstrip('_')
        ts=list(read_trees(path));r=audit_trees(ts,taxa=taxa,outgroups=[og],delimiter='_',reps=self.reps,seed=20261106)
        self.results['brassicaceae_primary']=r
        a=r['distance_calibration']['Delta_AC_fraction'];b=abs(float(row['Delta_AC_fraction']))
        self.check('Brassicaceae pruned path fraction',abs(a-b)<1e-12,a,b)
        self.check('Brassicaceae moment construction',r['six_moment']['compatible'],r['six_moment'],True)
        backbone=next(read_trees(self.source(f'{self.r5}/datasets/brassicaceae2020/extracted/Brassicaceae_phylogenomics-master/R_test_dir/phylonet171218.tre')))
        mapper=TaxonMapper(delimiter='_',known=backbone.tips());prepared=[]
        for t in ts:
            keep=[name for name in t.tips() if mapper(name) in backbone.tips()]
            prepared.append(mapper.tree(induced(t,keep)))
        scan=scan_trees(prepared,backbone,outgroups=[og],reps=self.reps,distances=False)
        self.results['brassicaceae_full_scan']=scan
        self.check('Brassicaceae six-core-species full scan',len(scan['rows'])==len(backbone.tips())-1,len(scan['rows']),len(backbone.tips())-1)
        for key in ('Delta_raw','Delta_B','Delta_AC'):
            self.check('Brassicaceae '+key,abs(r['distance_calibration'][key]-old['star_decomposition'][key])<1e-12,r['distance_calibration'][key],old['star_decomposition'][key])
    def jaltomata(self):
        refs=json.loads(self.source(f'{self.r5}/round5B/results/JALTOMATA_ALL_TRIPLETS.json').read_text());results=[];nfit=0;nwedge=0;old_lost=0;new_lost=0
        for ref in refs:
            # Relocate the recorded absolute source under the explicit project root.
            original=Path(ref['source']);rel=(original.parts[original.parts.index(R5)+1:] if R5 in original.parts else (original.parts[original.parts.index(self.r5)+1:] if self.r5 in original.parts else original.parts));path=self.source(str(Path(self.r5,*rel)))
            taxa=list(ref['taxa'].values());ts=list(read_trees(path));r=audit_trees(ts,taxa=taxa,rooted=True,reps=1000,seed=20261106)
            results.append(r);nwedge+=r['triplet']['open_wedge'];nfit+=r['six_moment']['compatible']
            q=r['distance_calibration'];oldq=ref['star_decomposition']
            for key in ('Delta_raw','Delta_B','Delta_AC','Delta_B_clock'):
                self.check(ref['dataset_id']+' '+key,abs(q[key]-oldq[key])<1e-12,q[key],oldq[key])
            a=oldq['bootstrap']['Delta_B_clock']['p_normal'];b=q['bootstrap']['Delta_B_clock']['p_two_sided']
            self.check(ref['dataset_id']+' bootstrap B-clock p',abs(a-b)<1e-10,b,a)
            old_lost+=a>=.05;new_lost+=b>=.05
        self.results['jaltomata']=results
        self.check('Jaltomata all source trios in wedge',nwedge==len(refs),nwedge,len(refs))
        self.check('Jaltomata all six-moment compatible',nfit==len(refs),nfit,len(refs))
        self.check('Jaltomata nominal significance loss',old_lost==new_lost,new_lost,old_lost)
    def honeycreeper(self):
        from .io_parsers import parse_newick
        p=self.source(f'{self.r5}/datasets/honeycreeper/extracted/ghostparser/ghostparser.txt')
        rp=self.source(f'{self.r5}/datasets/honeycreeper/extracted/SNaQ/SNaQ_Rep1/network_0_reticulations.out')
        candidates=[line for line in rp.read_text().splitlines() if line.startswith('(')]
        reference=parse_newick(min(candidates,key=lambda s:float(s.rsplit('=',1)[-1])).split(';')[0]+';')
        r=cross_c_table(p,reference.tips());self.results['honeycreeper']=r
        # Directive totals are external acceptance claims, deliberately confined
        # to this benchmark module, never fed into any estimator.
        self.check('Honeycreeper wedge rows',r['open_wedge_rows']==5433,r['open_wedge_rows'],5433)
        self.check('Honeycreeper all ghost-labelled rows in wedge',r['ghost_label_rows']==r['ghost_label_open_wedge_rows'],r['ghost_label_open_wedge_rows'],r['ghost_label_rows'])
        self.check('Honeycreeper ordered binary cross-C switch',r['ordered']['switching_rate']==1.,r['ordered'],1.)
        expected=list(csv.DictReader(self.source(f'{self.r8}/empirical_new_boundaries_20261108/results/E2_HONEYCREEPER_CONFLICT_DENOMINATORS.tsv').open(),delimiter='\t'))
        for row in expected:
            if row['stratum']!='all_mapped' or row['level']=='individual':continue
            got=r[row['level']][row['direction']]
            self.check('Honeycreeper '+row['level']+' '+row['direction']+' denominator',got['ghost_ever_pairs']==int(row['ghost_candidate_pairs']) and got['binary_label_switch_pairs']==int(row['cross_C_conflicting_pairs']),got,row)
    def panthera(self):
        refs=json.loads(self.source(f'{self.r5}/round5C_human/panthera/results/PANTHERA_DIAGNOSTICS.json').read_text());ref=next(r for r in refs if r['dataset_id'].endswith('combined'))
        roles=ref['primary_roles'];hap=roles[0]
        # Copy identities are in the archived comparison grid: recover the second
        # recipient sample from the same species prefix in the actual alignment.
        rows=[];loci=0
        from .io_parsers import read_alignments
        for part in ('odd','even'):
            path=self.source(f'{self.r5}/round5C_human/panthera/extracted/dataset/Panthera/data.{part}.nocat.txt')
            for _,seq in read_alignments(path):
                loci+=1
                # Primary A is lion, B leopard, C jaguar; source suffix labels do
                # not define theorem roles. The second lion shares suffix B2.
                suffix=hap.rsplit('^',1)[1][:-1]+'2'
                candidates=[n for n in seq if '^' in n and n.rsplit('^',1)[1]==suffix]
                if len(candidates)!=1:raise ValueError('Ambiguous second recipient in deposited Panthera mapping')
                row=sequence_multi_rows(seq,[hap,candidates[0]],roles[1],roles[2],roles[3])
                if row is not None:rows.append(row)
        r=bootstrap_mean(rows,1000,20261107);self.results['panthera']={'n_raw_loci':loci,'n_callable_double':len(rows),'squared_mismatch_screen':r,'scope':'Two deposited lion pseudohaploid individuals; finite-sequence proxy, not demonstrated exchangeable latent genealogical distances.'}
        expected=ref['T10_T12_sequence_screen']['I']
        self.check('Panthera double-copy point estimates',np.max(abs(np.array(r['estimate'])-expected['estimate']))<1e-16,r['estimate'],expected['estimate'])
        self.check('Panthera block-bootstrap p',abs(r['p_two_sided'][2]-expected['p_normal'][2])<1e-10,r['p_two_sided'][2],expected['p_normal'][2])
    def mcmc(self):
        initial=self.source(f'{self.r8}/empirical_new_boundaries_20261108/results/E4_DUALCHAIN_COMPLETE_AUDIT.json');final=self.source(f'{self.r8}/empirical_extension68_20261108/results/AMD_FINAL/05_results/EXTENSION68_ALL_RESULTS.json')
        merged=audit_summaries([initial,final]);tail=audit_summary(final);self.results['paired_mcmc']=merged
        self.check('Five bulk-stable rare-tail decision splits',tail['rare_tail_split_pairs']==5,tail['rare_tail_split_pairs'],5)
        rare=json.loads(self.source(f'{self.r8}/empirical_extension68_20261108/results/AMD_FINAL/05_results/RARE_REGION_BF_AUDIT.json').read_text())
        indexed={r['id']:r for r in merged['rows']}
        for r in rare['results']:
            row=indexed[r['id']];edge=next(e for e in row['phi'] if e['parameter']==r['parameter']);counts=[z['draws_below_eps'] for z in r['chains']]
            self.check('Long-tail counts '+r['id'],counts==edge['count_below_eps_by_chain'],edge['count_below_eps_by_chain'],counts)
        supported=[p for p in merged['paired'] if 'SINGLE_GRAPH_BF_ILLUSION_DETECTED' in p['flags']]
        self.check('Paired-network illusion actually detected',bool(supported),len(supported),'>0')
        # All convergence flags stay stage-specific; saved summaries are not raw replay.
        self.results['rate_model_shifts']=[]
        groups={}
        for r in merged['rows']:groups.setdefault((r.get('system'),r.get('loci'),r.get('model'),r.get('partition'),r.get('samples_config','single')), {})[r.get('mode')]=r
        for key,modes in groups.items():
            if 'default' not in modes or 'locusrate' not in modes:continue
            a=modes['default'].get('phi',[]);b=modes['locusrate'].get('phi',[])
            if len(a)==len(b)==1:
                self.results['rate_model_shifts'].append({'group':key,'locusrate_minus_default_phi':b[0]['posterior_mean_phi']-a[0]['posterior_mean_phi']})
        # Truly independent recomputation from local original draws.
        reference_path=self.source(f'{self.r8}/tools/dbpp_convergence.py')
        module=types.ModuleType('r8_readonly_numeric_reference');module.__file__=str(reference_path)
        # Execute definitions from text: importlib loaders could create a .pyc in
        # the protected source directory when invoked without PYTHONDONTWRITEBYTECODE.
        exec(compile(reference_path.read_text(),str(reference_path),'exec'),module.__dict__)
        pairs=[]
        cap=self.root/self.r8/'results/capability'
        for p in sorted(cap.glob('*chain1/chain.mcmc.txt')):
            q=p.parent.with_name(p.parent.name.replace('chain1','chain2'))/p.name
            if q.exists():pairs.append((p,q))
        deltas=[]
        short_windows=0;numeric_windows=0
        def compare(got,ref):
            nonlocal short_windows,numeric_windows
            if 'max_rhat' not in ref:
                short_windows+=1
                return got.get('reason')==ref.get('reason') and got['accepted']==ref['accepted']
            numeric_windows+=1
            a=[got['max_rhat'],got['min_ess_bulk']];b=[ref['max_rhat'],ref['min_ess_bulk']]
            for x,y in zip(a,b):
                deltas.append(0. if x==y else abs(x-y))
            return all(np.isfinite(d) and d<1e-10 for d in deltas)
        matches=[]
        for paths in pairs:
            got=audit_chains(paths)
            for name,refpaths in [('cumulative',paths)]:
                ref=module.diagnostics(refpaths);matches.append(compare(got[name],ref))
            traces=[module.read_trace(p) for p in paths];n=len(traces[0][1])
            # Reference reader is replaced in memory only to evaluate the recent
            # window with the very same published function, without temp traces.
            oldread=module.read_trace;mapping={str(p):(h,a[n//2:]) for p,(h,a) in zip(paths,traces)}
            module.read_trace=lambda p:mapping[str(p)]
            try:ref=module.diagnostics(paths)
            finally:module.read_trace=oldread
            matches.append(compare(got['recent_half'],ref))
        tiroot=self.root/self.r8/'results/TI_official_smoke_v2';ti=ti_audit([tiroot/'chain1',tiroot/'chain2']);self.results['TI']=ti
        for r in ti['nodes']:
            paths=r['chains'];ref=module.diagnostics(paths);matches.append(compare(r['cumulative'],ref))
        self.check('Local raw R8 numerical replay and short-window gates',bool(matches) and all(matches),{'capability_pairs':len(pairs),'numeric_windows':numeric_windows,'too_short_windows':short_windows,'max_finite_difference':max(deltas,default=0)},'<1e-10 where estimable; identical refusal for too few samples')
        self.check('All local TI quadrature nodes audited',len(ti['nodes'])==16,len(ti['nodes']),16)
        self.check('TI failure is retained',not ti['all_nodes_strict'],ti['all_nodes_strict'],False)

def reproduce(project_root,outdir=None,reps=2000,protection=True):
    project_package=Path(project_root).resolve()/('09_GhostScope_Toolkit_v1' if (Path(project_root).resolve()/'09_GhostScope_Toolkit_v1').exists() else 'ghostscope')
    out=Path(outdir) if outdir else project_package/'validation_outputs';out.mkdir(parents=True,exist_ok=True)
    replay=Replay(project_root,out,reps);start=time.monotonic()
    for name,method in [('simulation_controls',replay.positive_controls),('thuja',replay.thuja),('brassicaceae',replay.brassicaceae),('jaltomata',replay.jaltomata),('honeycreeper',replay.honeycreeper),('panthera',replay.panthera),('mcmc',replay.mcmc)]:replay.section(name,method)
    protection_record={'status':'NOT_RUN'}
    if protection:
        script=project_package/'tools/protect.py'
        is_public_archive=(Path(project_root).resolve()/'PUBLIC_ARCHIVE_SHA256.tsv').is_file() and not (Path(project_root).resolve()/'round8_456_boundaries_20261108/PROTECTED_BASELINE.json').is_file()
        if script.is_file() and (is_public_archive or (project_package/'validation_outputs/PREDEVELOPMENT_HASHES.json').is_file()):
            result=subprocess.run([sys.executable,str(script),'final'],capture_output=True,text=True,env={**__import__('os').environ,'PYTHONDONTWRITEBYTECODE':'1'})
            live_path=project_package/'validation_outputs'/('PROTECTION_PUBLIC_ARCHIVE_LIVE.json' if is_public_archive else 'PROTECTION_FINAL.json')
            if live_path.is_file():
                try:protection_record=json.loads(live_path.read_text())
                except Exception:protection_record={'status':'PARSE_ERROR','returncode':result.returncode}
            if result.returncode==0 and isinstance(protection_record,dict) and protection_record.get('baseline_entries',0)>0:
                if is_public_archive:
                    replay.check('Public archive manifest frozen inputs zero changes',not protection_record['baseline_changed_or_removed'],len(protection_record['baseline_changed_or_removed']),0,kind='PUBLIC_ARCHIVE_INTEGRITY')
                    replay.check('All replay-referenced archive inputs unchanged',not protection_record['session_changed_or_removed'],len(protection_record['session_changed_or_removed']),0,kind='PUBLIC_ARCHIVE_INTEGRITY')
                    replay.check('No unauthorized files added to frozen archive',not protection_record['session_added_outside_package'],len(protection_record['session_added_outside_package']),0,kind='PUBLIC_ARCHIVE_INTEGRITY')
                else:
                    replay.check('Historical baseline zero changes',not protection_record['baseline_changed_or_removed'],len(protection_record['baseline_changed_or_removed']),0,kind='PROTECTION')
                    replay.check('All predevelopment files unchanged',not protection_record['session_changed_or_removed'],len(protection_record['session_changed_or_removed']),0,kind='PROTECTION')
                    replay.check('No files added outside new package',not protection_record['session_added_outside_package'],len(protection_record['session_added_outside_package']),0,kind='PROTECTION')
            else:
                err_obs=protection_record if isinstance(protection_record,dict) and protection_record.get('status')!='NOT_RUN' else (result.stderr or result.stdout)
                replay.check('Public archive live hash verification' if is_public_archive else 'Protection verifier execution',False,err_obs,'Zero modified/missing frozen files and zero unauthorized additions',kind='PUBLIC_ARCHIVE_INTEGRITY' if is_public_archive else 'PROTECTION')
        else:replay.check('Archive or predevelopment protection manifest available',False,'Not present','PUBLIC_ARCHIVE_SHA256.tsv or PREDEVELOPMENT_HASHES.json',kind='PROTECTION')
    failures=[r for r in replay.checks if r['status']!='PASS'];execution=[r for r in failures if r['kind']=='EXECUTION']
    status='PASS' if not failures and protection else ('REPLAY_EXECUTION_ERRORS' if execution else 'BENCHMARK_MISMATCH')
    certificate={'status':status,'checks':replay.checks,'passed':len(replay.checks)-len(failures),'total':len(replay.checks),
                 'protection':protection_record,'input_sha256':replay.inputs,'timings':replay.results.get('timings'),
                 'seconds':time.monotonic()-start,'scope':'No historical input or manuscript is modified. Finite-500-locus power limits are reported transparently without tuning thresholds. 500k summaries are explicitly separate from raw-draw replay.'}
    (out/'VALIDATION_CERTIFICATE.json').write_text(json.dumps(clean(certificate),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    (out/'PAPER_EMPIRICAL_REPLAY.json').write_text(json.dumps(clean(replay.results),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    write_tsv(out/'PAPER_EMPIRICAL_SUMMARY.tsv',replay.checks);write_tsv(out/'POSITIVE_CONTROLS_REPLAY.tsv',replay.controls)
    return {'analysis':'paper_replay','status':status,'passed':certificate['passed'],'total':certificate['total'],
            'details':str(out/'VALIDATION_CERTIFICATE.json'),'scope':certificate['scope']}
