"""Command-line entry point. No files are written unless --outdir is given."""
import argparse
from importlib.resources import files
from pathlib import Path
import json
import sys
import time
import numpy as np
from . import __version__
from .analysis import audit_trees,scan_trees
from .io_parsers import read_trees,read_alignments,distance_matrix,neighbor_joining,TaxonMapper,load_tree
from .distance_engine import bootstrap_mean,sequence_multi_rows,survival_curvature,site_patterns,summarize_patterns,multi_haplotype_distances
from .grading import grade
from .mcmc_engine import audit_chains,paired_audit,audit_summary,audit_summaries,ti_audit
from .bpp_generator import make_bpp
from .report import card,export
from .topology_engine import cross_c_table

def csvlist(value):return [s.strip() for s in value.split(',') if s.strip()] if value else []

def common(p):
    p.add_argument('-o','--outdir',help='Write JSON, TSV and Markdown here; default is terminal only')

def tree_options(p):
    p.add_argument('-s','--species-tree',help='Rooted species tree: Newick text or a file')
    p.add_argument('--taxa',help='Comma-separated ordered ingroup tips along the ladder (recipient first)')
    p.add_argument('--recipient',help='Candidate recipient species')
    p.add_argument('--outgroup',help='Single outgroup or comma-separated monophyletic outgroup clade')
    p.add_argument('--rooted',action='store_true',help='Declare that input gene trees are already rooted')
    p.add_argument('--imap',help='Two columns: sample species; BPP suffix identifiers supported')
    p.add_argument('--taxon-delim',help='Split sample labels at this delimiter; exact species-tree labels are preserved')
    p.add_argument('--bootstrap',type=int,default=2000,help='Locus bootstrap replicates (default: 2000; minimum: 100)')
    p.add_argument('--seed',type=int,default=20261111)
    p.add_argument('--gtee-floor',type=float,default=.015,help='User tolerance for rooted-4 contrasts; not a calibrated error bound')
    common(p)

def parser():
    p=argparse.ArgumentParser(prog='ghostscope',description='Audit ghost-introgression evidence within explicit candidate models.',epilog='Tiers are conditional evidence routes, not biological truth labels. Try: ghostscope demo')
    p.add_argument('--version',action='version',version=__version__)
    subs=p.add_subparsers(dest='command',required=True)
    t=subs.add_parser('trees',help='Audit gene trees, with optional distance decomposition',epilog='Example: ghostscope trees -g genes.nwk --taxa A,D,B,C --outgroup O')
    t.add_argument('-g','--gene-trees',help='Newick/.gz, genealogies.jsonl, or directory (required unless --cross-c-table is given)')
    t.add_argument('--cross-c-table',help='Optional author-oriented ghostparser.tsv')
    t.add_argument('--haplotypes-a',help='Recipient sample tips; use with raw unmapped tree tips')
    t.add_argument('--taxon-b');t.add_argument('--taxon-c')
    t.add_argument('--assume-exchangeable',action='store_true',help='Assert same-species copies and common per-locus rate')
    t.add_argument('--latent-genealogies',action='store_true',help='Declare that distances are latent genealogical, not estimated sequence distances')
    tree_options(t)
    s=subs.add_parser('scan',help='Explore all candidate terminal recipients',epilog='Example: ghostscope scan -g genes.nwk -s species.nwk --outgroup O')
    s.add_argument('-g','--gene-trees',required=True);tree_options(s)
    q=subs.add_parser('seqs',help='Audit alignments; optional built-in JC69 + NJ',epilog='Example: ghostscope seqs -i loci/ --taxa A,D,B,C --outgroup O --nj-trees')
    q.add_argument('-i','--input',required=True)
    q.add_argument('--format',choices=['auto','fasta','phylip','bpp'],default='auto')
    q.add_argument('--nj-trees',action='store_true',help='Reconstruct each finite nonsaturated alignment using NJ')
    q.add_argument('--haplotypes-a',help='Comma-separated recipient sample names (>=2)')
    q.add_argument('--taxon-b');q.add_argument('--taxon-c')
    q.add_argument('--survival-window',help='Optional distance window LO:HI; exploratory, no automatic tier')
    tree_options(q)
    m=subs.add_parser('mcmc',help='Audit dual windows, rare tails and competing networks',epilog='Example: ghostscope mcmc --ghost g1.txt g2.txt --inflow i1.txt i2.txt')
    for flag in ('chains','ghost','inflow','outflow'):m.add_argument('--'+flag,nargs='+',metavar='TRACE')
    m.add_argument('--summary-json',nargs='+',help='Initial then continuation JSON files; later task records replace earlier ones');m.add_argument('--summary-tsv');m.add_argument('--ti-dirs',nargs=2)
    common(m)
    b=subs.add_parser('make-bpp',help='Generate 18 paired-network BPP v4.8 controls',epilog='Example: ghostscope make-bpp -s "((A,B),C);" --recipient A --sister B --out-clade C --seqfile data.bpp --imap imap.txt --nloci 100 -o bpp_suite')
    b.add_argument('-s','--species-tree',required=True);b.add_argument('--recipient',required=True);b.add_argument('--sister',required=True);b.add_argument('--out-clade',required=True)
    b.add_argument('--seqfile',required=True);b.add_argument('--imap',required=True);b.add_argument('--nloci',type=int,required=True)
    b.add_argument('-o','--outdir',required=True);b.add_argument('--seed',type=int,default=20261111)
    b.add_argument('--burnin',type=int,default=10000);b.add_argument('--nsample',type=int,default=20000);b.add_argument('--sampfreq',type=int,default=2)
    d=subs.add_parser('demo',help='Run compact positive/ambiguous input examples');common(d)
    au=subs.add_parser('audit',help='Compatibility alias for `demo` (--demo) or `trees`');au.add_argument('--demo',action='store_true');au.add_argument('--out',dest='outdir_alias');tree_options(au);au.add_argument('-g','--gene-trees')
    am=subs.add_parser('audit-mcmc',help='Compatibility alias for `mcmc`');
    for flag in ('chains','ghost','inflow','outflow'):am.add_argument('--'+flag,nargs='+',metavar='TRACE')
    am.add_argument('--summary-json',nargs='+');am.add_argument('--summary-tsv');am.add_argument('--ti-dirs',nargs=2);common(am)
    r=subs.add_parser('reproduce-paper',help='Read-only replay of available project benchmarks')
    r.add_argument('--project-root',default=str(Path.cwd().parent if Path.cwd().name in ('ghostscope','09_GhostScope_Toolkit_v1') else Path.cwd()))
    r.add_argument('--bootstrap',type=int,default=2000);r.add_argument('--skip-protection',action='store_true',help='Skip lengthy historical hashing; certificate remains incomplete')
    common(r)
    return p

def _kwargs(args):
    return dict(taxa=csvlist(args.taxa) or None,species_tree=args.species_tree,recipient=args.recipient,
                outgroups=csvlist(args.outgroup),imap=args.imap,delimiter=args.taxon_delim,reps=args.bootstrap,
                seed=args.seed,gtee_floor=args.gtee_floor,rooted=args.rooted)

def sequence_audit(args):
    tree=load_tree(args.species_tree) if args.species_tree else None
    mapper=TaxonMapper(args.imap,args.taxon_delim,tree.tips() if tree else csvlist(args.taxa)+csvlist(args.outgroup))
    trees=[];n=0;nonfinite=0;zero=0;multi=[];within=[];pattern_rows=[];role_rows=[[],[],[]];locus_lengths=[];missing_multi=0
    taxa=csvlist(args.taxa);outgroups=csvlist(args.outgroup);haps=csvlist(args.haplotypes_a)
    if not taxa and tree:
        from .topology_engine import spine
        recipient=args.recipient or next(n for n in tree.tips() if n not in outgroups)
        info=spine(tree,recipient,outgroups)
        taxa=[recipient]+[c[0] for c in info['sibling_clades']]
    for tag,seq in read_alignments(args.input,args.format):
        n+=1;locus_lengths.append(len(next(iter(seq.values()))))
        if haps:
            if not args.taxon_b or not args.taxon_c:raise ValueError('--haplotypes-a requires --taxon-b and --taxon-c')
            row=sequence_multi_rows(seq,haps,args.taxon_b,args.taxon_c,outgroups[0] if outgroups else None)
            if row is not None:multi.append(row)
            else:missing_multi+=1
        # Multiple individuals are retained for sequence statistics; NJ topology
        # auditing needs explicitly distinct species, so the first mapped sample
        # is NOT silently chosen as a representative.
        mapped={};collision=False
        for name,s in seq.items():
            label=mapper(name)
            if label in mapped:collision=True
            mapped[label]=s
        if not taxa and not tree:taxa=[x for x in mapped if x not in outgroups]
        three=([taxa[0],taxa[2],taxa[3]] if len(taxa)>=4 else taxa[:3]) if taxa else []
        if len(three)==3 and outgroups and not collision:
            roles=[three,[three[0],three[2],three[1]],[three[1],three[0],three[2]]]
            for i,role in enumerate(roles):
                pat=site_patterns(mapped,role+[outgroups[0]])
                if pat is not None:role_rows[i].append(pat)
        names,d,counts=distance_matrix(seq)
        if not np.isfinite(d).all():nonfinite+=1;continue
        if haps and set(haps)<=set(names):
            within.append(d[names.index(haps[0]),names.index(haps[1])])
        if np.max(d)<=1e-12:zero+=1
        if args.nj_trees:
            if collision:raise ValueError('NJ topology audit has multiple samples per mapped species; supply single representatives or omit --nj-trees')
            reconstructed=neighbor_joining(names,d);trees.append(reconstructed)
    if not n:raise ValueError('No alignment loci found')
    if args.nj_trees and trees:
        kw=_kwargs(args)
        if not tree and not kw['taxa']:kw['taxa']=taxa
        result=audit_trees(trees,**kw)
    else:result={'n_trees':0,'n_ingroup':len(taxa),'flags':[]}
    result.update(analysis='sequence_audit',n_loci=n,unconditional_jc_defined=nonfinite==0,
                  nonfinite_or_saturated_loci=nonfinite,finite_conditioning_loci=n-nonfinite,zero_variation_loci=zero,
                  finite_sequence_scope='JC69 is a nonlinear finite-site estimator. Saturated/missing-pair loci are reported and excluded only from the conditional NJ stream, not called negative evidence.',
                  alignment_sites=sum(locus_lengths),length_range=[min(locus_lengths),max(locus_lengths)])
    if role_rows[0]:result['site_patterns']=summarize_patterns(role_rows[0])
    result['three_role_pattern_screens']=[summarize_patterns(rows) for rows in role_rows if rows]
    if multi:
        stat=bootstrap_mean(multi,args.bootstrap,args.seed);negative=stat['ci95'][2][1]<0
        stat.update(code='MULTI_HAPLOTYPE_SUFFICIENT_ZONE_CANDIDATE' if negative else 'MULTI_HAPLOTYPE_IN_CONFOUNDED_CONE',
                    status='FINITE_SEQUENCE_PROXY_ONLY' if negative else 'SCOPE_ROUTING',tier_b_eligible=False,n=len(multi),missing_loci=missing_multi,
                    columns=['IB_plugin','IC_plugin','IB_minus_IC','IB_iid_site_debiased','IC_iid_site_debiased','difference_iid_site_debiased'],
                    scope='Squared mismatch-probability differences; exchangeability, shared locus rates and latent genealogical distance assumptions are not established by this input.')
        result['multi_haplotype']=stat
    if within:result['survival']=survival_curvature(within,[float(v) for v in args.survival_window.split(':')] if args.survival_window else None,args.bootstrap,args.seed)
    result['evidence']=grade(result)
    return result

def demo():
    bundle=files('ghostscope').joinpath('data/demo_bundle')
    # Multiplicities are compact input data, not precomputed audit results.
    from .topology_engine import six_leaf_counts,_freeze
    from collections import Counter
    data=json.loads(bundle.joinpath('positive_topologies.json').read_text())
    counts=Counter({_freeze(r['coordinate']):r['count'] for r in data['topology_counts']})
    positive={'n_trees':data['n_loci'],'n_ingroup':6,'six_leaf':six_leaf_counts(counts,data['taxa']),
              'fixture_scope':data['scope'],'flags':[]}
    positive['evidence']=grade(positive)
    negative=audit_trees(list(read_trees(bundle.joinpath('ambiguous.nwk'))),taxa=['A','B','C'],outgroups=['O'],reps=2000)
    from .mcmc_engine import audit_summary
    paired=audit_summary(bundle.joinpath('paired_summary.json'))
    return {'analysis':'demo','cases':{'Positive topology control':positive,'Ambiguous depth control':negative,'Paired-network trap':paired},'scope':'Illustrative microdata; benchmarks and held-out controls are separate.'}

def main(argv=None):
    p=parser();args=p.parse_args(argv);start=time.monotonic()
    try:
        if args.command=='trees':
            if not args.gene_trees and not args.cross_c_table:raise ValueError('trees requires --gene-trees or --cross-c-table')
            trees=list(read_trees(args.gene_trees)) if args.gene_trees else []
            if not trees and args.cross_c_table:
                ref=load_tree(args.species_tree).tips() if args.species_tree else []
                result={'cross_c':cross_c_table(args.cross_c_table,ref),'flags':['CROSS_C_LABEL_SENSITIVITY_SCREEN']}
                result['evidence']=grade(result)
                print(card(result))
                if args.outdir:export(result,args.outdir)
                return 0
            if args.haplotypes_a:
                if not args.taxon_b or not args.taxon_c:raise ValueError('Specify --taxon-b and --taxon-c')
                from .io_parsers import patristic
                haps=csvlist(args.haplotypes_a);names=haps+[args.taxon_b,args.taxon_c]
                matrices=[patristic(t,names) for t in trees if set(names)<=set(t.tips())]
                m=multi_haplotype_distances(matrices,names,haps,args.taxon_b,args.taxon_c,args.bootstrap,args.seed,args.latent_genealogies,args.assume_exchangeable)
                result={'n_trees':len(trees),'multi_haplotype':m,'flags':[]};result['evidence']=grade(result)
            else:result=audit_trees(trees,**_kwargs(args))
            if args.cross_c_table:result['cross_c']=cross_c_table(args.cross_c_table,load_tree(args.species_tree).tips() if args.species_tree else [])
        elif args.command=='scan':
            if not args.species_tree:raise ValueError('scan requires --species-tree')
            kw=_kwargs(args);kw.pop('taxa');kw.pop('recipient');st=kw.pop('species_tree')
            result=scan_trees(read_trees(args.gene_trees),st,**kw)
        elif args.command=='seqs':result=sequence_audit(args)
        elif args.command in ('mcmc','audit-mcmc'):
            alternatives=sum(bool(getattr(args,k)) for k in ('chains','summary_json','summary_tsv','ti_dirs'))+int(any(getattr(args,k) for k in ('ghost','inflow','outflow')))
            if alternatives!=1:raise ValueError('Choose one input mode: --chains, competing models, summary, or --ti-dirs')
            if args.summary_json or args.summary_tsv:result=audit_summaries(args.summary_json) if args.summary_json else audit_summary(args.summary_tsv)
            elif args.ti_dirs:result=ti_audit(args.ti_dirs)
            elif args.chains:result=audit_chains(args.chains)
            else:result=paired_audit({k:audit_chains(getattr(args,k)) for k in ('ghost','inflow','outflow') if getattr(args,k)})
        elif args.command=='make-bpp':
            result=make_bpp(args.species_tree,args.recipient,args.sister,args.out_clade,args.seqfile,args.imap,args.nloci,args.outdir,args.seed,args.burnin,args.nsample,args.sampfreq)
        elif args.command=='demo' or (args.command=='audit' and getattr(args,'demo',False)):
            if getattr(args,'outdir_alias',None) and not getattr(args,'outdir',None):args.outdir=args.outdir_alias
            result=demo()
        elif args.command=='audit':
            if getattr(args,'outdir_alias',None) and not getattr(args,'outdir',None):args.outdir=args.outdir_alias
            if not getattr(args,'gene_trees',None):result=demo()
            else:result=audit_trees(list(read_trees(args.gene_trees)),**_kwargs(args))
        else:
            from .benchmarks import reproduce
            result=reproduce(args.project_root,args.outdir,args.bootstrap,not args.skip_protection)
        result['elapsed_seconds']=time.monotonic()-start
        print(card(result))
        if args.outdir and args.command not in ('make-bpp','reproduce-paper'):export(result,args.outdir)
        return 1 if args.command=='reproduce-paper' and result['status']!='PASS' else 0
    except (ValueError,OSError,KeyError,IndexError) as error:
        p.exit(2,'ghostscope: '+str(error)+'\n')

if __name__=='__main__':raise SystemExit(main())
