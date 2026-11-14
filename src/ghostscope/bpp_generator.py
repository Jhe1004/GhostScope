"""Construct paired MSci graphs; writes controls but never starts an analysis."""
from pathlib import Path
import hashlib
import json
import re
import shlex
from .io_parsers import Node, clone, load_tree, newick, parse_newick, read_alignments, TaxonMapper

MODES={'default':('1','0'),'locusrate':('1','1 0 0 5 iid'),'clock2':('2 2 2 2 iid LN','0')}

def msci_network(tree,recipient,sister,out_clade,model):
    tree=clone(tree);tips=tree.tips()
    def is_binary(n):return not n.children or len(n.children)==2 and all(is_binary(c) for c in n.children)
    if not is_binary(tree):raise ValueError('BPP generation requires a fully bifurcating rooted species backbone')
    if not {recipient,sister,out_clade}<=set(tips) or len({recipient,sister,out_clade})!=3:raise ValueError('Specify three different species present in the backbone')
    if len(tips)<3 or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',t) for t in tips):raise ValueError('BPP generation requires >=3 simple unique species labels (letters/digits/underscore)')
    # These anchor roles define ((recipient,sister),out_clade), not any arbitrary trio.
    from .io_parsers import induced,canonical
    from .topology_engine import node
    if canonical(induced(tree,{recipient,sister,out_clade}))!=node(node(recipient,sister),out_clade):
        raise ValueError('Anchor roles must induce ((recipient,sister),out-clade) in the rooted backbone')
    used=set(tips);counter=[0]
    def fresh(pre):
        n=pre
        while n in used:counter[0]+=1;n=pre+str(counter[0])
        used.add(n);return n
    h=fresh('H')
    def ann(parent):return [f'[&phi=0.500000,tau-parent={parent}]']
    target=recipient if model=='ghost' else (sister if model=='inflow' else out_clade)
    donor=out_clade if model=='inflow' else sister
    if model not in ('ghost','inflow','outflow'):raise ValueError('Unknown MSci model')
    def rec(n):
        if not n.children:
            if model=='ghost' and n.name==recipient:return Node(h,annotations=ann('yes'))
            if model!='ghost' and n.name==target:return Node(h,children=[Node(n.name)],annotations=ann('yes'))
            if model!='ghost' and n.name==donor:return Node(children=[Node(h,annotations=ann('no')),Node(n.name)])
            return Node(n.name)
        return Node(children=[rec(c) for c in n.children])
    graph=rec(tree)
    if model=='ghost':graph=Node(children=[Node(h,children=[Node(recipient)],annotations=ann('yes')),graph])
    def label(n,isroot=False):
        for c in n.children:label(c)
        if n.children and n.name!=h:n.name=fresh('R' if isroot else 'S')
    label(graph,True);text=newick(graph)
    validate_msci(text,tips,model)
    return text

def validate_msci(text,species,model):
    # Generic Newick rejects duplicate tips; the MSci reference is a hybrid tip
    # sharing a label with an internal hybrid node, so no duplicate leaves occur.
    tree=parse_newick(text);occurrences={};observed=[]
    def walk(n):
        if n.annotations:occurrences.setdefault(n.name,[]).append(n)
        elif not n.children:observed.append(n.name)
        for c in n.children:walk(c)
    walk(tree)
    if sorted(observed)!=sorted(species) or len(occurrences)!=1:raise ValueError('MSci graph has incorrect species/hybrid occurrences')
    nodes=next(iter(occurrences.values()))
    if len(nodes)!=2 or sorted(len(n.children) for n in nodes)!=[0,1]:raise ValueError('MSci hybrid must have exactly two occurrences')
    parents=sorted('yes' if 'tau-parent=yes' in n.annotations[0] else 'no' for n in nodes)
    if parents!=(['yes','yes'] if model=='ghost' else ['no','yes']):raise ValueError('Incorrect tau-parent semantics')
    return True

def make_bpp(species_tree,recipient,sister,out_clade,seqfile,imap,nloci,outdir,seed=20261111,burnin=10000,nsample=20000,sampfreq=2):
    tree=load_tree(species_tree);species=tree.tips();seqfile=Path(seqfile).resolve();imap=Path(imap).resolve()
    if not seqfile.is_file() or not imap.is_file():raise ValueError('seqfile and imap must exist')
    if any(ch.isspace() for ch in str(seqfile)+str(imap)):raise ValueError('BPP v4.8 suite input paths must not contain whitespace; place or link inputs at whitespace-free paths')
    if min(nloci,nsample,sampfreq)<1 or burnin<0 or seed<=0:raise ValueError('Invalid BPP run lengths or seed')
    if nloci<2:raise ValueError('The three-rate-model suite requires at least two loci for BPP locusrate')
    mapper=TaxonMapper(imap);counts={s:0 for s in species};nblocks=0;normalize=False
    for _,seq in read_alignments(seqfile):
        local={s:0 for s in species}
        for name in seq:
            sp=mapper(name)
            if sp not in local:raise ValueError('Alignment sample not mapped to a backbone species: '+name)
            local[sp]+=1
            if '^' not in name or mapper.mapping.get(name.rsplit('^',1)[1])!=sp:normalize=True
        for sp in species:counts[sp]=max(counts[sp],local[sp])
        nblocks+=1
    if nblocks!=nloci:raise ValueError(f'nloci={nloci} differs from parsed alignment blocks={nblocks}')
    if any(v==0 for v in counts.values()):raise ValueError('Some backbone species have no mapped samples')
    outdir=Path(outdir).resolve()
    if outdir.exists() and any(outdir.iterdir()):raise ValueError('Refusing to overwrite a nonempty BPP suite directory')
    outdir.mkdir(parents=True,exist_ok=True);tasks=[];commands=['#!/bin/sh','set -eu','cd "$(dirname "$0")"']
    ctl_seq=str(seqfile);ctl_imap=str(imap);normalization={}
    if normalize:
        # BPP v4.8 requires ^individual tags even if ordinary FASTA labels have
        # an imap entry. Preserve source data and normalize a copy in the suite.
        with (outdir/'normalized.bpp').open('w') as stream:
            for _,seq in read_alignments(seqfile):
                stream.write(f'{len(seq)} {len(next(iter(seq.values())))}\n')
                for name,sequence in seq.items():
                    if name not in normalization:normalization[name]=('I'+str(len(normalization)+1),mapper(name))
                    ident,sp=normalization[name];stream.write(sp+'^'+ident+' '+sequence+'\n')
        (outdir/'normalized.imap').write_text(''.join(ident+' '+sp+'\n' for ident,sp in normalization.values()))
        (outdir/'NORMALIZED_INPUT_MAP.tsv').write_text('original\tBPP_label\tspecies\n'+''.join(name+'\t'+sp+'^'+ident+'\t'+sp+'\n' for name,(ident,sp) in normalization.items()))
        ctl_seq='../../../normalized.bpp';ctl_imap='../../../normalized.imap'
    networks={m:msci_network(tree,recipient,sister,out_clade,m) for m in ('ghost','inflow','outflow')}
    for model,network in networks.items():
        for mode,(clock,locusrate) in MODES.items():
            for chain in (1,2):
                relative=Path(model)/mode/('chain'+str(chain));folder=outdir/relative;folder.mkdir(parents=True)
                taskseed=(seed+len(tasks)*104729)%2147483646+1
                text=f'''seed = {taskseed}
seqfile = {ctl_seq}
Imapfile = {ctl_imap}
jobname = chain
speciesdelimitation = 0
speciestree = 0
species&tree = {len(species)} {' '.join(species)}
  {' '.join(str(counts[s]) for s in species)}
  {network}
usedata = 1
nloci = {nloci}
cleandata = 1
thetaprior = 3 0.06 e
tauprior = 3 0.008
phiprior = 1 1
clock = {clock}
locusrate = {locusrate}
finetune = 1
print = 1 0 0 0
burnin = {burnin}
sampfreq = {sampfreq}
nsample = {nsample}
threads = 1 1 1
'''
                (folder/'run.ctl').write_text(text)
                tasks.append({'model':model,'mode':mode,'chain':chain,'seed':taskseed,'control':str(relative/'run.ctl'),'network':network})
                commands.append('(cd '+shlex.quote(str(relative))+' && "${BPP:-bpp}" --cfile run.ctl)')
    for mode in MODES:
        args=' '.join('--'+m+' '+ ' '.join(shlex.quote(str(Path(m)/mode/('chain'+str(c))/'chain.mcmc.txt')) for c in (1,2)) for m in networks)
        commands.append('ghostscope mcmc '+args+' -o '+shlex.quote('audit_'+mode))
    script=outdir/'run_and_audit.sh';script.write_text('\n'.join(commands)+'\n');script.chmod(0o755)
    manifest={'controls':tasks,'networks':networks,'nloci':nloci,'sample_counts':counts,'inputs_normalized_for_BPP':normalize,
              'seqfile_sha256':hashlib.sha256(seqfile.read_bytes()).hexdigest(),'imap_sha256':hashlib.sha256(imap.read_bytes()).hexdigest(),
              'role_scope':'Ghost targets recipient. Inflow: out-clade -> sister. Outflow: sister -> out-clade. Default priors match the project plant baseline; choose scientifically appropriate priors before substantive inference.',
              'ghost_semantics':'Two tau-parent=yes occurrences encode an unsampled private route with its own attachment times; this is not a bidirectional exchange with two independent phi parameters.',
              'validation':'Internal graph roundtrip. BPP runtime parse validation is a separate development check, not performed automatically.'}
    (outdir/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return {'suite':str(outdir),'n_controls':len(tasks),'networks':networks,'next_command':'sh '+shlex.quote(str(script)),'scope':manifest['role_scope']}
