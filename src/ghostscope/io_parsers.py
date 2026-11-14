"""Small, strict tree/alignment parsers and a JC69 + NJ reconstruction path."""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import gzip
import json
import math
import re
import numpy as np

EPS = 1e-12

@dataclass(eq=False)
class Node:
    name: str = ''
    length: Optional[float] = None
    children: list = field(default_factory=list)
    annotations: list = field(default_factory=list)

    def tips(self):
        return [self.name] if not self.children else [x for c in self.children for x in c.tips()]

def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz' else path.open(encoding='utf-8')

def parse_newick(text):
    """Quoted labels (including escaped quotes), nested comments, scientific lengths."""
    text = text.strip()
    i = 0
    def spaces():
        nonlocal i
        while i < len(text) and text[i].isspace(): i += 1
    def comments():
        nonlocal i
        items = []
        spaces()
        while i < len(text) and text[i] == '[':
            start = i; depth = 0
            while i < len(text):
                ch = text[i]; i += 1
                if ch == '[': depth += 1
                if ch == ']':
                    depth -= 1
                    if depth == 0: break
            if depth: raise ValueError('Unterminated Newick comment')
            items.append(text[start:i]); spaces()
        return items
    def label():
        nonlocal i
        spaces()
        if i < len(text) and text[i] in "'\"":
            quote = text[i]; i += 1; chars = []
            while i < len(text):
                ch = text[i]; i += 1
                if ch == quote:
                    if i < len(text) and text[i] == quote:
                        chars.append(ch); i += 1; continue
                    return ''.join(chars)
                chars.append(ch)
            raise ValueError('Unterminated quoted Newick label')
        start = i
        while i < len(text) and text[i] not in ':,();[]' and not text[i].isspace(): i += 1
        return text[start:i]
    def subtree():
        nonlocal i
        spaces(); prefix = comments(); children = []
        if i < len(text) and text[i] == '(':
            i += 1
            children.append(subtree()); spaces()
            while i < len(text) and text[i] == ',':
                i += 1; children.append(subtree()); spaces()
            if i >= len(text) or text[i] != ')': raise ValueError('Expected closing Newick parenthesis')
            i += 1
        name = label(); ann = prefix + comments(); length = None
        if i < len(text) and text[i] == ':':
            i += 1; spaces(); start = i
            while i < len(text) and text[i] not in ',();[]' and not text[i].isspace(): i += 1
            try: length = float(text[start:i])
            except ValueError: raise ValueError('Invalid Newick branch length') from None
            if not math.isfinite(length): raise ValueError('Nonfinite Newick branch length')
            ann += comments()
        if not name and not children: raise ValueError('Empty Newick tip')
        return Node(name, length, children, ann)
    root = subtree(); spaces()
    if i < len(text) and text[i] == ';': i += 1
    spaces()
    if i != len(text): raise ValueError('Trailing Newick content')
    if len(root.tips()) != len(set(root.tips())): raise ValueError('Duplicate tip labels')
    return root

def quote_label(name):
    return "'" + name.replace("'", "''") + "'" if re.search(r"[\s:,();\[\]'\"]", name) else name

def newick(root):
    def rec(n):
        s = '(' + ','.join(rec(c) for c in n.children) + ')' if n.children else ''
        s += quote_label(n.name) + ''.join(n.annotations)
        if n.length is not None: s += ':' + format(n.length, '.15g')
        return s
    return rec(root) + ';'

def split_newicks(stream):
    buf = []; quote = None; depth = 0
    # Quote state handles doubled quotes because two consecutive toggles restore it.
    for chunk in stream:
        for ch in chunk:
            buf.append(ch)
            if quote:
                if ch == quote: quote = None
            elif ch in "'\"" and depth == 0: quote = ch
            elif ch == '[': depth += 1
            elif ch == ']': depth -= 1
            elif ch == ';' and depth == 0:
                yield ''.join(buf).strip(); buf = []
    if ''.join(buf).strip(): raise ValueError('Newick stream ends without semicolon')

def genealogy_tree(record):
    nodes = record['nodes']; visiting = set()
    def rec(i, length=None):
        if i in visiting: raise ValueError('Cycle in genealogy')
        visiting.add(i); row = nodes[i]
        result = Node(row.get('label', ''), length,
                      [rec(int(j), float(d)) for j, d in row['children']])
        visiting.remove(i)
        return result
    root = rec(int(record['root']))
    if len(root.tips()) != len(set(root.tips())): raise ValueError('Duplicate genealogy tips')
    return root

def read_trees(path):
    p = Path(path)
    if p.is_dir():
        files = sorted(x for x in p.iterdir() if x.is_file() and any(str(x).endswith(s) for s in ('.nwk','.tre','.tree','.trees','.txt','.nwk.gz','.tre.gz','.jsonl','.jsonl.gz')))
        candidates=[]
        for f in files:
            if f.suffix=='.txt':
                with open_text(f) as stream: beginning=stream.read(256).lstrip()
                if not beginning.startswith(('(','[','#NEXUS','#nexus')):continue
            candidates.append(f)
        files=candidates
        if not files: raise ValueError('No tree files found')
        for f in files: yield from read_trees(f)
    else:
        with open_text(p) as f:
            if str(p).endswith(('.jsonl', '.jsonl.gz')):
                for line in f:
                    if line.strip(): yield genealogy_tree(json.loads(line))
            else:
                first=f.readline()
                if first.lstrip().upper().startswith('#NEXUS'):
                    # Author PhyloNet files contain one named TREE per line. Numeric
                    # TRANSLATE blocks need a dedicated map and are rejected explicitly.
                    for line in f:
                        if line.lstrip().upper().startswith('TRANSLATE'):raise ValueError('NEXUS TRANSLATE blocks are not supported; export labelled Newick')
                        match=re.match(r'\s*tree\s+[^=]+=(.*)',line,re.I)
                        if match:yield parse_newick(match.group(1).strip())
                else:
                    import itertools
                    for text in split_newicks(itertools.chain([first],f)): yield parse_newick(text)

def load_tree(value):
    value=str(value)
    return parse_newick(value) if '(' in value else next(read_trees(value))

def clone(n):
    return Node(n.name, n.length, [clone(c) for c in n.children], list(n.annotations))

def induced(root, keep):
    keep = set(keep)
    if not keep <= set(root.tips()): raise ValueError('Requested tips absent from tree')
    def rec(n):
        if not n.children: return clone(n) if n.name in keep else None
        kids = [v for c in n.children for v in [rec(c)] if v is not None]
        if not kids: return None
        if len(kids) == 1:
            v = kids[0]
            v.length = None if v.length is None or n.length is None else v.length + n.length
            return v
        return Node(n.name, n.length, kids)
    result = rec(root); result.length = None
    return result

def canonical(root, eps=EPS):
    """Rooted unordered topology, contracting explicitly zero internal edges."""
    def rec(n):
        if not n.children: return n.name
        children = []
        for c in n.children:
            value = rec(c)
            if c.children and c.length is not None and c.length <= eps:
                children.extend(value if isinstance(value, tuple) else (value,))
            else: children.append(value)
        if len(children) == 1: return children[0]
        return tuple(sorted(children, key=repr))
    return rec(root)

def binary(topology):
    return isinstance(topology, str) or len(topology) == 2 and all(binary(t) for t in topology)

def unroot_key(root, eps=EPS):
    """Nontrivial bipartitions; same tie-breaking as the exact six-leaf certificate."""
    all_tips = set(root.tips()); splits = {}
    def rec(n):
        sub = {n.name} if not n.children else set().union(*(rec(c) for c in n.children))
        other = all_tips - sub
        if min(len(sub), len(other)) > 1:
            key = tuple(sorted(min([sub, other], key=lambda z:(len(z),sorted(z)))))
            old = splits.get(key, 0.)
            splits[key] = None if n.length is None or old is None else old + n.length
        return sub
    rec(root)
    return tuple(sorted(k for k,v in splits.items() if v is None or v > eps))

def patristic(root, names):
    """Full paths via LCA; never use pendant lengths as pruned distances."""
    paths = {}
    def walk(n, lineage, depth):
        path = lineage + [(n, depth)]
        if not n.children: paths[n.name] = path
        for c in n.children:
            if c.length is None: raise ValueError('Distance analysis requires branch lengths on every path')
            if c.length < -EPS: raise ValueError('Negative branch length in distance input')
            walk(c, path, depth + max(0., c.length))
    walk(root, [], 0.)
    out = np.zeros((len(names), len(names)))
    for i,a in enumerate(names):
        if a not in paths: raise ValueError('Missing tip: ' + a)
        for j in range(i):
            b = names[j]; lca = 0.
            for (na,da),(nb,db) in zip(paths[a], paths[b]):
                if na is not nb: break
                lca = da
            out[i,j] = out[j,i] = paths[a][-1][1] + paths[b][-1][1] - 2*lca
    return out

def reroot_outgroup(root, outgroups):
    """Root on a monophyletic outgroup edge; split its length exactly in half."""
    target = set(outgroups); all_tips = set(root.tips())
    if not target or not target < all_tips: raise ValueError('Outgroup must be a proper nonempty subset')
    graph = {}; descendants = {}
    def walk(n):
        graph.setdefault(n, [])
        sub = {n.name} if not n.children else set()
        for c in n.children:
            graph.setdefault(c, []); graph[n].append((c,c.length)); graph[c].append((n,c.length))
            sub |= walk(c)
        descendants[n] = sub
        return sub
    walk(root); edge = None
    for n,adj in graph.items():
        for c,length in adj:
            if c in n.children and (descendants[c] == target or all_tips - descendants[c] == target): edge = (n,c,length); break
        if edge: break
    if edge is None: raise ValueError('Outgroup is not monophyletic on a single edge')
    a,b,length = edge
    def orient(n, previous, edge_length):
        children = [orient(v,n,d) for v,d in graph[n] if v is not previous]
        if len(children) == 1:
            child = children[0]
            child.length = None if edge_length is None or child.length is None else edge_length + child.length
            return child
        return Node(n.name if not children else '', edge_length, children)
    half = None if length is None else length / 2
    return Node(children=[orient(a,b,half),orient(b,a,half)])

class TaxonMapper:
    def __init__(self, imap=None, delimiter=None, known=()):
        self.mapping = {}; self.delimiter = delimiter; self.known = set(known)
        if imap:
            with open_text(imap) as f:
                for line in f:
                    if not line.strip() or line.lstrip().startswith('#'): continue
                    fields = line.split()
                    if len(fields) != 2: raise ValueError('imap requires two columns: sample species')
                    a,b = fields
                    if a in self.mapping and self.mapping[a] != b: raise ValueError('Conflicting imap mapping')
                    self.mapping[a] = b
    def __call__(self, name):
        if name in self.mapping: return self.mapping[name]
        if '^' in name and name.rsplit('^',1)[1] in self.mapping: return self.mapping[name.rsplit('^',1)[1]]
        if name in self.known: return name
        if self.delimiter: return name.split(self.delimiter,1)[0]
        return name
    def tree(self, root):
        out = clone(root)
        def walk(n):
            if not n.children: n.name = self(n.name)
            for c in n.children: walk(c)
        walk(out)
        if len(out.tips()) != len(set(out.tips())): raise ValueError('Multiple samples per species: choose representatives explicitly for topology audits')
        return out

def read_alignments(path, fmt='auto'):
    p = Path(path)
    if p.is_dir():
        files = sorted(f for f in p.iterdir() if f.is_file() and any(str(f).endswith(s) for s in ('.fa','.fas','.fasta','.phy','.phylip','.bpp','.fa.gz','.fasta.gz','.phy.gz')))
        if not files: raise ValueError('No alignment files found')
        for f in files:
            for tag,seq in read_alignments(f, fmt): yield f.name + ':' + str(tag), seq
        return
    with open_text(p) as f: lines = [x.strip() for x in f if x.strip()]
    if not lines: raise ValueError('Empty alignment input')
    if lines[0].startswith('>'):
        seq = {}; name = None
        for line in lines:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in seq: raise ValueError('Duplicate FASTA label')
                seq[name] = ''
            elif name is None: raise ValueError('Missing FASTA header')
            else: seq[name] += ''.join(line.split()).upper()
        _check_sequences(seq); yield 0,seq
        return
    i = 0; locus = 0
    while i < len(lines):
        header = lines[i].split(); i += 1
        if len(header) != 2 or not all(x.isdigit() for x in header): raise ValueError('Expected PHYLIP/BPP nseq nsites header')
        n,L = map(int, header)
        if n < 2 or L < 1: raise ValueError('Invalid alignment dimensions')
        seq = {}; order = []; interleaved = None
        # Both relaxed and classic PHYLIP first blocks; sequential wrapped rows too.
        for _ in range(n):
            if i >= len(lines): raise ValueError('Truncated PHYLIP block')
            row = lines[i].split(); i += 1
            if len(row) < 2:
                raw = lines[i-1]; row = [raw[:10].strip(),raw[10:].replace(' ','')]
            name = row[0]; value = ''.join(row[1:]).upper()
            if name in seq or not name: raise ValueError('Duplicate/empty PHYLIP sample')
            if interleaved is None and len(value)<L and i<len(lines):
                interleaved=len(lines[i].split())>=2
            while not interleaved and len(value) < L and i < len(lines) and len(lines[i].split()) == 1:
                value += lines[i].upper(); i += 1
            seq[name] = value; order.append(name)
        while any(len(seq[name]) < L for name in order):
            for name in order:
                if i >= len(lines): raise ValueError('Truncated interleaved PHYLIP')
                row = lines[i].split(); i += 1
                if row and row[0] == name: row = row[1:]
                seq[name] += ''.join(row).upper()
        _check_sequences(seq, L); yield locus,seq; locus += 1

def _check_sequences(seq, length=None):
    lengths = {len(x) for x in seq.values()}
    if len(seq)<2 or len(lengths) != 1 or not next(iter(lengths)) or length is not None and lengths != {length}:
        raise ValueError('Alignment sequences have invalid or unequal lengths')

def distance_matrix(seq, metric='jc69'):
    names = list(seq); _check_sequences(seq)
    x = np.stack([np.frombuffer(seq[n].upper().encode('ascii'), dtype=np.uint8) for n in names])
    valid = np.isin(x, np.frombuffer(b'ACGT',dtype=np.uint8))
    d = np.zeros((len(names),len(names))); counts = np.zeros_like(d, dtype=int)
    for i in range(len(names)):
        for j in range(i):
            keep = valid[i] & valid[j]; n = int(keep.sum()); counts[i,j] = counts[j,i] = n
            p = float(np.mean(x[i,keep] != x[j,keep])) if n else float('nan')
            value = p if metric == 'p' else (-.75*np.log1p(-4*p/3) if n and p<.75 else float('nan'))
            d[i,j] = d[j,i] = value
    return names,d,counts

def neighbor_joining(names, distances):
    """Standard NJ. Nonpositive fitted edges are collapsed, never randomized."""
    d = np.array(distances, dtype=float, copy=True)
    if d.shape != (len(names),len(names)) or not np.isfinite(d).all() or not np.allclose(d,d.T) or (d < -EPS).any():
        raise ValueError('NJ requires a finite symmetric nonnegative distance matrix')
    if len(set(names)) != len(names) or len(names)<2: raise ValueError('NJ requires distinct tip names')
    nodes = [Node(n) for n in names]
    while len(nodes)>2:
        n = len(nodes); sums = d.sum(1); q = (n-2)*d-sums[:,None]-sums[None,:]
        np.fill_diagonal(q,np.inf); i,j = np.unravel_index(np.argmin(q),q.shape)
        if i>j: i,j=j,i
        li = .5*d[i,j] + (sums[i]-sums[j])/(2*(n-2)); lj = d[i,j]-li
        nodes[i].length = max(0.,float(li)); nodes[j].length = max(0.,float(lj))
        joined = Node(children=[nodes[i],nodes[j]])
        keep = [k for k in range(n) if k not in (i,j)]
        nd = np.zeros((n-1,n-1)); nd[:-1,:-1] = d[np.ix_(keep,keep)]
        nd[-1,:-1] = nd[:-1,-1] = np.maximum(0.,(d[i,keep]+d[j,keep]-d[i,j])/2)
        nodes = [nodes[k] for k in keep]+[joined]; d = nd
    nodes[0].length = nodes[1].length = max(0.,float(d[0,1]/2))
    return Node(children=nodes)
