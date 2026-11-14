"""Build/install the actual wheel and exercise the installed command, in isolation."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys
import time

P=Path(__file__).resolve().parents[1]
V=P/'validation_outputs'
ENV={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','TMPDIR':str(P/'.tmp'),
     'PIP_CACHE_DIR':str(P/'.cache'),'PIP_DISABLE_PIP_VERSION_CHECK':'1'}

def run(command,cwd,log,env=None):
    start=time.monotonic()
    r=subprocess.run(list(map(str,command)),cwd=cwd,env=env or ENV,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
    Path(log).write_text(r.stdout)
    return {'command':list(map(str,command)),'returncode':r.returncode,'seconds':time.monotonic()-start,'log':str(log)}

def main():
    python=P/'.venv/bin/python';cli=P/'.venv/bin/ghostscope'
    unit=run([sys.executable,'-m','unittest','discover','-s','tests','-v'],P,V/'UNIT_TESTS.log',dict(ENV,PYTHONPATH=str(P/'src')))
    text=(V/'UNIT_TESTS.log').read_text();match=re.search(r'Ran (\d+) tests',text)
    unit.update(tests=int(match.group(1)) if match else 0,skipped='skipped' in text,
                status='PASS' if unit['returncode']==0 and match and 'skipped' not in text else 'FAIL')
    (V/'UNIT_TEST_RESULT.json').write_text(json.dumps(unit,indent=2)+'\n')
    build=run([python,'-m','pip','wheel','.', '--no-deps','--no-build-isolation','-w','dist'],P,V/'WHEEL_BUILD.log')
    if build['returncode']:raise SystemExit('Wheel build failed; see WHEEL_BUILD.log')
    wheel=P/'dist/ghostscope-1.0.0-py3-none-any.whl'
    install=run([python,'-m','pip','install','--no-deps','--force-reinstall',wheel],P,V/'WHEEL_INSTALL.log')
    if install['returncode']:raise SystemExit('Wheel install failed; see WHEEL_INSTALL.log')
    runfolder=V/('cli_smoke_'+str(time.time_ns()));runfolder.mkdir()
    env=dict(ENV);env.pop('PYTHONPATH',None)
    rows=[]
    for command in (None,'trees','scan','seqs','mcmc','make-bpp','demo','reproduce-paper'):
        rows.append(run([cli]+([command] if command else [])+['--help'],runfolder,runfolder/((command or 'main')+'_help.log'),env))
    rows.append(run([cli,'demo','-o',runfolder/'demo'],runfolder,runfolder/'demo.log',env))
    data=P/'src/ghostscope/data/demo_bundle/ambiguous.nwk'
    rows.append(run([cli,'trees','-g',data,'--taxa','A,B,C','--outgroup','O','-o',runfolder/'trees'],runfolder,runfolder/'trees.log',env))
    source=P.parent/('08_Higher_Order_4to6_and_ArbitraryN_R8' if (P.parent/'08_Higher_Order_4to6_and_ArbitraryN_R8').exists() else 'round8_456_boundaries_20261108')
    dpath=source/'ingroup_4/4I-Exp0_r003';design=json.loads((dpath/'design.json').read_text());backbone=design['trees'][0]
    outgroups=[x[1] for x in design['labels'] if x[1] not in backbone.replace('(',' ').replace(')',' ').replace(',',' ').replace(';',' ').split()]
    outgroup=outgroups[0]
    rows.append(run([cli,'scan','-g',dpath/'genealogies.jsonl','-s',backbone,'--outgroup',outgroup,'--taxon-delim','_','-o',runfolder/'scan'],runfolder,runfolder/'scan.log',env))
    rows.append(run([cli,'seqs','-i',dpath/'loci','-s',backbone,'--outgroup',outgroup,'--taxon-delim','_','--nj-trees','-o',runfolder/'seqs'],runfolder,runfolder/'seqs.log',env))
    suite=V/'bpp_suite_v48'
    rows.append(run([cli,'mcmc','--chains',suite/'ghost/default/chain1/chain.mcmc.txt',suite/'ghost/default/chain2/chain.mcmc.txt','-o',runfolder/'mcmc'],runfolder,runfolder/'mcmc.log',env))
    mapping=[line.split() for line in (V/'bpp_input/imap.txt').read_text().splitlines() if line.strip()]
    species=list(dict.fromkeys(v[1] for v in mapping));a,b,c=species
    rows.append(run([cli,'make-bpp','-s',f'(({a},{b}),{c});','--recipient',a,'--sister',b,'--out-clade',c,
                     '--seqfile',V/'bpp_input/data.bpp','--imap',V/'bpp_input/imap.txt','--nloci','4','-o',runfolder/'bpp_suite'],runfolder,runfolder/'make_bpp.log',env))
    import_check=run([python,'-c',"import ghostscope,sys; print(ghostscope.__version__); print(ghostscope.__file__); print(sys.version)"],runfolder,runfolder/'installed_import.log',env)
    installed_text=(runfolder/'installed_import.log').read_text()
    record={'unit_tests':unit,'build':build,'install':install,'installed_import':import_check,
            'imported_from_site_packages':'site-packages' in installed_text,'cli_checks':rows,
            'all_commands_under_three_seconds':all(r['seconds']<3 for r in rows),
            'wheel_sha256':hashlib.sha256(wheel.read_bytes()).hexdigest(),'wheel_bytes':wheel.stat().st_size,
            'demo_bytes':sum(p.stat().st_size for p in (P/'src/ghostscope/data/demo_bundle').iterdir() if p.is_file()),
            'runtime_source_sha256':{str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((P/'src/ghostscope').rglob('*')) if p.is_file()},
            'scope':'Installed wheel checked outside the source import path. BPP smoke and scientific paper acceptance are separately recorded.'}
    record['status']='PASS' if unit['status']=='PASS' and all(r['returncode']==0 for r in rows) and record['imported_from_site_packages'] and record['all_commands_under_three_seconds'] else 'FAIL'
    (V/'ENGINEERING_VALIDATION.json').write_text(json.dumps(record,indent=2)+'\n')
    (V/'DEMO_OUTPUT.md').write_text((runfolder/'demo/report.md').read_text())
    print(json.dumps({'status':record['status'],'tests':unit['tests'],'cli_checks':len(rows),'wheel':str(wheel),'max_cli_seconds':max(r['seconds'] for r in rows)},indent=2))

if __name__=='__main__':main()
