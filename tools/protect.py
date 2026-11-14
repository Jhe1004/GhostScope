"""Read-only historical verification; all records are written inside this package."""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
OUT = PACKAGE / 'validation_outputs'

def digest(p):
    if p.is_symlink():
        return 'symlink:' + os.readlink(p)
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def hashes(paths):
    def one(key):
        try:
            return key, digest(ROOT / key)
        except (OSError, ValueError) as error:
            return key, 'ERROR:' + str(error)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        return dict(pool.map(one, paths))

def main():
    start = time.monotonic()
    if not (ROOT / 'round8_456_boundaries_20261108/PROTECTED_BASELINE.json').exists():
        import csv
        mode = sys.argv[1] if len(sys.argv) > 1 else 'final'
        pub_tsv = ROOT / 'PUBLIC_ARCHIVE_SHA256.tsv'
        if not pub_tsv.is_file():
            err = {'status': 'FAIL', 'mode': 'PUBLIC_ARCHIVE_LIVE_CHECK', 'reason': 'PUBLIC_ARCHIVE_SHA256.tsv not found at archive root'}
            print(json.dumps(err, ensure_ascii=False), file=sys.stderr, flush=True)
            raise SystemExit(1)
        dynamic_outputs = {
            '09_GhostScope_Toolkit_v1/validation_outputs/VALIDATION_CERTIFICATE.json',
            '09_GhostScope_Toolkit_v1/validation_outputs/PAPER_EMPIRICAL_REPLAY.json',
            '09_GhostScope_Toolkit_v1/validation_outputs/PAPER_EMPIRICAL_SUMMARY.tsv',
            '09_GhostScope_Toolkit_v1/validation_outputs/POSITIVE_CONTROLS_REPLAY.tsv',
            '09_GhostScope_Toolkit_v1/validation_outputs/PROTECTION_PUBLIC_ARCHIVE_LIVE.json',
            '08_Higher_Order_4to6_and_ArbitraryN_R8/certificates/ARBITRARY_N_SYMBOLIC_VALIDATION.json',
            '08_Higher_Order_4to6_and_ArbitraryN_R8/certificates/ARBITRARY_N_FINAL_VALIDATION.json',
        }
        expected = {}
        with pub_tsv.open('r', encoding='utf-8') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                rp = row['relative_path']
                if rp not in dynamic_outputs:
                    expected[rp] = row['sha256']
        current = hashes(expected)
        changed_or_removed = sorted(p for p, h in expected.items() if current.get(p) != h)
        replay_cert = OUT / 'VALIDATION_CERTIFICATE.json'
        replay_inputs = {}
        if replay_cert.is_file():
            try:
                replay_inputs = json.loads(replay_cert.read_text()).get('input_sha256', {})
            except Exception:
                replay_inputs = {}
        session_changed = sorted(p for p in replay_inputs if p in expected and current.get(p) != expected[p])
        allowed_all = set(expected) | dynamic_outputs | {'PUBLIC_ARCHIVE_SHA256.tsv'}
        added_outside = []
        for parent, dirs, names in os.walk(ROOT):
            dirs[:] = [d for d in dirs if d != '__pycache__' and not d.startswith('.')]
            for name in names:
                if name == '.DS_Store' or name.endswith('.pyc'):
                    continue
                rel_p = (Path(parent) / name).relative_to(ROOT).as_posix()
                if rel_p not in allowed_all:
                    added_outside.append(rel_p)
        added_outside.sort()
        ok = (len(expected) > 0) and (not changed_or_removed) and (not session_changed) and (not added_outside)
        record = {
            'status': 'PASS' if ok else 'FAIL',
            'verification_mode': 'PUBLIC_ARCHIVE_SHA256_LIVE_RECOMPUTATION',
            'phase': mode,
            'manifest_path': 'PUBLIC_ARCHIVE_SHA256.tsv',
            'manifest_sha256': digest(pub_tsv),
            'baseline_entries': len(expected),
            'baseline_changed_or_removed': changed_or_removed,
            'session_files': len(expected),
            'session_changed_or_removed': session_changed,
            'session_added_outside_package': added_outside,
            'seconds': time.monotonic() - start,
        }
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / 'PROTECTION_PUBLIC_ARCHIVE_LIVE.json').write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')
        print(json.dumps(record, ensure_ascii=False), flush=True)
        if not ok:
            raise SystemExit(1)
        return
    baseline_path = ROOT / 'round8_456_boundaries_20261108/PROTECTED_BASELINE.json'
    baseline = json.loads(baseline_path.read_text())
    current = hashes(baseline)
    changes = [p for p, h in baseline.items() if current[p] != h]
    mode = sys.argv[1]
    snap_path = OUT / 'PREDEVELOPMENT_HASHES.json'
    if mode == 'initial':
        if snap_path.exists():
            raise SystemExit('Refusing to overwrite predevelopment snapshot')
        paths = []
        for parent, dirs, files in os.walk(ROOT):
            if Path(parent) == ROOT:
                dirs[:] = [d for d in dirs if d != PACKAGE.name]
            paths.extend(str((Path(parent) / f).relative_to(ROOT)) for f in files)
        snapshot = dict(current)
        snapshot.update(hashes([p for p in paths if p not in snapshot]))
        snap_path.write_text(json.dumps(snapshot, ensure_ascii=False, sort_keys=True))
        session_changes = []
    else:
        before = json.loads(snap_path.read_text())
        snapshot = hashes(before)
        session_changes = [p for p, h in before.items() if snapshot[p] != h]
    additions=[]
    if mode!='initial':
        for parent,dirs,names in os.walk(ROOT):
            if Path(parent)==ROOT:dirs[:]=[d for d in dirs if d not in (PACKAGE.name,'manuscript_package_SciAdv')]
            additions.extend(str((Path(parent)/name).relative_to(ROOT)) for name in names if str((Path(parent)/name).relative_to(ROOT)) not in before)
    record = {'phase': mode, 'baseline_entries': len(baseline),
              'baseline_sha256': digest(baseline_path),
              'baseline_changed_or_removed': changes,
              'session_files': len(snapshot), 'session_changed_or_removed': session_changes,
              'session_added_outside_package':additions,
              'v2_files': sum(p.startswith('manuscript_package_SysBio_v2/') for p in snapshot),
              'v2_changed': [p for p in session_changes if p.startswith('manuscript_package_SysBio_v2/')],
              'seconds': time.monotonic() - start}
    (OUT / ('PROTECTION_' + mode.upper() + '.json')).write_text(json.dumps(record, ensure_ascii=False, indent=2))
    print(json.dumps(record, ensure_ascii=False), flush=True)

if __name__ == '__main__':
    main()
