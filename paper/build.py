"""Build the REVTeX manuscript, portable LaTeX archive, and provenance manifest."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'paper'
OUT = ROOT / 'output/pdf'
TMP = ROOT / 'tmp/latex'
STAGE = TMP / 'source'
TOOLCHAIN = json.loads((PAPER / 'toolchain.json').read_text())
VERSION = '0.1.1'
for directory in (OUT, TMP, STAGE, PAPER / 'tables'):
    directory.mkdir(parents=True, exist_ok=True)


def rows(relative):
    with (ROOT / relative).open() as stream:
        return list(csv.DictReader(stream))


def write_table(name, header, values):
    content = ['% Generated from committed CSV results by paper/build.py.',
               r'\begin{ruledtabular}', r'\begin{tabular}{'+'r'*len(header)+'}',
               ' & '.join(header)+r' \\', r'\hline']
    content.extend(' & '.join(row)+r' \\' for row in values)
    content += [r'\end{tabular}', r'\end{ruledtabular}']
    (PAPER / 'tables' / name).write_text('\n'.join(content)+'\n')


return_rows = [row for row in rows('results/e001/return_times.csv')
               if float(row['eta']) in [.5, .8, .9, .95, .99, .999]]
write_table('return.tex', [r'$\eta$', r'$x_a$', r'$P$'],
            [[f"{float(row['eta']):g}", f"{float(row['apocenter_over_b']):.6f}",
              f"{float(row['first_return_over_tdyn']):.6f}"] for row in return_rows])
cycle_rows = [row for row in rows('results/e002/cycles.csv') if float(row['eta']) in [.3, .8, .99]]
write_table('exposure.tex', [r'$\eta$', r'$I(P)$', r'$I(P)/P$', r'$C(P,1)/P$'],
            [[f"{float(row['eta']):g}", f"{float(row['cycle_density_exposure']):.6f}",
              f"{float(row['cycle_mean_density']):.6g}", f"{float(row['core_fraction_c_1']):.6g}"]
             for row in cycle_rows])
subprocess.run([sys.executable, str(PAPER / 'figures.py')], cwd=ROOT, check=True)
engine = shutil.which('tectonic')
if not engine:
    raise SystemExit('Install Tectonic 0.16.9, then rerun. See paper/README.md.')
engine_version = subprocess.check_output([engine, '--version'], text=True).strip()
if engine_version != 'Tectonic 0.16.9':
    raise SystemExit(f'Expected Tectonic 0.16.9; found {engine_version}.')
env = os.environ.copy()
env['SOURCE_DATE_EPOCH'] = str(int(datetime(2026, 10, 18, tzinfo=timezone.utc).timestamp()))
revtex_archive = ROOT / 'tmp/revtex-4.2f.zip'
if not revtex_archive.exists():
    subprocess.run(['curl', '--fail', '--silent', '--show-error', '--location',
                    TOOLCHAIN['revtex_archive_url'], '--output', str(revtex_archive)], check=True)
if hashlib.sha256(revtex_archive.read_bytes()).hexdigest() != TOOLCHAIN['revtex_archive_sha256']:
    raise SystemExit('REVTeX dependency checksum mismatch; do not use an unverified replacement.')
with zipfile.ZipFile(revtex_archive) as distribution:
    for item in distribution.infolist():
        if item.filename.startswith(('tex/latex/revtex/', 'bibtex/bst/revtex/')) and not item.is_dir():
            (STAGE / Path(item.filename).name).write_bytes(distribution.read(item))
for name in ['manuscript.tex', 'references.bib']:
    shutil.copyfile(PAPER / name, STAGE / name)
for folder in ['figures', 'tables']:
    shutil.copytree(PAPER / folder, STAGE / folder, dirs_exist_ok=True)
command = [engine, '--untrusted', '--keep-logs', '--keep-intermediates',
           '--outdir', str(TMP), 'manuscript.tex']
completed = subprocess.run(command, cwd=STAGE, env=env, text=True, capture_output=True)
(TMP / 'build-console.txt').write_text(completed.stdout+completed.stderr)
if completed.returncode:
    raise SystemExit(completed.stdout+completed.stderr)
log = (TMP / 'manuscript.log').read_text(errors='replace')
blg = (TMP / 'manuscript.blg').read_text(errors='replace')
for problem in [r'Overfull \\[hv]box', r'Citation .+ undefined', r'Reference .+ undefined',
                r'There were undefined', r'Missing character:', r'! LaTeX Error']:
    if re.search(problem, log):
        raise SystemExit(f'LaTeX validation failed: {problem}. Inspect tmp/latex/manuscript.log.')
style_diagnostics = ['Warning--jnrlst (dependency: not reversed) set 1']
entry_warnings = [line for line in blg.splitlines() if line.startswith('Warning--') and line not in style_diagnostics]
if entry_warnings or "I couldn't open" in blg:
    raise SystemExit('BibTeX warnings found; inspect tmp/latex/manuscript.blg.')
pdf = OUT / 'central-residence.pdf'
shutil.copyfile(TMP / 'manuscript.pdf', pdf)
reader = PdfReader(pdf)
text = '\n'.join(page.extract_text() for page in reader.pages)
for required in ['Venkata Hemanth Athota', 'bharadwaj.hemanth.002@gmail.com', 'Islam']:
    if re.sub(r'\s+', '', required) not in re.sub(r'\s+', '', text):
        raise AssertionError(f'Missing manuscript text: {required}')
for page in reader.pages:
    assert abs(float(page.mediabox.width)-612) < .1
    assert abs(float(page.mediabox.height)-792) < .1
for forbidden in ['BHRD-TN-001', 'Not peer reviewed', 'Published as a conference paper', 'DRAFT']:
    assert forbidden not in text, f'Unexpected publication label: {forbidden}'

archive = OUT / 'manuscript-source.zip'
source_files = [PAPER / 'manuscript.tex', PAPER / 'references.bib', PAPER / 'SOURCE_README.txt']
source_files += sorted((PAPER / 'figures').glob('*.pdf'))
source_files += sorted((PAPER / 'tables').glob('*.tex'))
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
    def add(name, content):
        info = zipfile.ZipInfo(name, date_time=(2026, 10, 18, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = 0o100644 << 16
        bundle.writestr(info, content)
    for path in source_files:
        add(path.relative_to(PAPER).as_posix(), path.read_bytes())
    add('manuscript.bbl', (TMP / 'manuscript.bbl').read_bytes())
    add('LICENSE', (ROOT / 'LICENSE').read_bytes())
    add('vendor/revtex-4.2f.zip', revtex_archive.read_bytes())

inputs = [ROOT / name for name in ['paper/build.py', 'paper/verify.py', 'paper/figures.py',
          'paper/requirements.txt', 'paper/toolchain.json', 'requirements-lock.txt', 'CITATION.cff', 'LICENSE',
          'experiments/config.json', 'results/run_all.json']]
inputs += source_files
inputs += sorted((ROOT / 'results').glob('e*/*.csv'))
inputs += sorted((ROOT / 'src').rglob('*.py'))
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {'paper': 'central-residence', 'version': VERSION, 'pages': len(reader.pages),
            'format': 'REVTeX 4.2 APS/PRD reprint; US Letter; two columns; txfonts',
            'engine': engine_version, 'source_date_epoch': env['SOURCE_DATE_EPOCH'],
            'revtex_log': next(line.strip() for line in log.splitlines() if line.startswith('Document Class: revtex')),
            'versions': {name: importlib.metadata.version(name) for name in ['matplotlib', 'numpy', 'scipy', 'pypdf']},
            'checks': {'undefined_references': 0, 'overfull_boxes': 0, 'bibliography_entry_warnings': 0,
                       'allowed_bst_diagnostics': [line for line in style_diagnostics if line in blg]},
            'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(inputs))},
            'output_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [pdf, archive]}}
(OUT / 'build.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(f'Built {len(reader.pages)}-page REVTeX manuscript and portable LaTeX source archive.')
