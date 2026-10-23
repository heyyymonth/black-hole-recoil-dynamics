"""Build the REVTeX manuscript, portable LaTeX archive, and provenance manifest."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
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
VERSION = '0.3.0-manuscript.5'
FIGURE_NAMES = ('return-time.pdf', 'core-residence.pdf', 'exposure-comparison.pdf',
                'encounter-proxy.pdf', 'snapshot-comparison.pdf', 'catalog-counts.pdf')


def validate_reused_figures(root):
    """Check saved figures and their declared scientific inputs without imports
    of the dynamics or figure renderer, numerical integration, or output writes.

    The previous build manifest is the reference. Manuscript text, citations,
    and build machinery may change; scientific code, recorded results,
    configurations, protocols, numerical dependencies, and plot code may not.
    """
    previous_bytes = (root / 'output/pdf/build.json').read_bytes()
    previous = json.loads(previous_bytes)
    bound = previous['input_sha256']
    scientific = {
        name: digest for name, digest in bound.items()
        if name.startswith(('src/', 'tests/', 'experiments/', 'results/'))
        or (name.startswith('docs/') and name.endswith('-protocol.md'))
        or name in ('requirements-lock.txt', 'requirements.txt', 'pyproject.toml',
                    'paper/figures.py')
    }
    required = {'requirements-lock.txt', 'paper/figures.py',
                'experiments/config.json', 'src/recoil_dynamics/plummer.py',
                'src/recoil_dynamics/__init__.py', 'results/run_all.json',
                'results/run_e003.json', 'results/run_e004.json',
                'results/run_e005.json'}
    if missing := required - scientific.keys():
        raise SystemExit(f'Reuse reference is missing scientific bindings: {sorted(missing)}')
    figures = {f'paper/figures/{name}': bound.get(f'paper/figures/{name}')
               for name in FIGURE_NAMES}
    actual_names = {path.name for path in (root / 'paper/figures').glob('*.pdf')}
    if actual_names != set(FIGURE_NAMES):
        raise SystemExit('Figure reuse requires exactly the six declared manuscript PDFs')
    for name, expected in {**scientific, **figures}.items():
        path = root / name
        if not expected or not path.is_file():
            raise SystemExit(f'Missing figure-reuse input or binding: {name}')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise SystemExit(f'Figure reuse rejected: changed input {name}; regenerate figures explicitly')
    for name in figures:
        figure = PdfReader(root / name, strict=True)
        if len(figure.pages) != 1:
            raise SystemExit(f'Expected one-page vector figure: {name}')
        page = figure.pages[0]
        if float(page.mediabox.width) <= 0 or float(page.mediabox.height) <= 0:
            raise SystemExit(f'Invalid figure dimensions: {name}')
    return {'mode': 'reuse-verified', 'physical_calculations_performed': False,
            'reference_manifest_sha256': hashlib.sha256(previous_bytes).hexdigest(),
            'reference_manuscript_version': previous['version'],
            'figure_sha256': figures, 'scientific_input_sha256': scientific}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reuse-figures', action='store_true',
                    help='Verify and reuse existing figures; do not run the numerical figure renderer')
parser.add_argument('--check-reuse-figures', action='store_true',
                    help='Only validate the figure-reuse guard; write no files and perform no build')
args = parser.parse_args()
figure_generation = {'mode': 'regenerated', 'physical_calculations_performed': True}
if args.reuse_figures or args.check_reuse_figures:
    figure_generation = validate_reused_figures(ROOT)
    print(f'Validated six saved figures and {len(figure_generation["scientific_input_sha256"])} '
          'scientific input hashes; no physical calculations performed.')
if args.check_reuse_figures:
    raise SystemExit(0)

title_match = re.search(r'\\title\{([^{}]+)\}', (PAPER/'manuscript.tex').read_text())
if title_match is None:
    raise SystemExit('Expected one plain-text manuscript title for provenance.')
TITLE = re.sub(r'\s+', ' ', title_match.group(1)).strip()
for directory in (OUT, TMP, STAGE, PAPER / 'tables'):
    directory.mkdir(parents=True, exist_ok=True)


def rows(relative):
    with (ROOT / relative).open() as stream:
        return list(csv.DictReader(stream))


def write_table(name, header, values, alignment=None):
    content = ['% Generated from committed results by paper/build.py.',
               r'\begin{ruledtabular}', r'\begin{tabular}{'+(alignment or 'r'*len(header))+'}',
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
encounter_rows = rows('results/e003/encounters.csv')
def scientific(value):
    mantissa, exponent = f'{float(value):.3e}'.split('e')
    return f'${mantissa}\\times10^{{{int(exponent)}}}$'
write_table('encounters.tex', [r'$\eta$', '$D$', '$Q$', r'$\delta$ (\%)', r'$\tau$'],
            [[f"{float(r['eta']):g}", f"{float(r['density_ratio']):.5f}",
              f"{float(r['encounter_ratio']):.5f}", f"{100*float(r['proxy_relative_error']):.2f}",
              scientific(r['tau_total'])] for r in encounter_rows])
population_names = {'all_systems': 'All systems', 'all_binaries': 'All binaries',
                    'single_black_holes': 'Single BHs', 'bh_bearing_binaries': 'BH-bearing binaries',
                    'binary_black_holes': 'BBHs'}
write_table('population.tex', ['Population', '$N$', r'$\mathrm{median}(r)/r_h$'],
            [[population_names[r['population']], f"{int(r['count']):,}",
              f"{float(r['median_radius_over_host_rh']):.5f}"]
             for r in rows('results/e004/population.csv')], alignment='lrr')
snapshot_rows = rows('results/e004/trajectories.csv')
finest = max(int(r['shell_count']) for r in snapshot_rows)
write_table('snapshot-return.tex', [r'$\eta_{\mathrm{Pl}}$', '$v_k$',
            r'$P_{\mathrm{sh}}$', r'$P_{\mathrm{Pl}}$', r'$P_{\mathrm{Pl}}/P_{\mathrm{sh}}$'],
            [[f"{float(r['eta_plummer']):g}", f"{float(r['kick_kms']):.4f}",
              f"{float(r['snapshot_return_myr']):.5f}", f"{float(r['plummer_return_myr']):.5f}",
              f"{1+float(r['plummer_return_relative_error']):.4f}"]
             for r in snapshot_rows if int(r['shell_count']) == finest])
catalog = json.loads((ROOT/'results/e005/summary.json').read_text())
catalog_counts = catalog['selected_bbh_counts']
write_table('catalog-counts.tex', ['Quantity', 'Count'], [
    ['Available model histories', f"{catalog['model_histories']:,}"],
    ['Validated saved states', f"{catalog['snapshot_records']:,}"],
    ['Selected near 12 Gyr', f"{catalog['selected_models']:,}"],
    ['Histories without age coverage', f"{len(catalog['models_without_age_coverage']):,}"],
    [r'Selected $N_{\rm BBH}$: min/median/max',
     '/'.join(f"{catalog_counts[key]:g}" for key in ['minimum', 'median', 'maximum'])],
    [r'Selected states with $N_{\rm BBH}=0$', f"{catalog['selected_zero_bbh_models']:,}"],
    [r'Selected states with $N_{\rm BBH}\leq2$', f"{catalog['selected_at_most_two_bbh_models']:,}"],
    [rf"Selected states with $N_{{\rm BBH}}\geq{catalog['minimum_bbh_count_screen']}$",
     f"{catalog['selected_meeting_count_screen_models']:,}"],
    [r'Maximum $N_{\rm BBH}$ over saved states', f"{catalog['saved_state_bbh_count_range'][1]:,}"],
], alignment='lr')
if not args.reuse_figures:
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
    source_archive = OUT / 'manuscript-source.zip'
    if source_archive.exists():
        with zipfile.ZipFile(source_archive) as bundle:
            revtex_archive.write_bytes(bundle.read('vendor/revtex-4.2f.zip'))
    else:
        subprocess.run(['curl', '--fail', '--silent', '--show-error', '--location',
                        '--connect-timeout', '15', '--max-time', '120', '--retry', '2',
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
reader = PdfReader(TMP / 'manuscript.pdf')
text = '\n'.join(page.extract_text() for page in reader.pages)
for required in ['Venkata Hemanth Athota', 'bharadwaj.hemanth.002@gmail.com',
                 'Companion access after black-hole recoil']:
    if re.sub(r'\s+', '', required) not in re.sub(r'\s+', '', text):
        raise AssertionError(f'Missing manuscript text: {required}')
if re.search(r'E00[1-5]', re.sub(r'\s+', '', text), flags=re.IGNORECASE):
    raise AssertionError('An internal experiment identifier is visible in the manuscript PDF')
for page in reader.pages:
    assert abs(float(page.mediabox.width)-612) < .1
    assert abs(float(page.mediabox.height)-792) < .1
for forbidden in ['BHRD-TN-001', 'Not peer reviewed', 'Published as a conference paper', 'DRAFT']:
    assert forbidden not in text, f'Unexpected publication label: {forbidden}'
shutil.copyfile(TMP / 'manuscript.pdf', pdf)

archive = OUT / 'manuscript-source.zip'
source_files = [PAPER / 'manuscript.tex', PAPER / 'references.bib',
                PAPER / 'SOURCE_README.txt', PAPER / 'CITATION.bib']
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
# Bind all reported experiments to their configurations, numerical source,
# tables and summaries. External file identities stay pinned in the snapshot
# and catalog manifests; the catalog input inventory is also directly bound.
for name in ['run_all', 'run_e003', 'run_e004', 'run_e005']:
    path = ROOT/'results'/f'{name}.json'
    inputs.append(path)
    record = json.loads(path.read_text())
    for group in ['source_sha256', 'output_sha256']:
        for relative in record.get(group, {}):
            inputs.append(ROOT/relative)
for experiment in ['e001', 'e002', 'e003', 'e004', 'e005']:
    inputs += sorted((ROOT/'results'/experiment).glob('*.csv'))
inputs += [ROOT/'src/recoil_dynamics/plummer.py', ROOT/'src/recoil_dynamics/__init__.py']
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {'paper': 'central-residence', 'version': VERSION, 'pages': len(reader.pages),
            'title': TITLE,
            'experiments': ['E001', 'E002', 'E003', 'E004', 'E005'],
            'format': 'REVTeX 4.2 APS/PRD reprint; US Letter; two columns; txfonts',
            'engine': engine_version, 'source_date_epoch': env['SOURCE_DATE_EPOCH'],
            'figure_generation': figure_generation,
            'revtex_log': next(line.strip() for line in log.splitlines() if line.startswith('Document Class: revtex')),
            'versions': {name: importlib.metadata.version(name) for name in ['matplotlib', 'numpy', 'scipy', 'pypdf']},
            'checks': {'undefined_references': 0, 'overfull_boxes': 0, 'bibliography_entry_warnings': 0,
                       'visible_experiment_identifiers': 0,
                       'allowed_bst_diagnostics': [line for line in style_diagnostics if line in blg]},
            'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(inputs))},
            'output_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [pdf, archive]}}
(OUT / 'build.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(f'Built {len(reader.pages)}-page REVTeX manuscript and portable LaTeX source archive.')
