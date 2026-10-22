"""Verify committed manuscript provenance, page dimensions, and vector fonts."""
from pathlib import Path
import hashlib
import json
import re
import zipfile
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
meta = json.loads((ROOT / 'output/pdf/build.json').read_text())
for group in ['input_sha256', 'output_sha256']:
    for relative, expected in meta[group].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f'{group} mismatch: {relative}; rebuild the paper')

version = re.fullmatch(r'(\d+\.\d+\.\d+)-manuscript\.(\d+)', meta['version'])
assert version is not None, 'Unsupported manuscript revision format'
tag = f'manuscript-v{version.group(1)}.{version.group(2)}'
citation_url = f'https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/{tag}'
with zipfile.ZipFile(ROOT / 'output/pdf/manuscript-source.zip') as archive:
    for name in ['CITATION.bib', 'SOURCE_README.txt']:
        relative = f'paper/{name}'
        assert relative in meta['input_sha256'], f'Unbound archive input: {relative}'
        assert archive.namelist().count(name) == 1, f'Missing or duplicate archive member: {name}'
        assert archive.read(name) == (ROOT / relative).read_bytes(), f'Archive content mismatch: {name}'
    citation = archive.read('CITATION.bib').decode()
    source_readme = archive.read('SOURCE_README.txt').decode()
    assert re.findall(r'revision\s+([^},\s]+)', citation) == [meta['version']], 'Citation revision mismatch'
    assert re.findall(r'\burl\s*=\s*\{([^{}]+)\}', citation) == [citation_url], 'Citation source tag mismatch'
    assert f'revision {meta["version"]}' in source_readme, 'Archive README revision mismatch'
    assert citation_url in source_readme, 'Archive README source tag mismatch'

fonts = set()
def verify_resources(resources):
    resources = resources.get_object()
    for ref in resources.get('/Font', {}).get_object().values() if '/Font' in resources else []:
        font = ref.get_object()
        assert font.get('/Subtype') != '/Type3', 'Bitmap/Type3 font found'
        for face in font.get('/DescendantFonts', [font]):
            face = face.get_object()
            descriptor = face.get('/FontDescriptor')
            assert descriptor, f'Missing font descriptor: {face.get("/BaseFont")}'
            descriptor = descriptor.get_object()
            assert any(key in descriptor for key in ['/FontFile', '/FontFile2', '/FontFile3']), 'Font is not embedded'
            fonts.add(str(face.get('/BaseFont')))
    if '/XObject' in resources:
        for ref in resources['/XObject'].get_object().values():
            obj = ref.get_object()
            assert obj.get('/Subtype') != '/Image', 'Raster figure or equation found'
            if '/Resources' in obj:
                verify_resources(obj['/Resources'])

reader = PdfReader(ROOT / 'output/pdf/central-residence.pdf')
assert len(reader.pages) == meta['pages']
text = '\n'.join(page.extract_text() for page in reader.pages)
assert not re.search(r'E00[1-5]', re.sub(r'\s+', '', text), re.IGNORECASE), 'Visible experiment identifier'
assert 'Companion access after black-hole recoil' in text.replace('\n', ' '), 'Incorrect paper title'
for required in ['Late-time availability', 'merger-conditioned accessibility', '7537', '180']:
    assert required.lower().replace(' ', '') in re.sub(r'\s+', '', text).lower(), required
for page in reader.pages:
    assert abs(float(page.mediabox.width)-612) < .1
    assert abs(float(page.mediabox.height)-792) < .1
    verify_resources(page['/Resources'])
print(f'Verified manuscript v{meta["version"]}: hashes, {len(reader.pages)} US Letter pages, '
      f'{len(fonts)} embedded vector fonts, no raster figures/equations, no visible experiment identifiers, '
      'and matching archive citation/README revisions and source tag.')
