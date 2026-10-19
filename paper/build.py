"""Render the repository's limited Markdown manuscript to a deterministic PDF."""
from pathlib import Path
from io import BytesIO
import hashlib
import importlib.metadata
import json
import re
from xml.sax.saxutils import escape

import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/pdf'
TMP = ROOT / 'tmp/pdfs/math'
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)
WIDTH = A4[0] - 108
FONTDIR = Path(matplotlib.get_data_path()) / 'fonts/ttf'
for name, file in [('Body', 'DejaVuSerif.ttf'), ('Body-Bold', 'DejaVuSerif-Bold.ttf'),
                   ('Body-Italic', 'DejaVuSerif-Italic.ttf'), ('Body-BoldItalic', 'DejaVuSerif-BoldItalic.ttf'),
                   ('Heading', 'DejaVuSans.ttf'), ('Heading-Bold', 'DejaVuSans-Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(FONTDIR / file)))
pdfmetrics.registerFontFamily('Body', normal='Body', bold='Body-Bold', italic='Body-Italic', boldItalic='Body-BoldItalic')
pdfmetrics.registerFontFamily('Heading', normal='Heading', bold='Heading-Bold')
styles = {
    'body': ParagraphStyle('body', fontName='Body', fontSize=9.5, leading=14, spaceAfter=8),
    'title': ParagraphStyle('title', fontName='Heading-Bold', fontSize=21, leading=27, spaceAfter=16),
    'h2': ParagraphStyle('h2', fontName='Heading-Bold', fontSize=13, leading=18, spaceBefore=8, spaceAfter=10, keepWithNext=True),
    'h3': ParagraphStyle('h3', fontName='Heading-Bold', fontSize=10.5, leading=15, spaceBefore=7, spaceAfter=6, keepWithNext=True),
    'caption': ParagraphStyle('caption', fontName='Body', fontSize=8, leading=11, spaceAfter=11),
    'cell': ParagraphStyle('cell', fontName='Body', fontSize=7.6, leading=11),
    'headcell': ParagraphStyle('headcell', fontName='Heading-Bold', fontSize=7.5, leading=10, textColor=colors.white),
}
matplotlib.rcParams.update({'mathtext.fontset': 'dejavuserif', 'font.family': 'DejaVu Serif'})


def equation(tex, size=10):
    key = hashlib.sha256(f'{size}:{tex}'.encode()).hexdigest()
    path = TMP / f'{key}.png'
    # Always render with current pinned libraries, rather than trust an old cache.
    buf = BytesIO()
    math_to_image('$'+tex+'$', buf, prop=FontProperties(size=size), dpi=220, format='png', color='#172b3a')
    path.write_bytes(buf.getvalue())
    with PILImage.open(path) as im:
        width, height = (v * 72 / 220 for v in im.size)
    return path, width, height


def inline(text):
    pattern = r'(\$[^$]+\$|\[[^\]]+\]\([^)]+\)|`[^`]+`|\*\*.+?\*\*|\*[^*]+\*)'
    parts = re.split(pattern, text)
    result = []
    for item in parts:
        if item.startswith('$') and item.endswith('$'):
            path, width, height = equation(item[1:-1], 9.5)
            result.append(f'<img src="{path}" width="{width:.3f}" height="{height:.3f}" valign="-3"/>')
        elif item.startswith('[') and '](' in item:
            label, url = re.match(r'\[([^\]]+)\]\(([^)]+)\)', item).groups()
            result.append(f'<link href="{escape(url)}" color="#14677d">{escape(label)}</link>')
        elif item.startswith('**'):
            result.append('<b>'+escape(item[2:-2])+'</b>')
        elif item.startswith('*'):
            result.append('<i>'+escape(item[1:-1])+'</i>')
        elif item.startswith('`'):
            result.append('<font name="Heading" size="8.5">'+escape(item[1:-1])+'</font>')
        else:
            result.append(escape(item))
    return ''.join(result)


def page(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#14677d'))
    canvas.setLineWidth(.7)
    canvas.line(54, A4[1]-39, A4[0]-54, A4[1]-39)
    canvas.setFont('Heading', 7)
    canvas.setFillColor(colors.HexColor('#465765'))
    canvas.drawString(54, A4[1]-29, 'BLACK-HOLE RECOIL DYNAMICS')
    canvas.drawRightString(A4[0]-54, A4[1]-29, 'BHRD-TN-001  /  v0.1.0')
    canvas.drawString(54, 30, 'Athota  |  Independent research  |  MIT license')
    canvas.drawRightString(A4[0]-54, 30, str(doc.page))
    canvas.restoreState()


lines = (ROOT / 'paper/manuscript.md').read_text().splitlines()
story = []
i = 0
while i < len(lines):
    line = lines[i].strip()
    if not line:
        i += 1
        continue
    if line == '<!-- pagebreak -->':
        story.append(PageBreak())
    elif line.startswith('### '):
        story.append(Paragraph(inline(line[4:]), styles['h3']))
    elif line.startswith('## '):
        story.append(Paragraph(inline(line[3:]), styles['h2']))
    elif line.startswith('# '):
        story.append(Paragraph(inline(line[2:]), styles['title']))
    elif line.startswith('$$') and line.endswith('$$'):
        path, width, height = equation(line[2:-2], 11)
        scale = min(1, WIDTH / width)
        story += [Spacer(1, 3), Image(str(path), width*scale, height*scale), Spacer(1, 12)]
    elif line.startswith('!['):
        match = re.match(r'!\[([^\]]*)\]\(([^)]+)\)', line)
        path = (ROOT / 'paper' / match.group(2)).resolve()
        with PILImage.open(path) as im:
            ratio = im.height / im.width
        story += [Spacer(1, 4), Image(str(path), WIDTH, WIDTH*ratio), Spacer(1, 7)]
    elif line.startswith('|'):
        rows = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            cells = [s.strip() for s in lines[i].strip().strip('|').split('|')]
            if not all(re.fullmatch(r':?-+:?', c) for c in cells):
                style = styles['headcell'] if not rows else styles['cell']
                rows.append([Paragraph(inline(c), style) for c in cells])
            i += 1
        table = Table(rows, colWidths=[WIDTH/len(rows[0])]*len(rows[0]), repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#14677d')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor('#f0f5f7'), colors.white]),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('LEFTPADDING', (0, 0), (-1, -1), 7), ('RIGHTPADDING', (0, 0), (-1, -1), 7),
            ('TOPPADDING', (0, 0), (-1, -1), 7), ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
        ]))
        story += [table, Spacer(1, 11)]
        continue
    else:
        paragraph = [lines[i]]
        while i+1 < len(lines) and lines[i+1].strip():
            if lines[i+1].startswith(('#', '$$', '|', '![', '<!--')):
                break
            i += 1
            paragraph.append(lines[i])
        text = '\n'.join(paragraph)
        text = text.replace('  \n', '<br/>').replace('\n', ' ')
        markup = inline(text).replace('&lt;br/&gt;', '<br/>')
        style = styles['caption'] if line.startswith('*Figure ') else styles['body']
        story.append(Paragraph(markup, style))
    i += 1

pdf = OUT / 'central-residence.pdf'
doc = SimpleDocTemplate(str(pdf), pagesize=A4, rightMargin=54, leftMargin=54,
                        topMargin=53, bottomMargin=49,
                        title='Central residence and density exposure of recoiling remnants in a Plummer potential',
                        author='Venkata Hemanth Athota', subject='Technical note BHRD-TN-001; version 0.1.0')
def invariant_canvas(*args, **kwargs):
    kwargs['invariant'] = 1
    return canvas.Canvas(*args, **kwargs)
doc.build(story, onFirstPage=page, onLaterPages=page, canvasmaker=invariant_canvas)
reader = PdfReader(pdf)
text = '\n'.join(p.extract_text() for p in reader.pages)
assert 'Venkata Hemanth Athota' in text
assert 'bharadwaj.hemanth.002@gmail.com' in text
assert len(reader.pages) == 6, f'Unexpected pagination: {len(reader.pages)}'
inputs = [ROOT / name for name in ['paper/manuscript.md', 'paper/build.py', 'paper/verify.py', 'paper/requirements.txt',
           'requirements-lock.txt', 'CITATION.cff', 'LICENSE', 'results/run_all.json']]
inputs += sorted((ROOT / 'results').glob('e*/*.csv'))
inputs += sorted((ROOT / 'results').glob('e*/*.png'))
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {'paper': 'BHRD-TN-001', 'version': '0.1.0', 'pages': len(reader.pages),
            'versions': {name: importlib.metadata.version(name) for name in ['reportlab', 'matplotlib', 'pillow', 'pypdf']},
            'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
            'output_sha256': {str(pdf.relative_to(ROOT)): sha(pdf)}}
(OUT / 'build.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(f'Built {pdf.relative_to(ROOT)} ({len(reader.pages)} pages); author and contact checked.')
