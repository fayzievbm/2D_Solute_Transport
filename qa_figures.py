"""Check export completeness and make contact sheets for visual inspection."""
from pathlib import Path
import json,hashlib,subprocess,sys,shutil,os
from PIL import Image,ImageOps,ImageDraw
from pypdf import PdfReader
R=Path(__file__).parent;F=R/'figures';Q=R/'figure_qa';Q.mkdir(exist_ok=True)
rendered=Q/'pdf_render';rendered.mkdir(exist_ok=True)
bundled=Path(sys.executable).parents[1]/'native'/'poppler'/'Library'/'bin'/'pdftoppm.exe'
poppler=os.environ.get('PDFTOPPM') or (str(bundled) if bundled.exists() else shutil.which('pdftoppm'))
assert poppler,'A PDF renderer is required for visual QA'
def render(pdf,name):
    subprocess.run([poppler,'-r','140','-png','-singlefile',str(pdf),str(rendered/name)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    return rendered/(name+'.png')
records=json.loads((F/'figure_data_manifest.json').read_text())
assert sorted(r['figure'] for r in records)==list(range(2,17))
rows=[]
for r in records:
    n=r['figure'];png=F/f'Figure_{n:02d}.png';pdf=F/f'Figure_{n:02d}.pdf'
    for ext in ['png','pdf']:assert hashlib.sha256((F/f'Figure_{n:02d}.{ext}').read_bytes()).hexdigest()==r['output_sha256'][ext]
    with Image.open(png) as im:
        assert min(im.size)>=1000
        dpi=im.info.get('dpi');assert dpi and min(dpi)>299
        dimensions=im.size
    reader=PdfReader(pdf);assert len(reader.pages)==1
    words=reader.pages[0].extract_text()
    assert words and ('Time' in words or 'x (m)' in words or 'Pressure' in words)
    render(pdf,f'Figure_{n:02d}')
    rows.append({'figure':n,'png_pixels':dimensions,'dpi':dpi,'pdf_pages':len(reader.pages),'export_checks_passed':True})
for batch,start in enumerate(range(2,17,5),1):
    sheet=Image.new('RGB',(1800,1800),'#e9edf0');draw=ImageDraw.Draw(sheet)
    for k,n in enumerate(range(start,min(start+5,17))):
        with Image.open(rendered/f'Figure_{n:02d}.png') as im:
            thumb=ImageOps.contain(im.convert('RGB'),(870,540))
        col=k%2;row=k//2;x=col*900+(900-thumb.width)//2;y=row*600+40
        sheet.paste(thumb,(x,y));draw.text((col*900+20,row*600+12),f'Figure {n}',fill='black')
    sheet.save(Q/f'contact_sheet_{batch}.png')
extra=Image.new('RGB',(1800,1600),'white')
for i,(source,name) in enumerate([(F/'Mass_balance.pdf','Mass_balance'),(R/'verification'/'Pressure_bounds.pdf','Pressure_bounds')]):
    path=render(source,name)
    with Image.open(path) as im:thumb=ImageOps.contain(im.convert('RGB'),(1760,760))
    extra.paste(thumb,((1800-thumb.width)//2,i*800+(800-thumb.height)//2))
extra.save(Q/'diagnostics_contact_sheet.png')
(Q/'export_audit.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
print('All 15 numbered PNG/PDF pairs verified and all PDFs rendered. Inspect PDF contact sheets and full-sized exports before delivery.')
