from pathlib import Path
import os,uuid,argparse,json,hashlib
R=Path(__file__).parent
parser=argparse.ArgumentParser();parser.add_argument('--name',required=True)
parser.add_argument('--level',type=int,default=2);parser.add_argument('--L',type=float,default=64.)
parser.add_argument('--dt',type=float,default=4.);parser.add_argument('--duration',type=float,default=54000.)
parser.add_argument('--pressure-rtol',type=float,default=1e-7);parser.add_argument('--principal-only',action='store_true')
args=vars(parser.parse_args());name=args.pop('name');args['all_cases']=not args.pop('principal_only')
if Path(name).name!=name:raise ValueError('Use a simple new run name')
if (R/'runs'/name).exists():raise FileExistsError('Refusing existing output')
cache=R.parent/'scn'/uuid.uuid4().hex[:10];cache.mkdir(parents=True);os.environ['NUMBA_CACHE_DIR']=str(cache)
os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('NUMBA_NUM_THREADS','2')
from corrected_solver import simulate
s=simulate(R/'runs'/name,**args)
s['code_hashes']['run.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
s['threads']={'OPENBLAS_NUM_THREADS':os.environ['OPENBLAS_NUM_THREADS'],'NUMBA_NUM_THREADS':os.environ['NUMBA_NUM_THREADS']}
(R/'runs'/name/'summary.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
