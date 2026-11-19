"""Plot independently verified pressure maxima through startup and 15 h."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).parent;out=R/'verification'
fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
for tol,style in [('1e-07','o-'),('1e-08','x--')]:
    d=json.loads((out/f'pressure_level2_tol{tol}.json').read_text());assert d['passed']
    ax.semilogx([r['time_s'] for r in d['rows']],[r['max_pressure_Pa']/1000 for r in d['rows']],style,label='Relative tolerance '+tol)
ax.axhline(500,color='#8c2d04',linewidth=1.3,label='Prescribed inlet pressure')
ax.axhline(100,color='#555555',linewidth=1,linestyle=':',label='Ambient pressure')
ax.set(xlabel='Physical time (s)',ylabel='Maximum interior pressure (kPa)',ylim=(90,520),
       title='Full-system implicit pressure integration: finest mesh')
ax.grid(alpha=.2);ax.legend(loc='lower right');fig.savefig(out/'Pressure_bounds.png',dpi=240);fig.savefig(out/'Pressure_bounds.pdf')
plt.close(fig)
