"""Generate proposed replacement figures from corrected outputs only."""
from pathlib import Path
import json,hashlib,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).parent
assessment=json.loads((ROOT/'assessment.json').read_text())
if not assessment['accepted']:raise RuntimeError('Convergence gate has not passed')
OUT=ROOT/'figures';OUT.mkdir(exist_ok=False)
BASE=ROOT/'runs'/assessment['selected_run']
p=json.loads((BASE/'summary.json').read_text())['parameters']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':14,
                     'axes.titlesize':14,'legend.fontsize':11,'pdf.fonttype':42,
                     'savefig.dpi':300,'axes.spines.top':False,'axes.spines.right':False})
sets={name:np.load(BASE/(name+'_snapshots.npz')) for name in ['constant','variable','exp_0.001','exp_0.0005']}
flow=np.load(BASE/'flow_snapshots.npz')
x=flow['x'];y=flow['y'];captions=[];provenance=[];used_sources={};source_hashes={}
labels={'c':r'$c/c_0$','ca':r'$10^3 c_a$','cp':r'$10^3 c_p$'}
scales={'c':1/p['c_inlet'],'ca':1000,'cp':1000}
def get(case,hours,f):
    d=flow if case=='flow' else sets[case]
    key=f'{int(round(hours*3600/p["dt"]))}_{f}'
    name='flow_snapshots.npz' if case=='flow' else case+'_snapshots.npz'
    used_sources.setdefault(name,set()).add(key)
    return d[key]
def save(fig,num,caption):
    fig.savefig(OUT/f'Figure_{num:02d}.png',bbox_inches='tight')
    fig.savefig(OUT/f'Figure_{num:02d}.pdf',bbox_inches='tight')
    sources={}
    for name,keys in used_sources.items():
        if name not in source_hashes:source_hashes[name]=hashlib.sha256((BASE/name).read_bytes()).hexdigest()
        sources[name]={'sha256':source_hashes[name],'array_keys_or_columns':sorted(keys)}
    provenance.append({'figure':num,'selected_run':assessment['selected_run'],'sources':sources,
        'caption':caption,'output_sha256':{ext:hashlib.sha256((OUT/f'Figure_{num:02d}.{ext}').read_bytes()).hexdigest() for ext in ('png','pdf')}})
    used_sources.clear()
    captions.append((num,caption));plt.close(fig)
def surf(ax,z,title,zlabel,xmax=2.,color=None,clim=None):
    ii=x<=xmax+1e-10
    X,Y=np.meshgrid(x[ii],y,indexing='ij')
    if color:ax.plot_wireframe(X,Y,z[ii],rstride=max(1,len(X)//22),cstride=max(1,len(y)//40),color=color,linewidth=.5)
    else:
        limits={} if clim is None else dict(vmin=clim[0],vmax=clim[1])
        ax.plot_surface(X,Y,z[ii],cmap='viridis',linewidth=0,antialiased=True,**limits)
    ax.set(xlabel='x (m)',ylabel='y (m)',zlabel=zlabel,title=title)
    ax.xaxis.label.set_fontsize(15);ax.yaxis.label.set_fontsize(15);ax.zaxis.label.set_fontsize(15)
    ax.xaxis.labelpad=18;ax.yaxis.labelpad=18;ax.zaxis.labelpad=14
    ax.view_init(elev=24,azim=-125);ax.set_box_aspect((1.5,1,.8))
    ax.tick_params(labelsize=14)

def figure3d(figsize):
    fig=plt.figure(figsize=figsize)
    # 3D label extents are incompletely handled by constrained_layout.
    # Reserve explicit gutters so one panel's z label cannot cover its neighbour.
    fig.subplots_adjust(left=.07,right=.99,bottom=.18,top=.91,wspace=.25)
    return fig
def profiles(axleft,axright,z,label,title='',panel_labels=('a','b')):
    for yloc in [0,.2,.4,.6,.8]:
        profile=np.array([np.interp(yloc,y,row) for row in z])
        axleft.plot(x,profile,label=f'y = {yloc:g} m')
    for xloc in [0,.2,.5,1,2,4]:
        profile=np.array([np.interp(xloc,x,z[:,j]) for j in range(len(y))])
        axright.plot(y,profile,label=f'x = {xloc:g} m')
    for ax,coord,panel in [(axleft,'x',panel_labels[0]),(axright,'y',panel_labels[1])]:
        ax.set(xlabel=f'{coord} (m)',ylabel=label,title=f'({panel}) {title}'.strip())
        ax.grid(alpha=.18);ax.legend(loc='best',fontsize=10)
    axleft.set_xlim(0,4);axright.set_xlim(0,p['H'])
# Figures 2-5: hydraulic fields.
fig=figure3d((10,4.8))
for k,h in enumerate([.5,5]):
    ax=fig.add_subplot(1,2,k+1,projection='3d')
    surf(ax,get('flow',h,'p')/1000,f'({chr(97+k)}) t = {h:g} h','Pressure (kPa)',xmax=4,clim=(100,500))
    ax.set_zlim(100,500)
save(fig,2,'Pressure surfaces at 0.5 h (a) and 5 h (b), with ambient downstream pressure 100 kPa and inlet pressure 500 kPa.')
fig,axs=plt.subplots(2,1,figsize=(9,5),layout='constrained')
for ax,h,letter in zip(axs,[.5,5],['a','b']):
    inside=x<=4+1e-10
    cs=ax.contourf(x[inside],y,get('flow',h,'p')[inside].T/1000,levels=np.linspace(100,500,21),cmap='viridis')
    ax.set(xlabel='x (m)',ylabel='y (m)',title=f'({letter}) t = {h:g} h')
fig.colorbar(cs,ax=axs,label='Pressure (kPa)',shrink=.9)
save(fig,3,'Pressure contours at 0.5 h (a) and 5 h (b). The displayed interval is 0-4 m; the downstream numerical boundary is at 64 m.')
for n,h in [(4,.5),(5,5)]:
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    profiles(*axs,get('flow',h,'p')/1000,'Pressure (kPa)',f't = {h:g} h')
    for ax in axs:ax.set_ylim(100,500)
    axs[1].legend(loc='lower left',ncol=3,fontsize=9)
    save(fig,n,f'Pressure cross sections at t = {h:g} h: longitudinal profiles at selected y coordinates (a) and transverse profiles at selected x coordinates (b).')
# Three surface panels.
for num,case,h in [(6,'constant',15),(8,'exp_0.001',15)]:
    fig=figure3d((13,4.8))
    for j,f in enumerate(['c','ca','cp']):
        ax=fig.add_subplot(1,3,j+1,projection='3d')
        surf(ax,get(case,h,f)*scales[f],f'({chr(97+j)}) t = {h} h',labels[f])
        ax.set_zlim(bottom=0)
    save(fig,num,f'Concentration surfaces at 15 h for {"constant inlet and constant dispersion" if case=="constant" else "exponential inlet with phi = 0.001 s^-1 and constant dispersion"}. The displayed x interval is 0-2 m; the computation extends to 64 m.')
for num,case in [(7,'constant'),(16,'variable')]:
    fig,axs=plt.subplots(3,2,figsize=(10,10),layout='constrained')
    for j,f in enumerate(['c','ca','cp']):
        profiles(*axs[j],get(case,15,f)*scales[f],labels[f],panel_labels=(chr(97+2*j),chr(98+2*j)))
        if j>0:
            for ax in axs[j]:ax.get_legend().remove()
    save(fig,num,f'Cross sections at 15 h for {case} dispersion and constant inlet concentration. Dissolved concentration c (a,b), adsorbed concentration ca (c,d), and adsorbed concentration cp (e,f). Left panels show longitudinal profiles; right panels show transverse profiles. The coordinate legends in (a) and (b) apply to all rows.')
for num,h,f in [(9,15,'c'),(10,2,'c'),(11,2,'ca'),(12,2,'cp')]:
    fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for j,case in enumerate(['exp_0.001','exp_0.0005']):
        profiles(*axs[j],get(case,h,f)*scales[f],labels[f],r'$\varphi = '+case[4:]+r'$ s$^{-1}$',panel_labels=(chr(97+2*j),chr(98+2*j)))
    ceiling=max(float(get(case,h,f).max())*scales[f] for case in ['exp_0.001','exp_0.0005'])
    for ax in axs.flat:ax.set_ylim(0,1.05*ceiling)
    save(fig,num,f'Cross sections of {f} at {h} h under exponential inlet histories with phi = 0.001 s^-1 (a,b) and 0.0005 s^-1 (c,d), using constant dispersion. Left panels show longitudinal profiles; right panels show transverse profiles.')
fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
for case,label in [('constant','Constant inlet'),('exp_0.001',r'$\varphi=0.001$ s$^{-1}$'),('exp_0.0005',r'$\varphi=0.0005$ s$^{-1}$')]:
    h=np.loadtxt(BASE/(case+'_history.csv'),delimiter=',',skiprows=1)
    used_sources.setdefault(case+'_history.csv',set()).update(['time_s','c_at_005_mid','ca_at_005_mid','cp_at_005_mid'])
    for ax,col,f in zip(axs,[7,8,9],['c','ca','cp']):
        ax.plot(h[:,0]/3600,h[:,col]*scales[f],label=label)
        ax.set(xlabel='Time (h)',ylabel=labels[f]);ax.grid(alpha=.18)
for j,ax in enumerate(axs):ax.set_title(f'({chr(97+j)})');ax.legend(fontsize=12)
save(fig,13,'Time evolution at (x,y) = (0.05,0.4) m for constant and exponential inlet histories. All histories are sampled at this same physical point.')
# Variable D at 2 and 15 h; and comparison of dispersion choices at 15 h.
for num,overlays in [(14,[('variable',2,'#238b45','2 h'),('variable',15,'#225ea8','15 h')]),
                     (15,[('constant',15,'#225ea8','Constant D'),('variable',15,'#d95f0e','Variable D')])]:
    fig=figure3d((13,4.8))
    for j,f in enumerate(['c','ca','cp']):
        ax=fig.add_subplot(1,3,j+1,projection='3d')
        for case,h,color,label in overlays:surf(ax,get(case,h,f)*scales[f],f'({chr(97+j)})',labels[f],color=color)
        ax.set_zlim(bottom=0)
        ax.legend([Line2D([0],[0],color=o[2]) for o in overlays],[o[3] for o in overlays],loc='upper left',fontsize=13)
    save(fig,num,('Variable-dispersion concentration surfaces under constant inlet concentration at 2 h (green) and 15 h (blue).' if num==14 else 'Comparison of constant anisotropic dispersion and spatially variable scalar dispersion under constant inlet concentration at 15 h. Blue: constant dispersion. Orange: variable dispersion.'))
# Additional verification visuals; not assigned manuscript figure numbers.
fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
for row,case in enumerate(['constant','variable']):
    h=np.loadtxt(BASE/(case+'_history.csv'),delimiter=',',skiprows=1);t=h[:,0]/3600
    ax=axs[row,0]
    for col,label in [(1,'Stored solute'),(2,'Cumulative input'),(4,'Cumulative degradation')]:
        ax.plot(t,h[:,col],label=label)
    ax.set(xlabel='Time (h)',ylabel='Integrated model amount',title=case.capitalize()+' dispersion')
    ax.legend();ax.grid(alpha=.18)
    ax=axs[row,1]
    ax.semilogy(t[1:],np.maximum(np.abs(h[1:,5])/np.maximum(h[1:,2],1e-300),1e-17))
    ax.set(xlabel='Time (h)',ylabel=r'$|R(t)|/I(t)$')
    ax.grid(alpha=.2)
fig.savefig(OUT/'Mass_balance.png',bbox_inches='tight');fig.savefig(OUT/'Mass_balance.pdf',bbox_inches='tight');plt.close(fig)
(OUT/'Mass_balance_caption.txt').write_text('Mass-balance diagnostics for constant and variable dispersion with constant inlet concentration. Left: total stored solute, cumulative boundary input, and cumulative degradation. Right: absolute cumulative balance residual R(t), divided by cumulative input I(t) at the same physical time. R=M-M0-I+O+G, with all phases included in M and actual numerical boundary fluxes used for I and O. Outflow is retained in the ledger; its values are negligible on the plotted scale. A floor of 1e-17 is used only for logarithmic display of vanishing residual ratios. Integrated amounts retain the supplied model concentration basis per unit transverse thickness.\n',encoding='utf-8')
(OUT/'Figure_captions.txt').write_text('\n\n'.join(f'Figure {n}. {text}' for n,text in sorted(captions))+'\n\nAll figures use the fixed 0.04 m inlet, the 64 m numerical domain, and the convergence-tested inlet-focused mesh. Pressure is integrated by an unsplit adaptive implicit BDF/NDF solver. A smooth geometric transport-step ramp is followed by a 2 s step after 800 s. Pressure is absolute pressure; c0=0.01, and the sorbed concentrations are multiplied by 1000 for display without changing their supplied concentration basis. No additional mass unit is assumed. Boundary values are algebraic and are excluded from storage integrals. See Verification_report.txt for measured spatial, temporal and domain changes and their limits.\n',encoding='utf-8')
with (OUT/'Figure_captions.txt').open('a',encoding='utf-8') as caption_file:
    caption_file.write('\nProfiles at symmetric y coordinates overlap: y=0 and 0.8 m, and y=0.2 and 0.6 m. All specified profiles are included.\n')
print(f'Created {len(captions)} convergence-accepted numerical figures and a mass-balance diagnostic.')
(OUT/'figure_data_manifest.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
gallery=''.join(f'<section id="figure-{n}"><h2>Figure {n}</h2><a href="Figure_{n:02d}.png"><img loading="lazy" src="Figure_{n:02d}.png" alt="{html.escape(caption,quote=True)}"></a><p>{html.escape(caption)}</p><p><a href="Figure_{n:02d}.pdf">Vector PDF</a> | <a href="Figure_{n:02d}.png">300 dpi PNG</a></p></section>' for n,caption in sorted(captions))
nav=' '.join(f'<a href="#figure-{n}">{n}</a>' for n in range(2,17))
(OUT/'Gallery.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Verified Figures 2-16</title><style>body{max-width:1100px;margin:32px auto;padding:0 20px;font:16px/1.6 system-ui;color:#17212b}img{width:100%;height:auto}section{border-top:1px solid #dbe2e8;padding:20px 0}a{color:#155c87}nav a{display:inline-block;padding:4px 10px}h1{font-size:28px}</style><body><h1>Verified Figures 2-16</h1><p>Replacement scientific figures from the accepted 15-hour calculations. Click an image to view the full-resolution PNG. Captions and verification details accompany this collection.</p><nav>'+nav+'</nav>'+gallery+'</body></html>',encoding='utf-8')
for z in [flow,*sets.values()]:z.close()

