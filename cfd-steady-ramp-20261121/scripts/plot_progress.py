"""Plot completed checkpoint data only, labelled as numerical iterations."""
from pathlib import Path
import json,csv,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tbm-ramp-plots')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/tbm-ramp-cache')
Path(os.environ['XDG_CACHE_HOME']).mkdir(parents=True,exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;O=R/'results';fig,axs=plt.subplots(2,2,figsize=(10,6))
for name,color in [('F01_upwind_long','#167ac6'),('F02_second_order_restore','#d76f14')]:
 D=O/name
 if not (D/'residuals.csv').exists():continue
 for ax,file,col,scale in [(axs[0,0],'residuals.csv','p',1),(axs[0,1],'residuals.csv','Ux',1),(axs[1,0],'boundary-history.csv','outlet_spread_percent',1),(axs[1,1],'internal-monitor.csv','low_fraction_0p5',100)]:
  q=list(csv.DictReader((D/file).open()));x=[float(r['iteration']) for r in q];y=[float(r[col])*scale for r in q]
  if file=='residuals.csv':ax.semilogy(x,y,c=color,label=name)
  else:ax.plot(x,y,c=color,label=name)
for ax,y in zip(axs.ravel(),['Initial p residual','Initial Ux residual','6-outlet spread / mean [%]','Monitor speed <0.5 m/s [%]']):ax.set(xlabel='Steady iteration (not physical time)',ylabel=y);ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Completed checkpoints only; ongoing numerical startup study',fontsize=11);fig.tight_layout();fig.savefig(R/'progress.png',dpi=150);plt.close(fig)
