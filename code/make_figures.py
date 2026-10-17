from pathlib import Path
import sys, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parents[1]; RES=ROOT/'result'; FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
sys.path.insert(0,str(Path(__file__).parent))
from dynamic_risk_surface_model import SCHEDULES, simulate

plt.rcParams.update({'font.size':10,'axes.titlesize':13,'axes.labelsize':10,'legend.fontsize':8})
def save(name):
    plt.tight_layout(); plt.savefig(FIG/name,dpi=300,bbox_inches='tight',facecolor='white'); plt.close()

# Figure 1 — snake-flow architecture
fig,ax=plt.subplots(figsize=(16,8.5)); ax.axis('off')
ax.text(.5,.965,'Dynamic Risk Surface Architecture for Bispecific-Antibody Toxicity',ha='center',fontsize=21,fontweight='bold')
ax.text(.5,.915,'Literature-parameterized computational proof-of-concept: no patient level calibration',ha='center',fontsize=13,style='italic')
boxes=[
(.04,.66,.20,.16,'1  Literature Anchors','Trials, reviews, real-world infection evidence'),
(.29,.66,.20,.16,'2  Dynamic Reserve Core','R, N, P, G: marrow and immune reserve'),
(.54,.66,.20,.16,'3  Treatment-Stress Layer','Normalized E(t), inflammation, microbial pressure'),
(.79,.66,.17,.16,'4  Hazard Layer','Cumulative infection hazard H(t)'),
(.72,.36,.24,.16,'5  Risk Surface','Frailty × dose intensity stability geometry'),
(.40,.36,.24,.16,'6  Latent Regimes','Controlled / strained / overloaded CTMC'),
(.08,.36,.24,.16,'7  Decision Layer','Supportive-care policies; hypothesis-generating')]
for x,y,w,h,t,s in boxes:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.012',facecolor='white',edgecolor='0.25',lw=1.7))
    ax.text(x+w/2,y+h*.66,t,ha='center',va='center',fontsize=12.3,fontweight='bold')
    ax.text(x+w/2,y+h*.34,s,ha='center',va='center',fontsize=9.7)
# snake arrows: top left->right, down, then bottom right->left
arrows=[((.24,.74),(.29,.74)),((.49,.74),(.54,.74)),((.74,.74),(.79,.74)),((.875,.66),(.84,.52)),((.72,.44),(.64,.44)),((.40,.44),(.32,.44))]
for a,b in arrows: ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='-|>',lw=1.7,connectionstyle='arc3,rad=0'))
ax.text(.5,.15,'Flow: evidence → biological state dynamics → treatment stress → cumulative hazard → risk geometry → latent regimes → decision hypotheses',ha='center',fontsize=10.5,style='italic')
save('Figure_1_Model_Architecture_REVISED.png')

# Figure 2 — patterned lines, markers, dose-event rug
D=pd.read_csv(RES/'Reference_Trajectory_Revised.csv'); fig,axs=plt.subplots(2,2,figsize=(14,9)); fig.suptitle('Reference Proof-of-Concept Dynamics Over 84 Days',fontsize=18,fontweight='bold')
styles=[('-', 'o'),('--','s'),('-.','^'),(':','D')]
ax=axs[0,0]
for (c,label),(ls,mk) in zip([('R','Marrow R'),('N','Neutrophil N'),('P','Platelet P'),('G','Immune G')],styles):
    ax.plot(D.day,D[c],ls=ls,marker=mk,markevery=70,ms=4,lw=1.7,label=label)
ax.axhline(.75,ls=(0,(5,3)),lw=1.1,label='illustrative recovery criterion'); ax.set_title('Reserve trajectories'); ax.set_xlabel('Day'); ax.legend(ncol=2)
ax=axs[0,1]
for (c,label),(ls,mk) in zip([('C','Inflammation C'),('B','Microbial pressure B'),('E','Normalized treatment-stress E(t)')],styles[:3]):
    ax.plot(D.day,D[c],ls=ls,marker=mk,markevery=70,ms=4,lw=1.7,label=label)
for d,a in SCHEDULES['reference_q14']: ax.axvline(d,alpha=.13,lw=.8,ls=':')
ax.set_title('Stress layer and structural exposure forcing'); ax.set_xlabel('Day'); ax.legend()
ax=axs[1,0]; ax.plot(D.day,D.stability,ls='-',marker='o',markevery=70,ms=4,lw=1.9,label='S(t)'); ax.axhline(.72,ls='--',label='computational boundary 0.72'); ax.axhline(.62,ls=':',label='computational boundary 0.62'); ax.set_title('Dynamic stability coordinate S(t)'); ax.set_xlabel('Day'); ax.legend()
ax=axs[1,1]; risk=1-np.exp(-D.H); ax.plot(D.day,risk,ls='-',marker='o',markevery=70,ms=4,lw=1.8,label='model-generated cumulative risk'); ax.plot(D.day,D.prob_controlled,ls='--',marker='s',markevery=70,ms=4,label='controlled'); ax.plot(D.day,D.prob_strained,ls='-.',marker='^',markevery=70,ms=4,label='strained'); ax.plot(D.day,D.prob_overloaded,ls=':',marker='D',markevery=70,ms=4,label='overloaded'); ax.set_title('Cumulative risk and latent regimes'); ax.set_xlabel('Day'); ax.legend(ncol=2)
save('Figure_2_Reference_Dynamics_REVISED.png')

# Figure 3 — surface + explicit contour interpretation
S=pd.read_csv(RES/'Dynamic_Risk_Surface_Revised.csv'); piv=S.pivot(index='frailty',columns='dose_intensity',values='infection_risk_output'); X,Y=np.meshgrid(piv.columns,piv.index); Z=piv.values
plt.figure(figsize=(10.5,7.5)); im=plt.pcolormesh(X,Y,Z,shading='auto',cmap='viridis'); levels=np.arange(.25,.601,.05); cs=plt.contour(X,Y,Z,levels=levels,colors='white',linewidths=1); plt.clabel(cs,fmt='%.2f'); plt.colorbar(im,label='Model-generated 84-day infection-risk output'); plt.xlabel('Relative dose intensity'); plt.ylabel('Frailty coefficient'); plt.title('Dynamic Infection-Vulnerability Surface',fontweight='bold',fontsize=16); plt.figtext(.5,.012,'Contour labels (e.g., 0.25, 0.30, 0.45) denote equal model-generated 84-day risk output; 0.25 = 25% under the model, not observed clinical incidence.',ha='center',fontsize=9.5); plt.subplots_adjust(bottom=.14); plt.savefig(FIG/'Figure_3_Dynamic_Risk_Surface_REVISED.png',dpi=300,bbox_inches='tight',facecolor='white'); plt.close()

# Figure 4 — scenarios with explicit NA meaning
sc=pd.read_csv(RES/'Scenario_Results_Revised.csv'); sc=sc[~sc.scenario.str.startswith('schedule_')].copy(); cols=['neutrophil_nadir','minimum_stability','recovery_day_N_ge_0_75','infection_risk_output','final_controlled_prob']; A=sc[cols].to_numpy(float); norm=(A-np.nanmin(A,axis=0))/(np.nanmax(A,axis=0)-np.nanmin(A,axis=0)+1e-12)
plt.figure(figsize=(13,7.5)); im=plt.imshow(norm,aspect='auto',cmap='viridis',alpha=.82); plt.colorbar(im,label='Column-normalized magnitude'); plt.xticks(range(5),['N nadir','Min stability','Recovery day','Infection risk','Controlled prob'],rotation=10); plt.yticks(range(len(sc)),[x.replace('_',' ').title() for x in sc.scenario]);
for i in range(len(sc)):
  for j in range(5):
    v=A[i,j]
    if np.isnan(v) and j==2:
        nad=A[i,0]; txt='NA*' if nad>=.75 else 'NA†'
    else: txt='NA' if np.isnan(v) else (f'{v:.1f}' if j==2 else f'{v:.3f}')
    plt.text(j,i,txt,ha='center',va='center',fontsize=9)
plt.title('Scenario Stress Test Across Mechanistic and Supportive-Care Perturbations',fontweight='bold',fontsize=14)
plt.figtext(.5,.012,'Recovery-day notation: NA* = N(t) never fell below 0.75, so recovery time is not applicable; NA† = N(t) fell below 0.75 but did not recover to 0.75 within 84 days. NA does not mean “no correlation”.',ha='center',fontsize=9.2)
plt.subplots_adjust(bottom=.15); plt.savefig(FIG/'Figure_4_Scenario_Stress_Test_REVISED.png',dpi=300,bbox_inches='tight',facecolor='white'); plt.close()

# Figure 5 — transparent sensitivity bars
pr=pd.read_csv(RES/'Global_Sensitivity_PRCC_Revised.csv'); q=pr[pr.outcome=='infection_risk_output'].sort_values('PRCC'); plt.figure(figsize=(10,7)); bars=plt.barh(q.parameter,q.PRCC,alpha=.42,edgecolor='black',linewidth=.9); plt.axvline(0,lw=1); plt.xlabel('PRCC with infection-risk output'); plt.title('Global Sensitivity of Infection Vulnerability',fontweight='bold',fontsize=16); plt.grid(axis='x',alpha=.15); save('Figure_5_Global_Sensitivity_PRCC_REVISED.png')

# Figure 6 — transparent histogram + gradient scatter
st=pd.read_csv(RES/'Stochastic_Perturbation_500_Revised.csv'); fig,axs=plt.subplots(1,2,figsize=(13,5)); axs[0].hist(st.infection_risk_output,bins=28,alpha=.38,edgecolor='black',linewidth=.7); axs[0].axvline(st.infection_risk_output.median(),ls='--',lw=1.4,label='median'); axs[0].legend(); axs[0].set_title('Stochastic risk distribution'); axs[0].set_xlabel('84-day model risk output'); axs[0].set_ylabel('Runs'); pts=axs[1].scatter(st.neutrophil_nadir,st.infection_risk_output,c=st.frailty,cmap='viridis',s=22,alpha=.58,edgecolors='none'); fig.colorbar(pts,ax=axs[1],label='Perturbed frailty coefficient'); axs[1].set_title('Cytopenia depth does not fully determine risk'); axs[1].set_xlabel('Neutrophil nadir'); axs[1].set_ylabel('84-day model risk output'); fig.suptitle('Stochastic Perturbation Ensemble (n=500)',fontweight='bold',fontsize=17); save('Figure_6_Stochastic_Uncertainty_REVISED.png')

# Figure 7 — transparent probability ribbons + lines, not solid stack blocks
fig,ax=plt.subplots(figsize=(11,6)); regime=[('prob_controlled','Controlled','-'),('prob_strained','Strained','--'),('prob_overloaded','Overloaded','-.')]
for c,label,ls in regime:
    ax.plot(D.day,D[c],ls=ls,lw=2,label=label); ax.fill_between(D.day,0,D[c],alpha=.08)
ax.set_ylim(0,1); ax.set_xlabel('Day'); ax.set_ylabel('Latent-state probability'); ax.set_title('Regime-Switching Toxicity State Probabilities',fontweight='bold',fontsize=16); ax.legend(ncol=3); ax.grid(alpha=.12); save('Figure_7_Regime_Switching_REVISED.png')

# Figure 8 — transparent bars with values
mc=pd.read_csv(RES/'Synthetic_Model_Comparison_AIC_BIC_Revised.csv').sort_values('delta_AIC'); plt.figure(figsize=(10,6)); bars=plt.bar(mc.model,mc.delta_AIC,alpha=.40,edgecolor='black',linewidth=.9); plt.axhline(2,ls='--',lw=1,label='ΔAIC = 2 reference');
for b,v in zip(bars,mc.delta_AIC): plt.text(b.get_x()+b.get_width()/2,b.get_height()+max(mc.delta_AIC.max()*.02,.2),f'{v:.2f}',ha='center',fontsize=9)
plt.ylabel('ΔAIC (relative to lowest AIC)'); plt.xticks(rotation=12,ha='right'); plt.title('Synthetic Structural Model Recovery',fontweight='bold',fontsize=16); plt.legend(); save('Figure_8_Model_Comparison_REVISED.png')

# Figure 9 — show solver robustness as deviations from reference, not two nearly-flat raw lines
rb=pd.read_csv(RES/'Numerical_Robustness_Revised.csv').copy(); ref=rb.loc[rb.setting=='reference'].iloc[0]
metrics9=[('infection_risk_output','Risk output'),('neutrophil_nadir','N nadir'),('minimum_stability','Min S')]
fig,axs=plt.subplots(1,2,figsize=(13,5.5)); x=np.arange(len(rb)); markers=['o','s','^']
for (col,label),mk in zip(metrics9,markers):
    dev=100*(rb[col]-ref[col])/ref[col]; axs[0].plot(x,dev,marker=mk,lw=1.7,label=label)
axs[0].axhline(0,ls='--',lw=1); axs[0].set_xticks(x,rb.setting); axs[0].set_ylabel('Relative deviation from reference (%)'); axs[0].set_title('Output deviation across solver settings'); axs[0].legend(); axs[0].grid(alpha=.15)
# second panel gives actual tolerances and absolute outputs
axs[1].axis('off'); rows=[]
for _,r in rb.iterrows(): rows.append([r.setting,f"{r.rtol:.0e}",f"{r.atol:.0e}",f"{r['dt']:.2f}",f"{r.infection_risk_output:.6f}",f"{r.neutrophil_nadir:.6f}",f"{r.minimum_stability:.6f}"])
tab=axs[1].table(cellText=rows,colLabels=['Setting','rtol','atol','dt','Risk','N nadir','Min S'],loc='center',cellLoc='center'); tab.auto_set_font_size(False); tab.set_fontsize(8.5); tab.scale(1.05,1.55); axs[1].set_title('Solver settings and absolute outputs',pad=14)
fig.suptitle('Numerical Robustness Across Solver Tolerances',fontweight='bold',fontsize=16); fig.text(.5,.015,'Near-zero deviations indicate that reported outputs are numerically stable to the tested tolerance/time-grid settings.',ha='center',fontsize=9.5); plt.subplots_adjust(bottom=.12); plt.savefig(FIG/'Figure_9_Numerical_Robustness_REVISED.png',dpi=300,bbox_inches='tight',facecolor='white'); plt.close()

# Figure 10 — large single-line workflow labels
fig,ax=plt.subplots(figsize=(16,7.5)); ax.axis('off'); items=[
('1  Published clinical anchors — trials, systematic reviews, real-world cohorts'),
('2  Literature-supported methods — ODEs, CTMC, LHS/PRCC, uncertainty and decision theory'),
('3  Provisional coefficients — kinetic rates, stress forcing and stability bands; not patient-fitted'),
('4  Model-generated outputs — trajectories, surfaces, scenario contrasts and uncertainty'),
('5  Future validation — registry calibration, competing risks, external validation and decision curves')]
for i,text in enumerate(items):
    y=.79-i*.145; ax.add_patch(FancyBboxPatch((.08,y),.84,.095,boxstyle='round,pad=.012,rounding_size=.012',facecolor='white',edgecolor='0.25',lw=1.7)); ax.text(.5,y+.0475,text,fontsize=11.7,fontweight='bold' if i in [0,4] else 'normal',ha='center',va='center')
    if i<4: ax.annotate('',xy=(.5,y-.035),xytext=(.5,y),arrowprops=dict(arrowstyle='-|>',lw=1.4))
ax.text(.5,.95,'Evidence and Parameter Provenance',ha='center',fontsize=21,fontweight='bold'); ax.text(.5,.055,'Integrity rule: literature anchors define plausibility; provisional coefficients generate hypotheses; clinical effect estimates require observed-data calibration.',ha='center',fontsize=10.5,style='italic'); save('Figure_10_Parameter_Provenance_REVISED.png')

# Figure 11 — patterned schedule lines
ss=pd.read_csv(RES/'Exposure_Schedule_Sensitivity.csv'); fig,axs=plt.subplots(1,2,figsize=(13,5.5)); style_map={'reference_q14':('-', 'o'),'weekly_after_stepup':('--','s'),'q21_after_stepup':('-.','^')}
for sch in SCHEDULES:
    t,Y,Ss,P,E=simulate(schedule=sch,dt=.1); ls,mk=style_map[sch]; axs[0].plot(t,E,ls=ls,marker=mk,markevery=70,ms=4,lw=1.7,label=sch); axs[1].plot(t,Y[:,1],ls=ls,marker=mk,markevery=70,ms=4,lw=1.7,label=sch)
axs[0].set_title('Normalized treatment-stress forcing'); axs[0].set_xlabel('Day'); axs[0].set_ylabel('E(t)'); axs[1].set_title('Schedule-sensitivity of neutrophil reserve'); axs[1].set_xlabel('Day'); axs[1].set_ylabel('N(t)'); axs[0].legend(); axs[1].legend(); axs[0].grid(alpha=.12); axs[1].grid(alpha=.12); fig.suptitle('Structural Exposure-Schedule Audit',fontweight='bold',fontsize=17); save('Figure_11_Exposure_Schedule_Audit.png')
print('Generated',len(list(FIG.glob('*.png'))),'figures in',FIG)
