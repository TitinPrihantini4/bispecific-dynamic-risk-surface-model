from pathlib import Path
import json, sys
import numpy as np, pandas as pd
from scipy.stats import qmc, rankdata
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(Path(__file__).parent))
from dynamic_risk_surface_model import BASE,simulate,simulate_fast,metrics,SCHEDULES
RES=ROOT/'result'; DATA=ROOT/'data'; RES.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)
rng=np.random.default_rng(20261017)

def run(p=None,dt=.1,schedule='reference_q14',rtol=1e-8,atol=1e-10):
    if dt>=.2 and rtol==1e-8 and atol==1e-10:
        t,Y,S,P,E=simulate_fast(p,dt,schedule)
    else:
        t,Y,S,P,E=simulate(p,dt,rtol,atol,schedule)
    return t,Y,S,P,E,metrics(t,Y,S,P)

def save_reference():
    t,Y,S,P,E,m=run(); df=pd.DataFrame(Y,columns=['R','N','P','G','C','B','H']); df.insert(0,'day',t); df['E']=E; df['stability']=S; df[['prob_controlled','prob_strained','prob_overloaded']]=P; df.to_csv(RES/'Reference_Trajectory_Revised.csv',index=False); return m

def scenarios():
    specs={'reference':({},'reference_q14'),'lower_dose_20pct':({'dose_intensity':.8},'reference_q14'),'gcsf_support':({'gcsf':.25},'reference_q14'),
    'antimicrobial_prophylaxis':({'prophylaxis':.06},'reference_q14'),'combined_adaptation':({'dose_intensity':.85,'gcsf':.25,'prophylaxis':.06},'reference_q14'),
    'high_frailty':({'frailty':.75},'reference_q14'),'high_marrow_injury':({'alpha_E':.095},'reference_q14'),'high_neutrophil_suppression':({'q_N':.135},'reference_q14'),
    'schedule_weekly':({},'weekly_after_stepup'),'schedule_q21':({},'q21_after_stepup')}
    rows=[]
    for name,(p,sch) in specs.items():
        *_,m=run(p,schedule=sch); rows.append({'scenario':name,'schedule':sch,**m})
    pd.DataFrame(rows).to_csv(RES/'Scenario_Results_Revised.csv',index=False)

def surface():
    rows=[]
    for f in np.linspace(.1,.9,17):
      for di in np.linspace(.7,1.2,21):
        *_,m=run({'frailty':float(f),'dose_intensity':float(di)},dt=.2); rows.append({'frailty':f,'dose_intensity':di,**m})
    pd.DataFrame(rows).to_csv(RES/'Dynamic_Risk_Surface_Revised.csv',index=False)

def prcc_for(X,y):
    Xr=np.column_stack([rankdata(X[:,j]) for j in range(X.shape[1])]); yr=rankdata(y); out=[]
    for j in range(X.shape[1]):
        others=[k for k in range(X.shape[1]) if k!=j]; A=np.column_stack([np.ones(len(Xr)),Xr[:,others]])
        rx=Xr[:,j]-A@np.linalg.lstsq(A,Xr[:,j],rcond=None)[0]; ry=yr-A@np.linalg.lstsq(A,yr,rcond=None)[0]
        out.append(np.corrcoef(rx,ry)[0,1])
    return out

def lhs():
    ranges={'alpha_E':(.045,.105),'q_N':(.06,.14),'q_G':(.02,.06),'lambda_B':(.010,.030),'h0':(.0008,.0022),'beta_N':(.9,2.2),'beta_G':(.8,2.0),'frailty':(.1,.9),'dose_intensity':(.75,1.20),'clear_B':(.07,.16)}
    names=list(ranges); sampler=qmc.LatinHypercube(d=len(names),seed=20261017); U=sampler.random(600); X=np.zeros_like(U)
    for j,n in enumerate(names): lo,hi=ranges[n]; X[:,j]=lo+(hi-lo)*U[:,j]
    rows=[]
    for vals in X:
        p=dict(zip(names,vals)); *_,m=run(p,dt=.25); rows.append({**p,**m})
    df=pd.DataFrame(rows); df.to_csv(RES/'LHS_600_Results_Revised.csv',index=False)
    pr=[]
    for outcome in ['infection_risk_output','minimum_stability','neutrophil_nadir']:
        vals=prcc_for(X,df[outcome].to_numpy())
        pr += [{'outcome':outcome,'parameter':n,'PRCC':v} for n,v in zip(names,vals)]
    pd.DataFrame(pr).to_csv(RES/'Global_Sensitivity_PRCC_Revised.csv',index=False)

def local_elasticity():
    *_,base=run(); rows=[]
    for n in ['alpha_E','q_N','q_G','lambda_B','h0','beta_N','beta_G','frailty','dose_intensity','clear_B']:
        x=BASE[n]; outs={}
        for sign in [-1,1]:
            *_,m=run({n:x*(1+.01*sign)},dt=.15); outs[sign]=m
        for o in ['infection_risk_output','minimum_stability','neutrophil_nadir']:
            elast=(outs[1][o]-outs[-1][o])/(.02*x)*(x/base[o])
            rows.append({'outcome':o,'parameter':n,'elasticity':elast})
    pd.DataFrame(rows).to_csv(RES/'Local_Elasticity_Revised.csv',index=False)

def stochastic():
    rows=[]
    sds={'alpha_E':.10,'q_N':.10,'q_G':.10,'lambda_B':.10,'h0':.10,'frailty':.12,'dose_intensity':.08}
    for i in range(500):
        p={}
        for n,sd in sds.items(): p[n]=max(1e-8,BASE[n]*(1+rng.normal(0,sd)))
        p['frailty']=min(max(p['frailty'],.02),.98); p['dose_intensity']=min(max(p['dose_intensity'],.55),1.35)
        *_,m=run(p,dt=.5); rows.append({'run':i+1,**p,**m})
    pd.DataFrame(rows).to_csv(RES/'Stochastic_Perturbation_500_Revised.csv',index=False)

def robustness():
    rows=[]
    for label,rt,at,dt in [('coarse',1e-6,1e-8,.2),('reference',1e-8,1e-10,.1),('tight',1e-10,1e-12,.05)]:
        *_,m=run(dt=dt,rtol=rt,atol=at); rows.append({'setting':label,'rtol':rt,'atol':at,'dt':dt,**m})
    pd.DataFrame(rows).to_csv(RES/'Numerical_Robustness_Revised.csv',index=False)

def model_comparison():
    t,Y,*_=simulate(dt=.1); obs_t=np.arange(0,85,7); true=np.interp(obs_t,t,Y[:,1]); y=true+rng.normal(0,.018,len(obs_t)); n=len(y)
    def score(pred,k,name):
        rss=float(np.sum((y-pred)**2)); aic=n*np.log(rss/n)+2*k; bic=n*np.log(rss/n)+k*np.log(n); return {'model':name,'n_parameters':k,'RSS':rss,'RMSE':np.sqrt(rss/n),'AIC':aic,'BIC':bic}
    # dynamic reserve: estimate multiplicative alpha_E and q_N around mechanistic model
    def resid(z):
        tt,YY,*_=simulate({'alpha_E':BASE['alpha_E']*np.exp(z[0]),'q_N':BASE['q_N']*np.exp(z[1])},dt=.2); return np.interp(obs_t,tt,YY[:,1])-y
    z=least_squares(resid,[0,0]).x; tt,YY,*_=simulate({'alpha_E':BASE['alpha_E']*np.exp(z[0]),'q_N':BASE['q_N']*np.exp(z[1])},dt=.2); pred=np.interp(obs_t,tt,YY[:,1]); rows=[score(pred,2,'Dynamic reserve ODE')]
    X=np.column_stack([np.ones(n),np.sin(2*np.pi*obs_t/28),np.cos(2*np.pi*obs_t/28)]); pred=X@np.linalg.lstsq(X,y,rcond=None)[0]; rows.append(score(pred,3,'Harmonic 28-day curve'))
    def biex(z): return z[0]+z[1]*np.exp(-obs_t/max(z[2],1e-3))+z[3]*np.exp(-obs_t/max(z[4],1e-3))
    zz=least_squares(lambda z:biex(z)-y,[.75,.1,10,.1,40],bounds=([0,-1,.1,-1,.1],[2,1,200,1,200])).x; rows.append(score(biex(zz),5,'Bi-exponential phenomenological'))
    d=pd.DataFrame(rows); d['delta_AIC']=d.AIC-d.AIC.min(); d.to_csv(RES/'Synthetic_Model_Comparison_AIC_BIC_Revised.csv',index=False)

def schedule_sensitivity():
    rows=[]
    for sch in SCHEDULES:
        *_,m=run(schedule=sch); rows.append({'schedule':sch,**m})
    pd.DataFrame(rows).to_csv(RES/'Exposure_Schedule_Sensitivity.csv',index=False)

def main():
    m=save_reference(); scenarios(); surface(); lhs(); local_elasticity(); stochastic(); robustness(); model_comparison(); schedule_sensitivity()
    with open(RES/'Analysis_Summary_Revised.json','w') as f: json.dump(m,f,indent=2)
    print(json.dumps(m,indent=2))
if __name__=='__main__': main()
