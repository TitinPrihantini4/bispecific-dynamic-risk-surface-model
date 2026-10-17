from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
HORIZON=84.0
INITIAL=np.array([0.92,0.95,0.94,0.92,0.08,0.07,0.0],float)
BASE={"rho_R":0.035,"alpha_E":0.070,"alpha_C":0.020,"k_N":0.085,"d_N":0.080,"q_N":0.095,
"k_P":0.055,"d_P":0.052,"q_P":0.040,"k_G":0.030,"d_G":0.028,"q_G":0.035,"c_E":0.050,
"c_B":0.030,"d_C":0.180,"lambda_B":0.018,"clear_B":0.110,"K_clear":0.40,"h0":0.00135,
"beta_N":1.55,"beta_G":1.35,"beta_B":1.10,"beta_C":0.75,"beta_F":0.90,"frailty":0.35,
"dose_intensity":1.0,"gcsf":0.0,"prophylaxis":0.0,"stress_decay":0.35}

# This is deliberately NOT a drug-specific PK regimen. It is a normalized treatment-stress forcing.
# Three step-up stress pulses are followed by repeated full-intensity pulses every 14 days.
# Alternative weekly and 21-day schedules are tested as structural schedule-sensitivity scenarios.
SCHEDULES={
    "reference_q14":[(0,0.20),(7,0.35),(14,0.50),(28,0.55),(42,0.55),(56,0.55),(70,0.55)],
    "weekly_after_stepup":[(0,0.20),(7,0.35),(14,0.50)]+[(d,0.55) for d in range(21,85,7)],
    "q21_after_stepup":[(0,0.20),(7,0.35),(14,0.50),(35,0.55),(56,0.55),(77,0.55)]
}

def exposure(t,p,schedule="reference_q14"):
    decay=p.get("stress_decay",0.18); di=p.get("dose_intensity",1.0)
    return sum(a*di*np.exp(-decay*(t-d)) for d,a in SCHEDULES[schedule] if t>=d)

def rhs(t,y,p,schedule="reference_q14"):
    R,N,P,G,C,B,H=y; E=exposure(t,p,schedule); F=p["frailty"]; ug=p["gcsf"]; up=p["prophylaxis"]
    return np.array([
        p["rho_R"]*(1-R)-p["alpha_E"]*E*R-p["alpha_C"]*C*R,
        p["k_N"]*R*(1+ug)-p["d_N"]*N-p["q_N"]*E*N,
        p["k_P"]*R-p["d_P"]*P-p["q_P"]*E*P,
        p["k_G"]*R-p["d_G"]*G-p["q_G"]*E*G,
        p["c_E"]*E+p["c_B"]*B-p["d_C"]*C,
        p["lambda_B"]*(1+F)-(p["clear_B"]*(N*G)/(p["K_clear"]+N*G)+up)*B,
        p["h0"]*np.exp(p["beta_N"]*max(0,1-N)+p["beta_G"]*max(0,1-G)+p["beta_B"]*B+p["beta_C"]*C+p["beta_F"]*F)
    ])

def stability(Y):
    R,N,P,G,C,B=Y[:,:6].T
    return 0.40*R+0.25*N+0.15*P+0.20*G-0.12*C-0.10*B

def regime_probabilities(t,S):
    p=np.zeros((len(t),3)); p[0]=[0.88,0.10,0.02]
    for i in range(1,len(t)):
        s=S[i-1]; dt=t[i]-t[i-1]
        a=0.010+0.14/(1+np.exp(25*(s-0.72)))
        b=0.055/(1+np.exp(-22*(s-0.68)))+0.008
        c=0.008+0.11/(1+np.exp(28*(s-0.62)))
        d=0.035/(1+np.exp(-24*(s-0.64)))+0.004
        Q=np.array([[-a,a,0],[b,-b-c,c],[0,d,-d]])
        p[i]=p[i-1]+dt*(p[i-1]@Q)
        p[i]=np.clip(p[i],0,None); p[i]/=p[i].sum()
    return p

def simulate(p=None,dt=.1,rtol=1e-8,atol=1e-10,schedule="reference_q14"):
    pp=BASE.copy(); pp.update({} if p is None else p)
    t=np.arange(0,HORIZON+1e-9,dt)
    sol=solve_ivp(lambda tt,yy:rhs(tt,yy,pp,schedule),(0,HORIZON),INITIAL,t_eval=t,rtol=rtol,atol=atol,method="RK45")
    if not sol.success: raise RuntimeError(sol.message)
    Y=sol.y.T; S=stability(Y); probs=regime_probabilities(t,S)
    E=np.array([exposure(x,pp,schedule) for x in t])
    return t,Y,S,probs,E

def metrics(t,Y,S,probs):
    N=Y[:,1]; H=Y[:,6]
    i=int(np.argmin(N)); j=int(np.argmin(S))
    after=np.where((t>=t[i])&(N>=0.75))[0]
    recovery=(float(t[after[0]]) if len(after) else np.nan) if N[i] < 0.75 else np.nan
    return {"neutrophil_nadir":float(N[i]),"neutrophil_nadir_day":float(t[i]),
            "marrow_nadir":float(Y[:,0].min()),"platelet_nadir":float(Y[:,2].min()),"immune_nadir":float(Y[:,3].min()),
            "minimum_stability":float(S[j]),"minimum_stability_day":float(t[j]),"recovery_day_N_ge_0_75":recovery,
            "infection_risk_output":float(1-np.exp(-H[-1])),"final_controlled_prob":float(probs[-1,0]),
            "final_strained_prob":float(probs[-1,1]),"final_overloaded_prob":float(probs[-1,2]),
            "days_below_0_72":float(np.trapezoid((S<0.72).astype(float),t)),"days_below_0_62":float(np.trapezoid((S<0.62).astype(float),t))}


def simulate_fast(p=None,dt=.25,schedule="reference_q14"):
    pp=BASE.copy(); pp.update({} if p is None else p)
    t=np.arange(0,HORIZON+1e-9,dt); Y=np.empty((len(t),7)); Y[0]=INITIAL
    for i in range(1,len(t)):
        ti=t[i-1]; y=Y[i-1]
        k1=rhs(ti,y,pp,schedule); k2=rhs(ti+dt/2,y+dt*k1/2,pp,schedule); k3=rhs(ti+dt/2,y+dt*k2/2,pp,schedule); k4=rhs(ti+dt,y+dt*k3,pp,schedule)
        Y[i]=y+dt*(k1+2*k2+2*k3+k4)/6
    S=stability(Y); probs=regime_probabilities(t,S); E=np.array([exposure(x,pp,schedule) for x in t])
    return t,Y,S,probs,E
