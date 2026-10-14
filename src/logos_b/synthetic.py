"""Generate invented profiles and observations, without reading source records."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .models import RemModel, WarModel

SEED = 20261012

def generate(out, n=4, seed=SEED):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    rem, war, observations = [], [], []
    for i in range(n):
        m = dict(ID=f'SYN_REM_{i+1:03d}', T_min=float(rng.uniform(8,16)),
                 rate_ug_min=float(rng.uniform(80,140)), k1_min=float(rng.uniform(.025,.06)),
                 k2_min=float(rng.uniform(.4,.9)), w1=float(rng.uniform(.001,.002)),
                 w2=float(rng.uniform(.006,.01)), ke0_min=float(rng.uniform(.15,.4)),
                 K_ng_ml=float(rng.uniform(5,12)), Hill=float(rng.uniform(1.2,2.2)),
                 E0=22., response_fraction=.75)
        rem.append(m)
        model = RemModel(m)
        t = np.linspace(0,80,49)
        cp, ce = model.values(t,np.ones(20),'cp'), model.values(t,np.ones(20),'ce')
        eeg = m['E0']*(1-m['response_fraction']*(ce/m['K_ng_ml'])**m['Hill']/(1+(ce/m['K_ng_ml'])**m['Hill']))
        for j, time in enumerate(t):
            observations.extend([dict(ID=m['ID'],time=float(time),kind='Cp',value=float(cp[j]*np.exp(rng.normal(0,.03))),unit='ng/mL'),
                                 dict(ID=m['ID'],time=float(time),kind='EEG',value=float(eeg[j]+rng.normal(0,.12)),unit='synthetic EEG unit')])
        m = dict(ID=f'SYN_WAR_{i+1:03d}',dose_mg=float(rng.uniform(40,70)),
                 k_h=float(rng.uniform(.012,.025)),scale_L_inv=float(rng.uniform(.09,.16)),
                 ka_h=1.2,lag_h=.4,kout_h=float(rng.uniform(.015,.035)),IC50_mg_L=float(rng.uniform(2,5)))
        war.append(m)
        model = WarModel(m)
        t = np.array([0,1,2,4,8,12,24,36,48,72,96,120,144.])
        cp, R = model.cp(t,[0],[1]), model.response(t)
        for j, time in enumerate(t):
            observations.extend([dict(ID=m['ID'],time=float(time),kind='Cp',value=float(cp[j]*np.exp(rng.normal(0,.025))),unit='mg/L'),
                                 dict(ID=m['ID'],time=float(time),kind='PCA',value=float(100*R[j]+(rng.normal(0,.25) if time else 0)),unit='percent')])
    (out/'profiles.json').write_text(json.dumps(dict(synthetic=True,seed=seed,Rem_models=rem,War_models=war),indent=2,allow_nan=False)+'\n')
    pd.DataFrame(observations).to_csv(out/'observations.csv',index=False,float_format='%.12g')
    rows=[]
    for i in range(12):
        baseline=float(rng.uniform(35,65))
        for course in range(1,5):
            value=baseline+30*np.log(course)+rng.normal(0,4)
            rows.append(dict(ID=f'SYN_MTX_{i+1:03d}',course=course,sampling_day=int(rng.integers(12,17)),MTX_PG2_5_nmol_L=float(value) if rng.random()>.08 else None))
    pd.DataFrame(rows).to_csv(out/'mtx_courses.csv',index=False,float_format='%.12g')
