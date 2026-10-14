"""Fresh generation, sequential estimation, matched optimization and verification."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from . import __version__
from .synthetic import generate, SEED
from .models import RemModel, WarModel
from .fitting import fit_rem_pk, fit_rem_pd, refit_war_baseline
from .continuation import solve_floor
from .benchmark import canonical, run_frontier


def save(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def run(out, n=4, seed=SEED):
    out=Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory must be empty to prevent stale results.')
    out.mkdir(parents=True,exist_ok=True)
    inputs=out/'data'; generate(inputs,n,seed)
    data=json.loads((inputs/'profiles.json').read_text())
    obs=pd.read_csv(inputs/'observations.csv')
    tables=out/'tables'; tables.mkdir()
    figures=out/'figures'; figures.mkdir()
    models={'Remifentanil':[], 'Warfarin':[]}; fitting=[]; checks=[]; results=[]; floor_rows=[]; controls=[]
    def check(name,passed,detail=None):
        checks.append(dict(check=name,passed=bool(passed),detail=detail))
    for m in data['Rem_models']:
        d=obs[obs.ID==m['ID']]; pk=d[d.kind=='Cp']; pdobs=d[d.kind=='EEG']
        w,k,pkfit=fit_rem_pk(pk.time.to_numpy(),pk.value.to_numpy(),m['rate_ug_min'],m['T_min'])
        pdfit=fit_rem_pd(pdobs.time.to_numpy(),pdobs.value.to_numpy(),w,k,m['rate_ug_min'],m['T_min'])
        ke,K,h=np.exp(pdfit.x[:3])
        fitted=dict(m,w1=float(w[0]),w2=float(w[1]),k1_min=float(k[0]),k2_min=float(k[1]),ke0_min=float(ke),K_ng_ml=float(K),Hill=float(h),E0=float(pdfit.x[3]),response_fraction=float(pdfit.x[4]))
        models['Remifentanil'].append(RemModel(fitted))
        fitting.append(dict(ID=m['ID'],drug='Remifentanil',PK_success=bool(pkfit.success),PD_success=bool(pdfit.success),PK_RMSE=float(np.sqrt(np.mean(pkfit.fun**2))),PD_RMSE=float(np.sqrt(np.mean(pdfit.fun**2)))))
    for m in data['War_models']:
        d=obs[obs.ID==m['ID']]; pk=d[(d.kind=='Cp')&(d.time>0)]
        def residual(z):
            temp=dict(m,k_h=float(np.exp(z[0])),scale_L_inv=float(np.exp(z[1])))
            cp=WarModel(temp).cp(pk.time.to_numpy(),[0],[1])
            return np.log(cp+.1)-np.log(pk.value.to_numpy()+.1)
        pkfit=least_squares(residual,np.log([.02,.12]),bounds=(np.log([.001,.001]),np.log([.1,2])),ftol=1e-10,xtol=1e-10,gtol=1e-10)
        fitted=dict(m,k_h=float(np.exp(pkfit.x[0])),scale_L_inv=float(np.exp(pkfit.x[1])))
        rows=[dict(id=m['ID'],dvid=2,time=float(row.time),dv=float(row.value)) for row in d[d.kind=='PCA'].itertuples()]
        R0=next(row['dv']/100 for row in rows if row['time']==0)
        (kout,K),pdfit=refit_war_baseline(fitted,rows,R0)
        fitted.update(kout_h=float(kout),IC50_mg_L=float(K))
        models['Warfarin'].append(WarModel(fitted,R0))
        fitting.append(dict(ID=m['ID'],drug='Warfarin',PK_success=bool(pkfit.success),PD_success=bool(pdfit['success']),PK_RMSE=float(np.sqrt(np.mean(pkfit.fun**2))),PD_RMSE=float(pdfit['rmse_pp'])))
    check('All sequential fits succeeded',all(r['PK_success'] and r['PD_success'] for r in fitting))
    save(out/'fitted_parameters.json',{drug:[model.m for model in ms] for drug,ms in models.items()})
    for drug,ms in models.items():
        for model in ms:
            print('Optimizing',drug,model.m['ID'],flush=True)
            if drug=='Remifentanil':
                base=model.objective(np.ones(model.n)); classes=[('R1',None)]
                classes += [(f'R{cap:g}',model.solve(cap=cap)) for cap in [1.25,1.5,2.]]
            else:
                base=model.objective([0],[1]);classes=[('single',None)]
                classes += [(f'R{spacing}',model.solve(spacing=spacing)) for spacing in [24,12,6]]
            previous,previous_gap=base,0.
            for name,r in classes:
                J=base if r is None else r['J']
                gap=0. if r is None else r.get('gap',r.get('sampled_relaxation_supporting_gap'))
                if gap is None: raise RuntimeError('Supporting bound failed')
                row=dict(ID=model.m['ID'],drug=drug,family=name,J=J,unit='min' if drug=='Remifentanil' else 'h',increment=J-previous,lower_increment=J-previous-previous_gap,gain_percent=100*(J/base-1),gap=gap)
                results.append(row)
                if r:
                    check(f'{model.m["ID"]} {name} solve',r['solver_success'])
                    check(f'{model.m["ID"]} {name} amount',r['amount_error']<1e-8)
                    check(f'{model.m["ID"]} {name} exposure',max(r['Cp_ratio'],r.get('Ce_ratio',1),r.get('tail_ratio',1))<=1+1e-7)
                    check(f'{model.m["ID"]} {name} independent ODE',r['independent_objective_error']<1e-6,r['independent_objective_error'])
                    check(f'{model.m["ID"]} {name} nested value',J>=previous-1e-7)
                    controls.append(dict(drug=drug,family=name,**r))
                previous,previous_gap=J,gap
            free=solve_floor(model,0,drug)
            for floor in [.05,.25]:
                continuation=solve_floor(model,floor,drug)
                delta=free['J']-continuation['J']
                bound=delta-continuation['supporting_gap']-free['independent_objective_error']-continuation['independent_objective_error']
                floor_rows.append(dict(ID=model.m['ID'],drug=drug,floor=floor,increment=delta,lower_increment=bound,positive_lower_bound=bound>1e-6))
                check(f'{model.m["ID"]} floor {floor} nested',delta>=-1e-7)
                controls.append(dict(drug=drug,ID=model.m['ID'],family='continuation',floor=floor,free=free,continuation=continuation))
    pd.DataFrame(fitting).to_csv(tables/'fit_diagnostics.csv',index=False)
    res=pd.DataFrame(results);res.to_csv(tables/'nested_classes.csv',index=False)
    pd.DataFrame(floor_rows).to_csv(tables/'continuation.csv',index=False)
    save(out/'controls.json',controls)
    c=canonical();save(out/'canonical.json',c)
    check('Canonical analytical objective',abs(c['J']-28.826261662949886)<1e-8,c['J'])
    check('Canonical amount',abs(c['B']-12.613903616328)<1e-9)
    check('Delayed injury peak',abs(c['delayed_injury_peak']-.668740304976)<1e-9)
    run_frontier(out)
    frontier=json.loads((out/'B_frontier_final_verification.json').read_text())
    check('Dose-frontier independent checks',frontier['all_passed'])
    mtx=pd.read_csv(inputs/'mtx_courses.csv');pivot=mtx.pivot(index='ID',columns='course',values='MTX_PG2_5_nmol_L')
    pairs=(pivot[4]-pivot[1]).dropna();rng=np.random.default_rng(seed+1)
    interval=np.quantile(np.median(rng.choice(pairs.to_numpy(),size=(3000,len(pairs)),replace=True),axis=1),[.025,.975])
    save(out/'mtx_summary.json',dict(synthetic=True,available_per_course=mtx.groupby('course').MTX_PG2_5_nmol_L.count().to_dict(),paired_n=len(pairs),median_change=float(pairs.median()),bootstrap_95=interval.tolist(),interpretation='Descriptive synthetic persistence; no stopping estimand'))
    for drug,ms in models.items():
        fig,ax=plt.subplots(figsize=(7,4))
        for ID,g in res[res.drug==drug].groupby('ID',sort=False):
            ax.plot(g.family,g.gain_percent,'o-',label=ID)
        ax.set(ylabel='Synthetic relative proxy gain (%)',title=f'{drug}: matched amount')
        ax.legend(fontsize=7);fig.tight_layout();fig.savefig(figures/f'{drug}_nested.png',dpi=160,metadata={"Software": None});plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4));t=np.linspace(0,4,500)
    C=np.where(t<=c['cessation'],3.,3*np.exp(-(t-c['cessation'])))
    axs[0].plot(t,C,label='Flexible');axs[0].plot(t,c['B']/4*(1-np.exp(-t)),label='Constant');axs[0].set(xlabel='Time',ylabel='C',title='Quadratic theoretical reference');axs[0].legend()
    tt=np.linspace(0,12,500);axs[1].plot(tt,(np.exp(-.2*tt)-np.exp(-tt))/.8);axs[1].axhline(.3,color='red',ls='--');axs[1].set(xlabel='Time',ylabel='Z',title='Separate unit-bolus delayed injury')
    fig.tight_layout();fig.savefig(figures/'quadratic_and_delayed_injury.png',dpi=160,metadata={"Software": None});plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4))
    for _,g in mtx.groupby('ID'):ax.plot(g.course,g.MTX_PG2_5_nmol_L,'o-',alpha=.5)
    ax.set(xlabel='Course (categorical order)',ylabel='Synthetic MTX-PG2–5 (nmol/L)',title='Synthetic biomarker persistence');ax.set_xticks([1,2,3,4]);fig.tight_layout();fig.savefig(figures/'mtx_persistence.png',dpi=160,metadata={"Software": None});plt.close(fig)
    save(out/'verification.json',dict(all_passed=all(x['passed'] for x in checks),checks=checks))
    hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).parent.glob('*.py'))}
    save(out/'run_manifest.json',dict(source_sha256=source_hashes,version=__version__,synthetic=True,seed=seed,profiles_per_drug=n,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=f"{platform.system()}-{platform.machine()}",sha256=hashes))
    if not all(x['passed'] for x in checks):raise RuntimeError('Verification failed; see verification.json')
    print('Completed:',out,flush=True)


def main():
    p=argparse.ArgumentParser(description='Reproduce the synthetic LOGOS finite-course workflow.')
    p.add_argument('--output',type=Path,required=True);p.add_argument('--profiles',type=int,default=4);p.add_argument('--seed',type=int,default=SEED)
    a=p.parse_args()
    if a.profiles<2:p.error('--profiles must be >=2')
    run(a.output,a.profiles,a.seed)

if __name__=='__main__':main()
