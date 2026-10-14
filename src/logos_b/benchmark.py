"""Quadratic analytical reference and injury-constrained dose frontier."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize, linprog, minimize_scalar, brentq
from scipy.integrate import solve_ivp, quad
from .continuation import check

def canonical():
 a=4.; T=4.; c=3.; L=float(np.arccosh(4/3)); tau=T-L; B=c*(1+tau)
 g=lambda C:4*C-C*C/2
 J=quad(lambda t:g(c if t<=tau else c*np.exp(-(t-tau))),0,T,points=[tau],epsabs=1e-12)[0]
 Jref=quad(lambda t:g(B/T*(1-np.exp(-t))),0,T,epsabs=1e-12)[0]
 eps=.5
 Le=brentq(lambda le:(eps+(a-eps)/np.cosh(le))*(1+T-le)+eps*le-B,1e-8,T)
 ce=eps+(a-eps)/np.cosh(Le); te=T-Le
 Jfloor=quad(lambda t:g(ce if t<=te else eps+(ce-eps)*np.exp(-(t-te))),0,T,points=[te],epsabs=1e-12)[0]
 peak=np.log(5)/.8
 return dict(B=B,plateau=c,cessation=tau,coast=L,J=J,Jref=Jref,Jfloor=Jfloor,gain_percent=100*(J/Jref-1),redistribution_value=Jfloor-Jref,stopping_value=J-Jfloor,delayed_injury_at_0_2=(np.exp(-.04)-np.exp(-.2))/.8,delayed_injury_peak=(np.exp(-.2*peak)-np.exp(-peak))/.8,delayed_injury_peak_time=peak)

def run_frontier(OUT):
 checks=[]
 def check(name,error,tol):
  checks.append(dict(check=name,error=float(error),tolerance=tol,passed=bool(np.isfinite(error) and error<=tol)))
 def save(name,data):
  (OUT/(name+'.json')).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
 T=4.;rho=.2;n=80;edges=np.linspace(0,T,n+1);dt=T/n;cap=2.5;slack=1e-6
 q,w=np.polynomial.legendre.leggauss(12);ts=np.concatenate([(a+b)/2+(b-a)/2*q for a,b in zip(edges[:-1],edges[1:])]);weights=np.tile(w*dt/2,n)
 def cmatrix(t):
  t=np.atleast_1d(t);s=np.maximum(t[:,None]-edges[:-1],0);e=np.maximum(t[:,None]-edges[1:],0);return np.column_stack([np.exp(-t),np.exp(-e)-np.exp(-s)])
 def zatom(t):t=np.asarray(t);return (np.exp(-rho*t)-np.exp(-t))/(1-rho)
 def zstep(t):t=np.maximum(np.asarray(t),0);return ((-np.expm1(-rho*t))/rho-(-np.expm1(-t)))/(1-rho)
 def zmatrix(t):
  t=np.atleast_1d(t);return np.column_stack([zatom(t),zstep(t[:,None]-edges[:-1])-zstep(t[:,None]-edges[1:])])
 M=cmatrix(ts);amount=np.r_[1,np.full(n,dt)]
 def benefit(x):c=M@x;return float(weights@(4*c-c*c/2)),M.T@(weights*(4-c))
 def injury_peak(x):
  dz=lambda t:float(((cmatrix([t])-rho*zmatrix([t]))@x).item())
  stops=np.r_[edges,20];roots=[brentq(dz,a,b,xtol=1e-12) for a,b in zip(stops[:-1],stops[1:]) if dz(a)*dz(b)<0];pts=np.unique(np.r_[stops,roots]);values=zmatrix(pts)@x;return float(values.max()),float(pts[np.argmax(values)]),roots,dz(20)
 # Full constant family: analytic objective plus independent full-tail peak factor.
 I1=T+np.expm1(-T);I2=T+2*np.expm1(-T)-.5*np.expm1(-2*T)
 unit=np.r_[0,np.ones(n)];unitpeak,_,_,_=injury_peak(unit);safe_rate=cap/unitpeak;reference_max=4*safe_rate*I1-.5*safe_rate**2*I2
 referenceJ=lambda u:4*u*I1-.5*u*u*I2
 start=np.r_[1.,np.full(n,.5)]
 frontier=[];frontier_controls=[]
 for target in [10.,11.,12.,13.]:
  grid=np.unique(np.r_[edges,np.linspace(0,20,401)]);x=start.copy()
  for iteration in range(20):
   Z=zmatrix(grid);res=minimize(lambda x:(float(amount@x),amount),x,jac=True,bounds=[(0,20)]*(n+1),constraints=[{'type':'ineq','fun':lambda x:benefit(x)[0]-target,'jac':lambda x:benefit(x)[1]},{'type':'ineq','fun':lambda x:cap-slack-Z@x,'jac':lambda x:-Z},{'type':'ineq','fun':lambda x:12.613903616328-amount@x,'jac':lambda x:-amount}],method='SLSQP',options={'ftol':2e-12,'maxiter':2000});x=res.x
   if not res.success: raise RuntimeError(res.message)
   peak,tpeak,roots,tail_deriv=injury_peak(x)
   if peak<=cap-slack+1e-9:break
   grid=np.unique(np.r_[grid,roots])
  x*=1+1e-10  # Positive target slack, with injury slack preserved.
  peak,tpeak,roots,tail_deriv=injury_peak(x)
  J,grad=benefit(x);dose=float(amount@x)
  # Concavity gives a necessary linearized target condition; relaxing path to grid yields a lower dose bound.
  lp=linprog(amount,A_ub=np.vstack([zmatrix(grid),-grad,amount]),b_ub=np.r_[np.full(len(grid),cap),-(target-J+grad@x),12.613903616328],bounds=[(0,20)]*(n+1),method='highs');assert lp.success
  state=np.r_[x[0],0.,0.];odepeak=0.
  for j,(left,right) in enumerate(zip(np.r_[edges,20][:-1],np.r_[edges,20][1:])):
   rate=x[j+1] if j<n else 0.
   sol=solve_ivp(lambda t,y:[-y[0]+rate,y[0]-rho*y[1],4*y[0]-y[0]**2/2 if j<n else 0.],(left,right),state,method='DOP853',rtol=2e-12,atol=1e-13,dense_output=True,max_step=.1);assert sol.success;state=sol.y[:,-1]
   deriv=lambda t:sol.sol(t)[0]-rho*sol.sol(t)[1];pts=[left,right]
   if deriv(left)*deriv(right)<0:pts.append(brentq(deriv,left,right))
   odepeak=max(odepeak,float(np.max(sol.sol(pts)[1])))
  cDose=T*brentq(lambda u:referenceJ(u)-target,0,safe_rate) if target<=reference_max else None
  check('frontier declining future tail '+str(target),max(0,tail_deriv),0.)
  check('frontier available cap '+str(target),max(0,dose-12.613903616328),1e-8)
  check('frontier target '+str(target),max(0,target-J),1e-7);check('frontier injury '+str(target),max(0,peak-cap),1e-9);check('frontier ODE J '+str(target),abs(state[2]-J),1e-7);check('frontier ODE peak '+str(target),abs(odepeak-peak),1e-8);check('frontier lower dose bound '+str(target),max(0,lp.fun-dose),1e-7)
  row=dict(target=target,candidate_J=J,minimum_dose_candidate=dose,dose_lower_bound=float(lp.fun),finite_family_dose_gap=float(dose-lp.fun),peak_Z=peak,injury_margin=cap-peak,peak_time=tpeak,tail_derivative_at20=tail_deriv,independent_J_error=abs(state[2]-J),independent_peak_error=abs(odepeak-peak),constant_feasible=cDose is not None,constant_minimum_dose=cDose,dose_saving_percent=100*(1-dose/cDose) if cDose else None,solver_success=bool(res.success))
  frontier.append(row);frontier_controls.append(dict(**row,initial_atom=float(x[0]),rates=x[1:].tolist(),edges=edges.tolist()));print('B dose frontier target',target,'dose',dose,'constant',cDose,flush=True)
 save('B_dose_frontier',dict(objective='Minimum administered amount at a fixed target; distinct from value at a common available cap',benefit_horizon=4,full_injury_tail='Analytic continuous extrema through time20; thereafter zero input and negative two-exponential derivative with no future maximum',input_family='Initial atom plus 80 constant bins; rates bounded by20',injury_cap=cap,constant_safe_rate=safe_rate,constant_family_maximum_J=reference_max,targets=frontier));save('B_dose_frontier_controls',frontier_controls);pd.DataFrame(frontier).to_csv(OUT/'B_dose_frontier.csv',index=False)
 save('B_frontier_final_verification',dict(all_passed=all(c['passed'] for c in checks),checks=checks))
 assert all(c['passed'] for c in checks),[c for c in checks if not c['passed']]
 print('All B extensions complete; checks',len(checks),flush=True)
