"""Sequential fits to generated observations; fixed oral absorption and lag."""
import numpy as np
from scipy import optimize, special
from .models import WarModel, quad_segments
from scipy.optimize import least_squares

def step_response(t,w,k,ke=None):
 t=np.maximum(np.asarray(t,float),0);out=np.zeros_like(t)
 for wi,ki in zip(w,k):
  if ke is None:out+=wi*(-np.expm1(-ki*t))/ki
  else:
   if abs(ke-ki)<1e-6*max(ke,ki): out+=wi*((1-np.exp(-ki*t))/ki-t*np.exp(-ki*t))
   else:out+=wi*ke/(ke-ki)*((-np.expm1(-ki*t))/ki-(-np.expm1(-ke*t))/ke)
 return out

def infusion(t,w,k,rate,T,ke=None):
 t=np.asarray(t);return rate*(step_response(t,w,k,ke)-step_response(t-T,w,k,ke))

def decode_pk(p):
 # ordered positive rates avoid component permutation
 return np.exp(p[:2]),np.array([np.exp(p[2]),np.exp(p[2])+np.exp(p[3])])

def fit_rem_pk(t,y,rate,T):
 floor=.1;lo=np.log([1e-7,1e-7,.001,.005]);hi=np.log([1,1,.5,5])
 def res(p):
  w,k=decode_pk(p);return np.log(np.maximum(infusion(t,w,k,rate,T),1e-9)+floor)-np.log(y+floor)
 fits=[]
 for slow,fast in [(.02,.6),(.05,1.5),(.008,.2)]:
  start=np.log([.001,.008,slow,fast]);fits.append(optimize.least_squares(res,start,bounds=(lo,hi),max_nfev=400,ftol=1e-10,xtol=1e-10,gtol=1e-10))
 best=min(fits,key=lambda a:np.dot(a.fun,a.fun));w,k=decode_pk(best.x)
 return w,k,best

def sigmoid(C,p):
 ke,K,h,E0,frac=np.exp(p[0]),np.exp(p[1]),np.exp(p[2]),p[3],p[4]
 eff=special.expit(h*(np.log(np.maximum(C,1e-15))-np.log(K)))
 return E0*(1-frac*eff)

def fit_rem_pd(t,y,w,k,rate,T,direct=False):
 def res(p):
  C=infusion(t,w,k,rate,T,None if direct else np.exp(p[0]));return sigmoid(C,p)-y
 lo=np.array([np.log(.01),np.log(.1),np.log(.3),8,.1]);hi=np.array([np.log(5),np.log(200),np.log(8),35,1])
 E0=np.clip(np.median(y[t<=.5]) if np.any(t<=.5) else np.max(y),9,34)
 fits=[]
 for ke in [.15,.6,1.5]:
  p=np.array([np.log(ke),np.log(15),np.log(2),E0,.8]);fits.append(optimize.least_squares(res,p,bounds=(lo,hi),max_nfev=300,ftol=1e-9,xtol=1e-9,gtol=1e-9))
 best=min(fits,key=lambda a:np.dot(a.fun,a.fun));return best

def refit_war_baseline(m,rows,Rinit):
 obs=[r for r in rows if r['id']==m['ID'] and r['dvid']==2 and r['dv'] is not None and r['time']>0];tt=np.array([r['time'] for r in obs]);yy=np.array([r['dv']/100 for r in obs]);idx=np.argsort(tt);tt=tt[idx];yy=yy[idx]
 def resid(z):return WarModel(m,Rinit,np.exp(z[0]),np.exp(z[1])).response(tt)-yy
 starts=[np.log([m['kout_h'],m['IC50_mg_L']]),np.log([.02,2]),np.log([.05,5])]
 results=[least_squares(resid,z,bounds=(np.log([1e-4,1e-3]),np.log([1,1000])),ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=300) for z in starts];best=min(results,key=lambda r:sum(r.fun**2))
 return np.exp(best.x),{'n_fit':len(tt),'rmse_pp':float(100*np.sqrt(np.mean(best.fun**2))),'success':bool(best.success),'bound_distance':float(min(best.x-np.log([1e-4,1e-3])).item())}

