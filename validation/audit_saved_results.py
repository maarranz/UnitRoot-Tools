#!/usr/bin/env python3
"""Read-only Stage 2C research audit; no production UnitRoot-Tools code.

Requires NumPy, SciPy and released ARCH 8.0.0. Reads all saved Stage 2B
results; writes only an explicitly requested --report path. Does not import
cross_validate.py, require R/development ARCH, or change upstream packages.
"""
from pathlib import Path
import json,csv,hashlib,collections,math,argparse
import numpy as np
from scipy.linalg import lstsq, solve_triangular, qr
parser=argparse.ArgumentParser(description='Independently audit saved Stage 2 evidence; no upstream packages are modified.')
parser.add_argument('--report', type=Path)
args=parser.parse_args()
if not __debug__:
 raise RuntimeError('Run without -O: this research audit uses assertions')
ROOT=Path(__file__).resolve().parents[1]; P=ROOT/'validation/results'
def read(n):
 with (P/n).open(newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
all_hashes={p.name:sha(p) for p in sorted(P.iterdir()) if p.is_file()}
js={p.name:json.loads(p.read_text()) for p in P.glob('*.json')}
manifest=js['manifest.json']
assert manifest['tolerance']=={'absolute':1e-8,'selected_lags':'exact'}
assert js['python_release.csv.metadata.json']['arch_runtime_version']=='8.0.0'
assert js['python_development.csv.metadata.json']['arch_unitroot_sha256']==manifest['arch_source_sha256']
assert {x['name']:x['version'] for x in js['r_metadata.json']['packages']}['bootUR']=='1.0.5'
for n,h in manifest['file_sha256'].items():assert sha(P/n)==h,n
for n,h in manifest['script_sha256'].items():assert sha(ROOT/'validation'/n)==h,n
C=read('comparisons.csv'); assert len(C)==1188
for r in C:
 a,b=float(r['a']),float(r['b']);assert np.isfinite([a,b]).all()
 delta=abs(a-b); assert delta==float(r['abs_difference'])
 assert (delta<=1e-8)==(r['within_1e_8']=='True')
 assert r['expected_equal'] in ['True','False']
expected=[r for r in C if r['expected_equal']=='True'];assert len(expected)==840
assert all(r['within_1e_8']=='True' for r in expected)
for g,s in js['summary.json'].items():
 z=[r for r in C if r['group']==g]
 actual=dict(comparisons=len(z),max_abs_difference=max(float(r['abs_difference']) for r in z),differences_over_1e_8=sum(r['within_1e_8']=='False' for r in z),expected_failures=sum(r['expected_equal']=='True' and r['within_1e_8']=='False' for r in z))
 assert s==actual,g
S=read('series.csv');Y={k:np.array([float(r[k]) for r in S if r[k]]) for k in S[0]}
for k,v in Y.items():assert np.isfinite(v).all() and len(v)==manifest['lengths'][k]
# Independently reconstruct recurrences without importing the validation harness.
rng=np.random.Generator(np.random.PCG64(20261027));E=rng.normal(size=(5,501))
walk=np.add.accumulate(E[0,1:]); neg=np.add.accumulate(E[1,1:]-.8*E[1,:-1])
A=np.zeros(501);T=np.zeros(501)
for i in range(1,501): A[i]=.7*A[i-1]+E[2,i];T[i]=.8*T[i-1]+E[3,i]
G=dict(rw=walk[300:],rw_negma=neg[300:],ar07=A[301:],trend_ar08=T[301:]+2+.03*np.arange(1,201),hetero_rw=np.add.accumulate(E[4,1:201]*np.r_[np.ones(100),np.full(100,3)]))
e=np.random.Generator(np.random.PCG64(457)).normal(size=65);G['mic_boundary']=np.add.accumulate(e[1:]-.8*e[:-1])
for k in G:assert np.array_equal(G[k],Y[k]),k
# All raw rows, dimensions, uniqueness, and effective samples.
raw={n:read(n) for n in ['python_release.csv','python_development.csv','r_results.csv']}
for n,nn in [('python_release.csv',216),('python_development.csv',120),('r_results.csv',252)]:
 assert len(raw[n])==nn
 keys=[]
 for r in raw[n]:
  k=(r['series'],r['backend'],r['test'],r['trend'],r['lags'],r['variant']);keys.append(k)
  p=float(r['lags']);assert p==int(p) and p>=0
  N=len(Y[r['series']]);test=r['test']
  obs=N if test=='KPSS' else N-1 if test.startswith('PP_') else N-1-int(p)
  assert int(r['nobs'])==obs,(n,k)
  for fld in ['statistic','pvalue','cv1','cv5','cv10','estimate','criterion_value','reported_lags']:
   if r.get(fld):assert math.isfinite(float(r[fld]))
  if r['pvalue']:assert 0<=float(r['pvalue'])<=1
 assert len(keys)==len(set(keys)),n
# New QR-based calculation, explicit row indices, rather than harness mic_curve.
def D(N,tr):
 count={'n':0,'c':1,'ct':2,'ctt':3}[tr];return np.column_stack([np.linspace(1/N,1,N)**j for j in range(count)]) if count else np.empty((N,0))
def fit(X,y):
 q,r=qr(X,mode='economic');b=solve_triangular(r,q.T@y);e=y-X@b
 return b,e

def residual(y,tr,gls=False):
 z=D(len(y),tr)
 if not z.shape[1]:return y.copy()
 yz=y.copy();zz=z.copy()
 if gls:
  a=1-({'c':7,'ct':13.5}[tr])/len(y)
  yz[1:]=y[1:]-a*y[:-1];zz[1:]=z[1:]-a*z[:-1]
 b,_=fit(zz,yz);return y-z@b

candidates=read('mic_candidates.csv');assert len(candidates)==1620
errors=collections.defaultdict(float);curves={}
for r in candidates:
 id,tr,K,p,v=r['series'],r['trend'],int(r['K']),int(r['p']),r['variant']
 y=Y[id];N=len(y);one=v.startswith('one_step')
 u=y if one else residual(y,tr,gls=(v=='NP_GLS'))
 # 1-based response rows t=K+2,...,N; convert each index explicitly.
 ix=np.arange(K+1,N);resp=u[ix]-u[ix-1];level=u[ix-1]
 xcols=[level]+[u[ix-j]-u[ix-j-1] for j in range(1,p+1)]
 z=D(len(ix),tr) if one else np.empty((len(ix),0))
 X=np.column_stack(xcols+[z]);beta,e=fit(X,resp);n=N-1-K
 if v.endswith('bootUR'):ss=level[1:]@level[1:]
 elif one:ss=(level-z@fit(z,level)[0])@(level-z@fit(z,level)[0]) if z.shape[1] else level@level
 else:ss=level@level
 s2=e@e/n;tau=beta[0]**2*ss/s2
 val=dict(n=n,beta=beta[0],sigma2=s2,level_ss=ss,tau=tau,MAIC=np.log(s2)+2*(p+tau)/n,MBIC=np.log(s2)+np.log(n)*(p+tau)/n)
 for fld,num in val.items():
  diff=abs(float(r[fld])-num);errors[fld]=max(errors[fld],float(diff));assert diff<=1e-8,(id,tr,K,p,v,fld,diff)
 curves.setdefault((id,tr,K,v),[]).append((p,val))
# Resolve all saved MIC selections against independent QR calculation.
selection_checks=0
for n in raw:
 for r in raw[n]:
  test=r['test']
  if not(test.startswith('MAIC') or test=='MBIC') or r['variant'].endswith('scaled') and not r['variant'].endswith('unscaled'):continue
  id,tr=r['series'],r['trend'];K=int(r['variant'].split('_')[0][1:]) if r['backend']=='ARCH_dev' else int(r['variant'].split('_')[2][1:])
  variant='one_step_ARCH' if test=='MAIC_ADF' else 'PQ_full' if test=='MAIC_DFGLS' else 'one_step_bootUR' if r['variant'].startswith('one_step') else 'PQ_bootUR'
  crit='MBIC' if test=='MBIC' else 'MAIC';p,z=min(curves[(id,tr,K,variant)],key=lambda pz:(pz[1][crit],pz[0]))
  assert int(float(r['lags']))==p,(n,id,tr,K,test)
  if r.get('criterion_value'):assert abs(float(r['criterion_value'])-z[crit])<1e-8
  if r.get('reported_lags'):
   extra={'n':0,'c':1,'ct':2}[tr] if r['variant'].startswith('one_step') else 0
   assert int(float(r['reported_lags']))==p+extra
  selection_checks+=1
# Independent OLS t statistics and KPSS from saved input, plus PP reconstructions.
stat_err=collections.defaultdict(float)
for r in raw['python_release.csv']+raw['r_results.csv']:
 if r['variant']!='fixed':continue
 id,tr=r['series'],r['trend'];y=Y[id];N=len(y);p=int(r['lags']);test=r['test']
 if test in ['ADF','DFGLS']:
  v=residual(y,tr,True) if test=='DFGLS' else y
  ix=np.arange(p+1,N);resp=v[ix]-v[ix-1]
  X=np.column_stack([v[ix-1]]+[v[ix-j]-v[ix-j-1] for j in range(1,p+1)]+[D(len(ix),tr) if test=='ADF' else np.empty((len(ix),0))])
  b,e=fit(X,resp);q,R=qr(X,mode='economic');invR=solve_triangular(R,np.eye(R.shape[0]));cov=invR@invR.T;stat=b[0]/np.sqrt(e@e/(len(ix)-X.shape[1])*cov[0,0])
 elif test=='KPSS':
  u=residual(y,tr);L=u@u/N
  for h in range(1,p+1):L+=2*(1-h/(p+1))*(u[h:]@u[:-h])/N
  stat=np.sum(np.cumsum(u)**2)/N**2/L
 elif test.startswith('PP_'):
  n=N-1;lhs=y[1:];level=y[:-1];z=D(n,tr);X=np.column_stack([level,z]);b,u=fit(X,lhs)
  q,R=qr(X,mode='economic');invR=solve_triangular(R,np.eye(R.shape[0]));c=invR@invR.T;g=u@u/n;L=g
  for h in range(1,p+1):L+=2*(1-h/(p+1))*(u[h:]@u[:-h])/n
  se=np.sqrt(u@u/(n-X.shape[1])*c[0,0])
  lev=lhs if r['backend']=='R' else level;levres=lev-z@fit(z,lev)[0];M=levres@levres/n**2
  if r['backend']=='R' and tr=='ct':M*=1-n**-2
  stat=(np.sqrt(g/L)*(b[0]-1)/se-(L-g)/(2*np.sqrt(L*M))) if test=='PP_tau' else n*(b[0]-1)-(L-g)/(2*M)
 else:raise ValueError(test)
 err=abs(stat-float(r['statistic']));stat_err[test]=max(stat_err[test],float(err));assert err<1e-8,(test,id,tr,p,err)
# Ensure direct package comparison rows have correct raw provenance.
flat=sum(raw.values(),[])
def lookup(backend,r,test=None):
 x=[a for a in flat if a['series']==r['series'] and a['backend']==backend and a['test']==(test or r['group']) and a['trend']==r['trend'] and a['lags']==r['lags'] and a['variant']=='fixed'];assert len(x)==1;return x[0]
for r in C:
 if r['group'] in ['ADF','KPSS','DFGLS','PP_tau','PP_rho']:
  left=lookup('ARCH',r);right=lookup('Statsmodels' if r['group'] in ['ADF','KPSS'] else 'R',r)
  assert float(left[r['field']])==float(r['a']) and float(right[r['field']])==float(r['b'])
# boundary deltas exactly match omitted first squared level contribution.
boundary=[]
for p,a in curves[('mic_boundary','c',4,'PQ_full')]:
 b=dict(curves[('mic_boundary','c',4,'PQ_bootUR')])[p];u=residual(Y['mic_boundary'],'c')
 for criterion,C0 in [('MAIC',2),('MBIC',np.log(59))]:assert abs((a[criterion]-b[criterion])-C0/59*a['beta']**2*u[4]**2/a['sigma2'])<1e-12
 boundary.append(dict(p=p,MAIC_full=a['MAIC'],MAIC_bootUR=b['MAIC'],MBIC_full=a['MBIC'],MBIC_bootUR=b['MBIC']))
# Check the intended production composition using only released public APIs.
import arch
from arch.unitroot import ADF, DFGLS
assert arch.__version__=='8.0.0'
integration_errors=[]
for (id,tr,K,v),curve in curves.items():
 if v!='PQ_full':continue
 for method in ['MAIC','MBIC']:
  p,detail=min(curve,key=lambda pz:(pz[1][method],pz[0]))
  for cls in ([ADF,DFGLS] if tr!='n' else [ADF]):
   y=Y[id];f=cls(y,trend=tr,lags=p)
   transformed=residual(y,tr,True) if cls is DFGLS else y
   ix=np.arange(p+1,len(y));response=transformed[ix]-transformed[ix-1]
   X=np.column_stack([transformed[ix-1]]+[transformed[ix-j]-transformed[ix-j-1] for j in range(1,p+1)]+[D(len(ix),tr) if cls is ADF else np.empty((len(ix),0))])
   b,e=fit(X,response);q,R=qr(X,mode='economic');invR=solve_triangular(R,np.eye(R.shape[0]))
   variance=e@e/(len(ix)-X.shape[1]);stat=b[0]/np.sqrt(variance*(invR@invR.T)[0,0])
   assert f.lags==p and f.nobs==len(y)-1-p
   error=abs(f.stat-stat);assert error<1e-8
   integration_errors.append(float(error))
assert len(integration_errors)==120
report=dict(audit_script_sha256=sha(Path(__file__)),released_arch_public_composition_checks=len(integration_errors),released_arch_composition_max_error=max(integration_errors),files_read_sha256=all_hashes,total_comparisons=len(C),expected_equalities=len(expected),expected_failures=0,raw_rows={k:len(v) for k,v in raw.items()},candidate_rows=len(candidates),qr_candidate_max_errors=dict(errors),raw_selection_checks=selection_checks,qr_statistic_max_errors=dict(stat_err),boundary=boundary,metadata_files=list(js))
if args.report:
 args.report.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
