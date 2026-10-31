"""Independently reproduce the four saved ADF references; writes only a new report.

Uses exact binary64 values converted to Decimal and 80/120-digit elimination.
Normal equations here are high-precision reference arithmetic, not a production
binary64 fitting method. No production helper or ARCH estimator is imported.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext
import argparse, csv, hashlib, json
import numpy as np

def decimal_reference(X,Y,precision):
 with localcontext() as ctx:
  ctx.prec=precision;n,k=X.shape;a=[[D(float(x)) for x in row] for row in X];y=[D(float(v)) for v in Y]
  G=[[sum(row[i]*row[j] for row in a) for j in range(k)] for i in range(k)];rhs=[sum(row[i]*v for row,v in zip(a,y)) for i in range(k)]
  A=[G[i]+[rhs[i]]+[D(int(i==j)) for j in range(k)] for i in range(k)]
  for i in range(k):
   pivot=max(range(i,k),key=lambda j:abs(A[j][i]));A[i],A[pivot]=A[pivot],A[i];d=A[i][i];A[i]=[v/d for v in A[i]]
   for j in range(k):
    if i==j:continue
    f=A[j][i];A[j]=[v-f*c for v,c in zip(A[j],A[i])]
  beta=[row[k] for row in A];inverse=[row[k+1:] for row in A];res=[v-sum(x*b for x,b in zip(row,beta)) for row,v in zip(a,y)];rss=sum(v*v for v in res);variance=rss/D(n-k);cov=[[variance*v for v in row] for row in inverse]
  return dict(beta=np.array([float(v) for v in beta]),rss=float(rss),cov=np.array([[float(v) for v in row] for row in cov]),rss_decimal=str(rss))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args(); root=args.root
    manifest=json.loads((root/'validation/results/manifest.json').read_text())
    path=root/'validation/results/series.csv'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==manifest['file_sha256']['series.csv']
    with path.open() as f: rows=list(csv.DictReader(f))
    references=json.loads((root/'validation/scaled_svd/adf_ctt_references.json').read_text())['cases']
    output=[]
    for case in references:
        name,p=case['series'],case['lags']
        y=np.array([float(row[name]) for row in rows[:manifest['lengths'][name]]])
        delta=np.diff(y); n=len(delta)-p; t=np.arange(1,n+1,dtype=float)
        X=np.column_stack([y[p:-1]]+[delta[p-j:len(delta)-j] for j in range(1,p+1)]+[np.ones(n),t,t*t]); Y=delta[p:]
        refs=[decimal_reference(X,Y,precision) for precision in (80,120)]
        statistics=[float(a['beta'][0]/np.sqrt(a['cov'][0,0])) for a in refs]
        assert statistics[0]==statistics[1]==case['reference_stat']
        output.append(dict(series=name,lags=p,reference_stat=statistics[1],precision_80_120_equal=True,
            design_sha256=hashlib.sha256(X.astype('<f8').tobytes()).hexdigest(),
            response_sha256=hashlib.sha256(Y.astype('<f8').tobytes()).hexdigest(),
            beta=refs[1]['beta'].tolist(),covariance=refs[1]['cov'].tolist(),rss=refs[1]['rss_decimal']))
    with args.report.open('x') as f: json.dump(output,f,indent=2); f.write('\n')
    print('Four independent 80/120-digit ADF references reproduced exactly')

if __name__=='__main__': main()
