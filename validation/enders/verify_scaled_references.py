"""Independent Decimal references for fixed Enders regression designs.

Reads saved binary64 inputs, never imports production fitting code. Normal
 equations use 80/120-digit arithmetic only as an independent reference.
Run --report to reproduce and check the saved references; --create writes a
new reference fixture only and refuses to overwrite any existing file.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
import argparse, hashlib, json
import numpy as np


def reference(x, y, precision):
    with localcontext() as ctx:
        ctx.prec=precision
        n,k=x.shape
        rows=[[D(float(v)) for v in row] for row in x]
        response=[D(float(v)) for v in y]
        if not k:
            rss=sum(v*v for v in response)
            return dict(params=[],covariance=[],rss=str(rss),tvalues=[],df_resid=n)
        gram=[[sum(row[i]*row[j] for row in rows) for j in range(k)] for i in range(k)]
        rhs=[sum(row[i]*v for row,v in zip(rows,response)) for i in range(k)]
        matrix=[gram[i]+[rhs[i]]+[D(int(i==j)) for j in range(k)] for i in range(k)]
        for i in range(k):
            pivot=max(range(i,k),key=lambda j:abs(matrix[j][i]))
            matrix[i],matrix[pivot]=matrix[pivot],matrix[i]
            divisor=matrix[i][i]; matrix[i]=[v/divisor for v in matrix[i]]
            for j in range(k):
                if j!=i:
                    factor=matrix[j][i];matrix[j]=[v-factor*a for v,a in zip(matrix[j],matrix[i])]
        beta=[row[k] for row in matrix]
        residual=[v-sum(a*b for a,b in zip(row,beta)) for row,v in zip(rows,response)]
        rss=sum(v*v for v in residual); variance=rss/D(n-k)
        covariance=[[variance*v for v in row[k+1:]] for row in matrix]
        tvalues=[beta[i]/covariance[i][i].sqrt() for i in range(k)]
        return dict(params=[float(v) for v in beta],covariance=[[float(v) for v in row] for row in covariance],rss=str(rss),tvalues=[float(v) for v in tvalues],df_resid=n-k)


def designs(y, trend, p):
    delta=np.diff(y);response=delta[p:];n=len(response)
    columns=[y[p:-1]]+[delta[p-j:len(delta)-j] for j in range(1,p+1)]+[np.ones(n)]
    names=['Level.L1']+[f'Diff.L{j}' for j in range(1,p+1)]+['const']
    if trend=='ct':columns.append(np.arange(1,n+1,dtype=float));names.append('trend')
    x=np.column_stack(columns)
    removed=('Level.L1','trend' if trend=='ct' else 'const')
    restricted=[i for i,name in enumerate(names) if name not in removed]
    return x,response,names,restricted,list(range(1,len(names)))


def generate(root):
    inputs_path=root/'validation/scaled_svd/adf_inputs.json'
    inputs=json.loads(inputs_path.read_text())['inputs'];records=[]
    for key in ('y','positive','negative'):
        y=np.array([float.fromhex(v) for v in inputs[key]])
        for trend in ('c','ct'):
            for p in (0,1,4):
                x,response,names,restricted,auxiliary=designs(y,trend,p)
                fits={}
                for label,columns in [('unrestricted',list(range(x.shape[1]))),('restricted',restricted),('auxiliary',auxiliary)]:
                    a=reference(x[:,columns],response,80);b=reference(x[:,columns],response,120)
                    for field in ('params','covariance','tvalues','df_resid'):assert a[field]==b[field]
                    assert float(a['rss'])==float(b['rss'])
                    b['names']=[names[j] for j in columns]
                    b['design_sha256']=hashlib.sha256(x[:,columns].astype('<f8').tobytes()).hexdigest();fits[label]=b
                with localcontext() as ctx:
                    ctx.prec=120
                    rss_u=D(fits['unrestricted']['rss']);rss_r=D(fits['restricted']['rss'])
                    phi=float((rss_r-rss_u)/2/(rss_u/D(fits['unrestricted']['df_resid'])))
                records.append(dict(key=key,trend=trend,lags=p,nobs=len(response),response_sha256=hashlib.sha256(response.astype('<f8').tobytes()).hexdigest(),fits=fits,phi=phi,precision_80_120_binary64_equal=True))
    return dict(provenance='Independent exact binary64-to-Decimal elimination at 80 and 120 digits; no production estimator imports. Phi uses unrounded 120-digit RSS and unrestricted df.',input_sha256=hashlib.sha256(inputs_path.read_bytes()).hexdigest(),cases=records)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--report',type=Path);p.add_argument('--create',type=Path);args=p.parse_args()
    if (args.report is None)==(args.create is None):p.error('Use exactly one of --report or --create')
    result=generate(args.root)
    if args.report:
        expected=json.loads((args.root/'validation/enders/scaled_references.json').read_text());assert result==expected
    with (args.create or args.report).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('18 independent Enders references agree at 80/120 digits')

if __name__=='__main__':main()
