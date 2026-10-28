"""Portable per-cell Monte Carlo summaries, retaining all exact-law aliases."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta,norm
from .runner import load_rows
from .engine import METHODS,ALPHAS,SPECS


def wilson(k,n):
    if not n:return None,None
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (0. if k==0 else float(max(0,mid-half))), (1. if k==n else float(min(1,mid+half)))


def exact_interval(k,n,confidence=.975):
    if not n:return None,None
    a=(1-confidence)/2
    return float(beta.ppf(a,k,n-k+1)) if k else 0.,float(beta.ppf(1-a,k+1,n-k)) if k<n else 1.


def csv_write(path,rows):
    path=Path(path)
    with gzip.open(path,'wt',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def summarize(source,destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    manifest,rows=load_rows(source)
    if not manifest['complete']:raise ValueError('Final pilot summaries require a complete run')
    groups={c['id']:[] for c in manifest['catalogue']}
    for row in rows:groups[row['cell_id']].append(row)
    rates=[];selection=[];branches=[];paths=[];pairs=[];margins=[];failures=[];warnings=[]
    primary={('enders_original','enders_conservative'),('det_aic','det_bic'),('enders_original','adf_ct'),('enders_conservative','adf_ct'),('det_aic','adf_ct'),('det_bic','adf_ct')}
    for cell in manifest['catalogue']:
        sample=groups[cell['id']]
        assert [r['replicate'] for r in sample]==list(range(manifest['reps']))
        n=len(sample)
        D=np.array([r['decisions'] for r in sample]);S=np.array([r['selected'] for r in sample])
        fields=dict(cell_id=cell['id'],T=cell['law']['T'],rho=cell['rho'],intercept=cell['intercept'],is_null=cell['is_null'],aliases=json.dumps(cell['aliases'],sort_keys=True,separators=(',',':')))
        for i,alpha in enumerate(ALPHAS):
            for j,method in enumerate(METHODS):
                valid=D[:,i,j]>=0;nv=int(valid.sum());k=int((D[:,i,j]==1).sum());nf=n-nv
                lo,hi=wilson(k,nv);rate=k/nv if nv else None
                rates.append(dict(fields,alpha=alpha,method=method,attempted=n,valid=nv,failures=nf,rejections=k,
                                  rejection_rate=rate,mcse=float(np.sqrt(rate*(1-rate)/nv)) if nv else None,
                                  ci95_low=lo,ci95_high=hi,attempted_lower=k/n,attempted_upper=(k+nf)/n))
                for code,spec in enumerate(SPECS):
                    count=int((S[:,i,j]==code).sum());lo,hi=wilson(count,nv)
                    p=count/nv if nv else None
                    selection.append(dict(fields,alpha=alpha,method=method,specification=spec,count=count,valid=nv,frequency=p,mcse=float(np.sqrt(p*(1-p)/nv)) if nv else None,ci95_low=lo,ci95_high=hi))
            for ja,jb in itertools.combinations(range(len(METHODS)),2):
                a,b=METHODS[ja],METHODS[jb]
                if (b,a) in primary:
                    ja,jb=jb,ja
                    a,b=b,a
                A=D[:,i,ja];B=D[:,i,jb];valid=(A>=0)&(B>=0);nv=int(valid.sum())
                d=A[valid]-B[valid];delta=float(d.mean()) if nv else None
                p10=int(((A==1)&(B==0)).sum());p01=int(((A==0)&(B==1)).sum())
                if (a,b)==('enders_original','enders_conservative'):
                    assert p01==0
                    lo,hi=wilson(p10,nv);interval='Wilson for nested extra-rejection event'
                elif nv:
                    l10,u10=exact_interval(p10,nv);l01,u01=exact_interval(p01,nv)
                    lo,hi=l10-u01,u10-l01;interval='Conservative paired interval from simultaneous exact discordance bounds'
                else:lo=hi=None;interval='unavailable'
                lower=np.where(A<0,0,A)-np.where(B<0,1,B)
                upper=np.where(A<0,1,A)-np.where(B<0,0,B)
                is_primary=(a,b) in primary or (b,a) in primary
                pairs.append(dict(fields,alpha=alpha,method_a=a,method_b=b,primary_5pct=is_primary and alpha==.05,
                                  common_valid=nv,discordant_a_only=p10,discordant_b_only=p01,difference=delta,
                                  mcse=float(d.std(ddof=1)/np.sqrt(nv)) if nv>1 else None,
                                  ci95_low=lo,ci95_high=hi,interval_method=interval,
                                  attempted_difference_lower=float(lower.mean()),attempted_difference_upper=float(upper.mean())))
            for variant in ['original','conservative']:
                observed=[r['enders'][str(alpha)].get(variant) for r in sample]
                observed=[r for r in observed if r is not None]
                counts=Counter((r['selected'],r['terminal'],r['reject'],r['steps'][-1]['next_action'].startswith('terminate: retain constant')) for r in observed)
                for (spec,terminal,reject,policy),count in sorted(counts.items()):
                    lo,hi=wilson(count,len(observed));p=count/len(observed)
                    branches.append(dict(fields,alpha=alpha,variant=variant,specification=spec,terminal=terminal,reject=reject,policy_stop=policy,count=count,valid=len(observed),frequency=p,mcse=float(np.sqrt(p*(1-p)/len(observed))),ci95_low=lo,ci95_high=hi))
                counts=Counter('|'.join(f"{s['specification']}:{s['test']}:{'R' if s['reject'] else 'NR'}" for s in r['steps']) for r in observed)
                for key,count in sorted(counts.items()):paths.append(dict(fields,alpha=alpha,variant=variant,path=key,count=count,valid=len(observed),frequency=count/len(observed)))
        for criterion in ['aic','bic']:
            gaps=np.array([r['ic'][criterion]['gap'] for r in sample if criterion in r['ic']])
            values=dict(fields,criterion=criterion,valid=len(gaps),exact_ties=int((gaps==0).sum()),gap_below_0_1=int((gaps<.1).sum()))
            for name,q in [('min',0),('p01',.01),('p10',.1),('median',.5),('p90',.9),('p99',.99),('max',1)]:values[name]=float(np.quantile(gaps,q)) if len(gaps) else None
            values['mean']=float(gaps.mean()) if len(gaps) else None
            margins.append(values)
        for row in sample:
            failures.extend(dict(fields,replicate=row['replicate'],**f) for f in row['failures'])
            warnings.extend(dict(fields,replicate=row['replicate'],**v) for v in row['warnings'])
    tables=dict(rejection_rates=rates,selection_frequencies=selection,enders_terminals=branches,enders_paths=paths,paired_differences=pairs,ic_margins=margins)
    for name,table in tables.items():csv_write(destination/(name+'.csv.gz'),table)
    (destination/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
    (destination/'warnings.json').write_text(json.dumps(warnings,indent=2)+'\n')
    done=[json.loads(p.read_text()) for p in (Path(source)/'shards').glob('*.done.json')]
    timing=[]
    for c in manifest['catalogue']:
        chunks=[d for d in done if d['cell_id']==c['id']]
        timing.append(dict(cell_id=c['id'],T=c['law']['T'],compute_seconds=sum(d['compute_seconds'] for d in chunks),elapsed_task_seconds=sum(d['elapsed_seconds'] for d in chunks)))
    csv_write(destination/'cell_timings.csv.gz',timing)
    report=dict(complete=True,identity_hash=manifest['identity_hash'],cells=len(groups),null_cells=sum(c['is_null'] for c in manifest['catalogue']),replications_per_cell=manifest['reps'],cell_replications=len(rows),method_alpha_decisions=int(np.prod(D.shape[1:])*len(rows)),
                failure_records=len(failures),warning_records=len(warnings),equivalence_checks=sum(r['equivalence_checks'] for r in rows),
                nesting_violations=0,source=str(Path(source).resolve()),table_rows={name:len(table) for name,table in tables.items()},
                caution='Exploratory pilot; per-cell Monte Carlo intervals are pointwise. No overall nominal-size claim, no main experiment, no size-adjusted calibration.')
    (destination/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();print(json.dumps(summarize(args.source,args.output),indent=2))


if __name__=='__main__':main()
