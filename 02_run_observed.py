#!/usr/bin/env python3
"""Audited AS reassignment Core/standard-IS BayPass 3.1 workflow.

No biological or user-approval claim is made by the numerical Omega gate.
The predesignated Core seed 41001 supplies fixed Omega and Pi-Beta means.
Paired leave-one-population covariates are separate IS tests, not a joint model.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, itertools, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('pod_runner', ROOT/'scripts/12_run_final_pod.py')
POD = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(POD)
CASES = ('full8', 'strict8', 'no_TYA', 'no_MYO')
SEEDS = {'core': (41001,41002,41003), 'gea': (51001,51002,51003)}
MCMC = dict(npilot=20,pilotlength=500,burnin=5000,nval=1000,thin=20,nthreads=1)
GRID = dict(minbeta=-0.3,maxbeta=0.3,nbetagrid=201)
sha, stamp, atomic_json, atomic_tsv = POD.sha, POD.stamp, POD.atomic_json, POD.atomic_tsv

def info(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size}

def case_paths(case):
    base=ROOT/'data/cases'/case
    files={k:base/v for k,v in dict(counts='genotypes.geno',locus_map='locus_map.tsv',
        order='population_order.tsv',environment='environment.cov',config='config.json').items()}
    config=json.loads(files['config'].read_text())
    order=POD.population_order(files['order'])
    counts=POD.numeric_matrix(files['counts']); env=POD.numeric_matrix(files['environment'])
    if config['population_order'] != order or counts.shape[1] != 2*len(order):
        raise ValueError(f'{case}: population order/count shape mismatch')
    if env.shape != (len(config['covariates']),len(order)) or np.any(np.std(env,axis=1)==0):
        raise ValueError(f'{case}: environmental shape/variation mismatch')
    if np.any(counts<0) or np.any(counts!=np.floor(counts)):
        raise ValueError(f'{case}: counts are not nonnegative integers')
    loci=POD.tsv(files['locus_map'])
    if [int(r['MRK']) for r in loci] != list(range(1,len(counts)+1)):
        raise ValueError(f'{case}: locus-map MRK differs from counts')
    return files,config,order,counts,loci

def run_job(case,phase,seed):
    files,config,order,counts,loci=case_paths(case)
    directory=ROOT/'results/02_observed'/case/phase;directory.mkdir(parents=True,exist_ok=True)
    prefix=directory/f'ASR_{case}_{phase}_s{seed}'
    binary=ROOT/'software/vendor/g_baypass'
    if sha(binary)!=POD.BINARY_SHA256: raise ValueError('BayPass binary hash mismatch')
    command=[str(binary),'-countdatafile',str(files['counts']),'-outprefix',str(prefix),'-seed',str(seed)]
    for key,value in MCMC.items(): command += ['-'+key,str(value)]
    inputs={k:info(v) for k,v in files.items()}
    inputs['binary']=info(binary)
    omega=ROOT/'results/02_observed'/case/'core'/f'ASR_{case}_core_s41001_mat_omega.out'
    beta=omega.with_name(f'ASR_{case}_core_s41001_summary_beta_params.out')
    if phase=='gea':
        gate=ROOT/'results/02_observed'/case/'omega_review/CORE_NUMERICAL_REVIEW.json'
        if json.loads(gate.read_text())['status']!='PASS': raise ValueError('Omega numerical gate not passed')
        a,b=POD.beta_params(beta)
        command += ['-efile',str(files['environment']),'-omegafile',str(omega),'-setpibetapar','-betapiprior',str(a),str(b)]
        for key,value in GRID.items():command += ['-'+key,str(value)]
        if config['nocovscaling']:command += ['-nocovscaling']
        inputs.update(omega=info(omega),beta_params=info(beta),numerical_gate=info(gate))
    record={'case':case,'phase':phase,'seed':seed,'command':command,'inputs':inputs}
    complete=prefix.with_name(prefix.name+'.completed.json')
    if complete.exists():
        done=json.loads(complete.read_text())
        if any(done.get(k)!=v for k,v in record.items()) or done['return_code']!=0:
            raise ValueError(f'Cannot resume changed command/input: {complete}')
        for value in done['outputs'].values():
            if info(value['path'])!=value:raise ValueError(f'Completed output has changed: {value["path"]}')
        return done
    old=list(directory.glob(prefix.name+'*'))
    if old:raise ValueError(f'Refusing to overwrite existing partial outputs: {old}')
    running=prefix.with_name(prefix.name+'.running.json')
    started=stamp()
    with running.open('x') as f:json.dump({**record,'started_utc':started},f,indent=2)
    log=prefix.with_name(prefix.name+'.stdout.log')
    print('START',prefix.name,flush=True)
    with log.open('x') as f:
        f.write('COMMAND_JSON '+json.dumps(command)+'\n');f.flush()
        result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,
            env={**os.environ,'OMP_NUM_THREADS':'1','LC_ALL':'C'})
        f.flush();os.fsync(f.fileno())
    if result.returncode:raise RuntimeError(f'BayPass failed {result.returncode}: {log}')
    if phase=='core':
        mat=POD.numeric_matrix(str(prefix)+'_mat_omega.out')
        if mat.shape!=(len(order),len(order)):raise ValueError('Core Omega dimension mismatch')
        POD.beta_params(Path(str(prefix)+'_summary_beta_params.out'))
    else:
        actual=POD.read_is(Path(str(prefix)+'_summary_betai_reg.out'),len(counts))
        if sorted(actual)!=list(range(1,len(config['covariates'])+1)):raise ValueError('IS covariate count mismatch')
    outputs={p.name:info(p) for p in directory.glob(prefix.name+'*') if p.is_file() and p!=running}
    done={**record,'started_utc':started,'completed_utc':stamp(),'return_code':0,
          'outputs':outputs,'runner_sha256':sha(Path(__file__))}
    atomic_json(complete,done);running.unlink()
    print('COMPLETE',prefix.name,flush=True)
    return done

def write_manifest(case,phase,done):
    files,config,order,counts,loci=case_paths(case)
    rows=[]
    for result in sorted(done,key=lambda r:r['seed']):
        command=result['command'];prefix=command[command.index('-outprefix')+1]
        row={'phase':phase,'seed':result['seed'],'status':'COMPLETED','return_code':0,
             'outprefix':prefix,'command_json':json.dumps(command),'population_order':','.join(order),
             'n_loci':len(counts),'n_populations':len(order)}
        for key,filekey in [('count','counts'),('population_order','order'),('environment','environment'),
                            ('environment_population_order','order')]:
            row[key+'_file']=str(files[filekey]);row[key+'_sha256']=sha(files[filekey])
        if phase=='gea':
            for key,inputkey in [('omega','omega'),('beta_prior_source','beta_params')]:
                row[key+'_file']=result['inputs'][inputkey]['path'];row[key+'_sha256']=result['inputs'][inputkey]['sha256']
        rows.append(row)
    atomic_tsv(ROOT/'results/02_observed'/case/phase/'command_manifest.tsv',rows)

def review_omega(case):
    files,config,order,counts,loci=case_paths(case)
    directory=ROOT/'results/02_observed'/case/'core'
    paths=[directory/f'ASR_{case}_core_s{s}_mat_omega.out' for s in SEEDS['core']]
    matrices=[np.loadtxt(p,ndmin=2) for p in paths]
    minimum_eigenvalue=min(float(np.linalg.eigvalsh(m).min()) for m in matrices)
    symmetric=all(np.allclose(m,m.T,atol=1e-8,rtol=0) for m in matrices)
    pairs=[];ix=np.triu_indices(len(order),1)
    for i,j in itertools.combinations(range(3),2):
        a,b=matrices[i],matrices[j]
        pairs.append({'seed1':SEEDS['core'][i],'seed2':SEEDS['core'][j],
            'offdiagonal_correlation':float(np.corrcoef(a[ix],b[ix])[0,1]),
            'relative_frobenius':float(np.linalg.norm(a-b)/((np.linalg.norm(a)+np.linalg.norm(b))/2))})
    corr=min(x['offdiagonal_correlation'] for x in pairs);dist=max(x['relative_frobenius'] for x in pairs)
    passed=symmetric and minimum_eigenvalue>0 and np.isfinite(corr) and corr>=.9 and dist<=.25
    result={'status':'PASS' if passed else 'FAIL','case':case,'population_order':order,
        'selected_seed_prespecified':41001,'symmetric':symmetric,'minimum_eigenvalue':minimum_eigenvalue,
        'minimum_offdiagonal_correlation':corr,'maximum_relative_frobenius':dist,'pairs':pairs,
        'thresholds':{'minimum_correlation':.9,'maximum_relative_frobenius':.25},
        'source_files':[info(p) for p in paths],
        'interpretation':'Operational numerical reproducibility gate only; not formal MCMC convergence or user approval.'}
    out=ROOT/'results/02_observed'/case/'omega_review';out.mkdir(parents=True,exist_ok=True)
    atomic_json(out/'CORE_NUMERICAL_REVIEW.json',result)
    if not passed:raise ValueError(f'{case}: predeclared Omega reproducibility gate failed: {result}')
    print('OMEGA_PASS',case,corr,dist,flush=True)

def summarize(case):
    files,config,order,counts,loci=case_paths(case)
    directory=ROOT/'results/02_observed'/case/'gea'
    runs=[POD.read_is(directory/f'ASR_{case}_gea_s{s}_summary_betai_reg.out',len(counts)) for s in SEEDS['gea']]
    for cov in sorted(runs[0]):
        rows=[]
        for i,locus in enumerate(loci):
            bfs=[x[cov][i,0] for x in runs];betas=[x[cov][i,1] for x in runs]
            rows.append({**locus,'case':case,'covariate_index':cov,
                **{f'BFdB_run{j+1}':bfs[j] for j in range(3)},
                'BFdB_median':float(np.median(bfs)),'BFdB_range':float(np.ptp(bfs)),
                'Beta_is_median':float(np.median(betas)),'Beta_sign_allele':'REF_first_count'})
        atomic_tsv(directory/f'ASR_{case}.cov{cov}.observed.tsv',rows)
    print('OBSERVED_SUMMARY',case,len(loci),len(runs[0]),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cases',nargs='+',choices=CASES,default=list(CASES))
    p.add_argument('--phase',choices=['core','gea','all'],default='all')
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--validate-only',action='store_true')
    args=p.parse_args()
    if not 1<=args.workers<=8:p.error('workers must be 1..8')
    for case in args.cases:case_paths(case)
    if args.validate_only:
        print(json.dumps({'status':'INPUT_SCHEMA_VALIDATED','cases':args.cases}));return
    phases=['core','gea'] if args.phase=='all' else [args.phase]
    for phase in phases:
        if phase=='gea':
            for case in args.cases:review_omega(case)
        done={case:[] for case in args.cases}
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futures={ex.submit(run_job,case,phase,seed):case for case in args.cases for seed in SEEDS[phase]}
            for future in as_completed(futures):done[futures[future]].append(future.result())
        for case in args.cases:
            write_manifest(case,phase,done[case])
            if phase=='core':review_omega(case)
            else:summarize(case)

if __name__=='__main__':main()
