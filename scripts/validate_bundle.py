#!/usr/bin/env python3
"""Check the compact paper bundle against its definitions and independent records."""
import argparse,hashlib,json,math,platform,re,sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
CANONICAL='bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30'

def require(condition,message):
    if not condition:raise AssertionError(message)

def crossings(t,y):
    idx=np.where(np.diff(np.sign(y))!=0)[0]
    return np.array([t[i]-y[i]*(t[i+1]-t[i])/(y[i+1]-y[i]) for i in idx if y[i+1]!=y[i]])

def validate():
    summary={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'cases':[]}
    checksums=json.loads((ROOT/'provenance/checksums.json').read_text())
    for name,digest in checksums.items():
        p=ROOT/name;require(p.is_file(),f'Missing bundle file: {name}')
        require(hashlib.sha256(p.read_bytes()).hexdigest()==digest,f'Checksum mismatch: {name}')
    summary['checksummed_files']=len(checksums)
    inputs=list((ROOT/'cases').iterdir());require(len(inputs)==19,'Expected 19 paper cases')
    for case in inputs:
        for field in ('U','Us','rho','p','Ksl','voidfraction'):
            require((case/'CFD/0'/field).is_file(),f'Missing initial field: {case.name}/{field}')
    sweep=sorted((ROOT/'data/cases').iterdir());require(len(sweep)==15,'Expected 15 sweep cases')
    for case in sweep:
        df=pd.read_csv(case/'time_series.tsv',sep='\t');t=df.time.to_numpy();gap=df.z_gap_large_minus_small.to_numpy()
        require(np.isfinite(t).all() and np.isfinite(gap).all(),f'Nonfinite height data: {case.name}')
        require(np.all(np.diff(t)>=0),f'Unsorted time: {case.name}')
        require(abs(t[-1]-40)<1e-9,f'Incomplete curve: {case.name} ends at {t[-1]} s')
        require(np.max(abs(df.z_mean_large-df.z_mean_small-gap))<1e-12,f'Height identity: {case.name}')
        require((df.count_small>0).all() and (df.count_large>0).all(),f'Missing species: {case.name}')
        params=json.loads((ROOT/'cases'/case.name/'case_params.json').read_text())
        shots=list((case/'snapshots').glob('*.npz'))
        for shot in shots:
            with np.load(shot,allow_pickle=False) as f:
                xyz=f['xyz'];typ=f['types'];rad=f['radii'];time=float(f['time']);box=f['box']
            require(xyz.shape==(len(typ),3) and rad.shape==typ.shape,f'Snapshot shapes: {shot.name}')
            require(set(np.unique(typ))=={1,2} and np.isfinite(xyz).all(),f'Snapshot species: {case.name}')
            row=df.iloc[int(np.argmin(abs(t-time)))];require(abs(float(row.time)-time)<1e-8,f'Snapshot timestamp: {case.name}')
            for k,species,key in ((1,'small','small_radius'),(2,'large','large_radius')):
                mask=typ==k
                require(mask.sum()==row['count_'+species],f'Snapshot count: {case.name}/{species}')
                require(abs(xyz[mask,2].mean()-row['z_mean_'+species])<1e-11,f'Snapshot height: {case.name}/{species}')
                require(np.max(abs(rad[mask]-params[key]))<1e-12,f'Snapshot radius: {case.name}/{species}')
            require(box.shape==(3,2),f'Snapshot box: {case.name}')
        require(bool(shots)==('ang30' in case.name),f'Snapshot coverage: {case.name}')
        summary['cases'].append({'case':case.name,'frames':len(df),'time_end_s':float(t[-1]),
                                 'dz_min_mm':float(gap.min()*1e3),'t_min_s':float(t[gap.argmin()]),'snapshots':len(shots)})
    df=pd.read_csv(ROOT/'data/cases'/CANONICAL/'time_series.tsv',sep='\t')
    independent=pd.read_csv(ROOT/'data/canonical_mean_heights.csv')
    require(len(df)==len(independent),'Canonical independent record length')
    require(np.max(abs(df.time-independent.time))<1e-9,'Canonical independent times')
    require(np.max(abs(df.z_gap_large_minus_small-independent.dz))<1e-12,'Canonical independent heights')
    zc=crossings(df.time.to_numpy(),df.z_gap_large_minus_small.to_numpy());tm=float(df.time.iloc[df.z_gap_large_minus_small.argmin()])
    t1=float(zc[zc<tm][-1]);t2=float(zc[zc>tm][0])
    metadata=json.loads((ROOT/'provenance/original_figure06_metadata.json').read_text())
    require(abs(t1-metadata['t1_s'])<1e-9 and abs(t2-metadata['t2_s'])<1e-9,'Canonical phase boundaries')
    summary['canonical_crossings_s']=[t1,t2]
    for time,key in ((1.5,'F'),(28.5,'G')):
        profile=pd.read_csv(ROOT/'data/cases'/CANONICAL/f'fluid_profile_{time:g}.csv')
        require(len(profile)==16 and np.isfinite(profile.to_numpy()).all(),'Fluid profile shape')
        require(abs(profile.ux_m_s.max()*1e3-metadata['max_ux_mm_s'][key])<1e-8,'Fluid profile maximum')
    controls=[]
    for kind in ('baseline','settling_only','shear_only','equal_settling_stokes'):
        name=f'mech_{kind}_dS200_dL600_phi0p1_ang30_seed1';d=pd.read_csv(ROOT/'data/mechanisms'/f'{name}_meanz.csv')
        require(abs(d.time.iloc[-1]-40)<1e-9,f'Incomplete mechanism: {kind}')
        require(np.max(abs(d.zbar_large-d.zbar_small-d.dz))<1e-12,f'Mechanism height identity: {kind}')
        controls.append({'control':kind,'frames':len(d),'dz_min_mm':float(d.dz.min()*1e3),'dz_final_mm':float(d.dz.iloc[-1]*1e3)})
    name='mech_equal_settling_stokes_dS200_dL600_phi0p1_ang30_seed1';case=ROOT/'cases'/name
    p=json.loads((case/'case_params.json').read_text());rho=re.search(r'internalField\s+uniform\s+([^;]+);',(case/'CFD/0/rho').read_text())
    expected=(p['rhop']-float(rho[1]))*(4/3)*math.pi*p['large_radius']**3*abs(p['gz'])*(1-(p['small_radius']/p['large_radius'])**2)
    actual=float(re.search(r'^fix\s+addlarge\s+glarge\s+addforce\s+0\.0\s+0\.0\s+(\S+)',(case/'DEM/in.liggghts_setting').read_text(),re.M)[1])
    require(abs(actual/expected-1)<1e-12,'Stokes-control body force')
    summary['mechanism_controls']=controls;summary['stokes_force_N']=actual
    summary['status']='passed';return summary

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--report',type=Path);args=parser.parse_args()
    result=validate()
    if args.report:args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(f"PASS: {result['checksummed_files']} files, 19 case inputs, 15 complete curves, snapshots, fluid profiles, and 4 mechanism controls")
if __name__=='__main__':main()
