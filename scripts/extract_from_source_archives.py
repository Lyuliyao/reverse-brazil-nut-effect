"""Extract paper inputs and minimal data from uncompressed case archives."""
import argparse, csv, io, json, re, sys, tarfile, time, shutil, hashlib
from pathlib import Path, PurePosixPath
import numpy as np
import h5py

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root', required=True, type=Path, help='CFDEM_NRE root containing CFDEM_PROJECT')
parser.add_argument('--out', required=True, type=Path, help='Destination for the compact extraction')
args=parser.parse_args()
ROOT=args.source_root.resolve()
PROJECT=ROOT/'CFDEM_PROJECT'
OUT=args.out.resolve()
OUT.mkdir(parents=True,exist_ok=True)
names=[]
for phi in ('0p05','0p1','0p2','0p3','0p4'):
    for angle,large in ((20,300),(30,300),(30,600)):
        names.append(f'bne_vf{phi}_rs100um_rl{large}um_rho2470_cor0p1_mu0p5_ang{angle}')
for kind in ('baseline','settling_only','shear_only','equal_settling_stokes'):
    names.append(f'mech_{kind}_dS200_dL600_phi0p1_ang30_seed1')
provenance=[]

def relative(member,name):
    parts=PurePosixPath(member.name).parts
    if name not in parts: return None
    return PurePosixPath(*parts[parts.index(name)+1:])

def save_frame(tf, members, name, target, label=None):
    best=None
    for m in members:
        if not m.name.endswith('/forcebalance.h5'): continue
        seg=PurePosixPath(m.name).parent.name
        if label and seg!=label: continue
        # All frame metadata is small. Only the selected particle frame is read.
        with h5py.File(tf.extractfile(m),'r') as f:
            ts=f['frames']['timestep']*1e-7
            k=int(np.argmin(abs(ts-target)));err=abs(float(ts[k])-target)
            if best is None or err<best[0]: best=(err,m,k,float(ts[k]))
    if best is None: raise ValueError('No packed frame for '+name)
    err,m,k,actual=best
    with h5py.File(tf.extractfile(m),'r') as f:
        fr=f['frames'][k];start=int(fr['row_start']);count=int(fr['row_count'])
        a=f['particles'][start:start+count]
        cols=[x.decode() if isinstance(x,bytes) else str(x) for x in f['particles'].attrs['columns']]
        ix={x:i for i,x in enumerate(cols)}
        box=f['box_bounds'][0]
    dest=OUT/'data/cases'/name/'snapshots'/f't{actual:.8g}.npz'
    dest.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest,xyz=a[:,[ix['x'],ix['y'],ix['z']]],types=a[:,ix['type']].astype(int),radii=a[:,ix['radius']],box=box,time=actual)
    return {'path':str(dest.relative_to(OUT)),'archive_member':m.name,'frame_index':k,'requested_time_s':target,'actual_time_s':actual,'row_count':count}

def complete_heights(tf, members):
    rows=[];sources=[]
    for m in sorted((m for m in members if m.name.endswith('/forcebalance.h5')), key=lambda m:m.name):
        with h5py.File(tf.extractfile(m),'r',rdcc_nbytes=16*1024*1024) as f:
            columns=[v.decode() if isinstance(v,bytes) else str(v) for v in f['particles'].attrs['columns']]
            it=columns.index('type');iz=columns.index('z')
            for frame in f['frames'][:]:
                start=int(frame['row_start']);n=int(frame['row_count']);a=f['particles'][start:start+n]
                types=a[:,it].astype(int);z=a[:,iz];zs=float(z[types==1].mean());zl=float(z[types==2].mean());step=int(frame['timestep'])
                rows.append(dict(time=step*1e-7,dem_timestep=step,count_small=int((types==1).sum()),count_large=int((types==2).sum()),z_mean_small=zs,z_mean_large=zl,z_gap_large_minus_small=zl-zs))
        sources.append(m.name)
    rows.sort(key=lambda r:r['time'])
    return rows,sources

for name in names:
    print('START '+name,flush=True)
    candidates=[PROJECT/(name+'.tar'),PROJECT/'brazil_nut/cases'/(name+'.tar'),PROJECT/'mechanism_campaign_20260611/cases'/(name+'.tar')]
    archive=next(p for p in candidates if p.exists())
    with tarfile.open(archive,'r:') as tf:
        members=tf.getmembers()
        print('INDEXED '+str(len(members)),flush=True)
        rec={'case':name,'archive':str(archive),'archive_bytes':archive.stat().st_size,'files':[],'snapshots':[]}
        for m in members:
            if not m.isfile(): continue
            rel=relative(m,name)
            if rel is None: continue
            s=str(rel)
            is_input=(len(rel.parts)==1 and (s.endswith('.json') or s in ('Allrun.sh','parDEMrun.sh','parCFDDEMrun.sh','check_force_balance.py')))
            is_input=is_input or (len(rel.parts)==3 and rel.parts[0]=='CFD' and rel.parts[1] in ('0','constant','system'))
            is_input=is_input or (len(rel.parts)==2 and rel.parts[0]=='DEM' and rel.name.startswith('in.'))
            is_series=s=='analysis/time_series.tsv'
            if is_input or is_series:
                dest=OUT/('cases' if is_input else 'data/cases')/name/(rel if is_input else Path('time_series.tsv'))
                dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes(tf.extractfile(m).read())
                rec['files'].append({'path':str(dest.relative_to(OUT)),'archive_member':m.name,'bytes':m.size})
        if name.startswith('bne_'):
            ts_path=OUT/'data/cases'/name/'time_series.tsv'
            rows=list(csv.DictReader(ts_path.open(),delimiter='\t'))
            if float(rows[-1]['time']) < 40.0-1e-9:
                rows,sources=complete_heights(tf,members)
                with ts_path.open('w',newline='') as output:
                    writer=csv.DictWriter(output,fieldnames=list(rows[0]),delimiter='\t')
                    writer.writeheader();writer.writerows(rows)
                rec['height_derivation']={'method':'Unweighted species arithmetic means from all stored HDF5 particle frames; sorted by time','archive_members':sources}

            times=np.array([float(r['time']) for r in rows]);dz=np.array([float(r['z_gap_large_minus_small']) for r in rows])
            rec['series']={'n_frames':len(rows),'time_end_s':float(times[-1]),'dz_min_m':float(dz.min()),'time_min_s':float(times[dz.argmin()])}
            if 'ang30' in name:
                for t in (float(times[dz.argmin()]),float(times[-1])):
                    rec['snapshots'].append(save_frame(tf,members,name,t))
            if name=='bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30':
                for label,target in [('t0_to_1',.005),('t9_to_10',10.),('t27_to_28',28.),('t39_to_40',40.)]:
                    rec['snapshots'].append(save_frame(tf,members,name,target,label))
                for t in (1.5,28.5):
                    m=next(m for m in members if m.name.endswith(f'/CFD/{t:g}/U'))
                    text=tf.extractfile(m).read().decode()
                    match=re.search(r'internalField\s+nonuniform\s+List<vector>\s*\n(\d+)\s*\n\(',text)
                    body=text[match.end():];body=body[:body.index('\n)\n')]
                    u=np.array([[float(x) for x in v.split()] for v in re.findall(r'\(([-0-9.eE+ ]+)\)',body)])
                    prof=u.reshape(16,32,32,3).mean(axis=(1,2))
                    dest=OUT/'data/cases'/name/f'fluid_profile_{t:g}.csv'
                    with dest.open('w') as f:
                        w=csv.writer(f);w.writerow(['z_m','ux_m_s','uy_m_s','uz_m_s'])
                        w.writerows(np.column_stack(((np.arange(16)+.5)/16*.005,prof)))
                    rec['files'].append({'path':str(dest.relative_to(OUT)),'archive_member':m.name,'derivation':'32x32 horizontal arithmetic mean, 16 vertical layers'})
        provenance.append(rec)
        (OUT/'export_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('DONE '+name,flush=True)
print('EXPORT COMPLETE',flush=True)

# The original generator copies identical CFD/0 fields into every case.
# Restore folders emptied in the archive only after checking that identity.
canonical='bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30'
template=OUT/'cases'/canonical/'CFD/0'
restoration=[]
for case in (OUT/'cases').iterdir():
    existing=list((case/'CFD/0').glob('*'))
    if existing:
        for field in existing:
            if field.name in ('U','Us','rho','p','Ksl','voidfraction'):
                assert field.read_bytes()==(template/field.name).read_bytes(), 'Initial-field templates differ'
for case in (OUT/'cases').iterdir():
    if list((case/'CFD/0').glob('*')): continue
    shutil.copytree(template,case/'CFD/0',dirs_exist_ok=True)
    restoration.append({'case':case.name,'source_case':canonical,'method':'Common initial-field template, as in the original generator','fields':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in template.iterdir()}})
(OUT/'initial_field_restoration.json').write_text(json.dumps(restoration,indent=2)+'\n')

for case in (OUT/'cases').iterdir():
    for folder in ('constant','system'):
        if not list((case/'CFD'/folder).glob('*')):
            assert 'ang30' in case.name, 'Common template applies to the 30-degree cases'
            shutil.copytree(OUT/'cases'/canonical/'CFD'/folder,case/'CFD'/folder,dirs_exist_ok=True)
            for record in restoration:
                if record['case']==case.name:
                    record.setdefault('additional_restored_files',{}).update({str(f.relative_to(case)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (case/'CFD'/folder).iterdir()})
(OUT/'initial_field_restoration.json').write_text(json.dumps(restoration,indent=2)+'\n')
