"""Compare all numeric scientific results; ignore environment-dependent metadata."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

def compare(a,b,where='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),where
        for key in a:compare(a[key],b[key],where+'.'+key)
    elif isinstance(a,list):
        assert len(a)==len(b),where
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f'{where}[{i}]')
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        np.testing.assert_allclose(a,b,rtol=1e-9,atol=1e-10,err_msg=where)
    else:assert a==b,(where,a,b)

p=argparse.ArgumentParser();p.add_argument('first',type=Path);p.add_argument('second',type=Path);args=p.parse_args()
files=[x.relative_to(args.first) for x in args.first.rglob('*') if x.is_file() and x.suffix in ('.json','.csv') and x.name!='run_manifest.json']
other={x.relative_to(args.second) for x in args.second.rglob('*') if x.is_file() and x.suffix in ('.json','.csv') and x.name!='run_manifest.json'}
assert set(files)==other,'Scientific output inventories differ'
for rel in sorted(files):
    a,b=args.first/rel,args.second/rel
    if rel.suffix=='.json':compare(json.loads(a.read_text()),json.loads(b.read_text()),str(rel))
    else:
        x,y=pd.read_csv(a),pd.read_csv(b)
        pd.testing.assert_frame_equal(x,y,check_exact=False,rtol=1e-9,atol=1e-10)
print(f'PASS: {len(files)} scientific files agree at rtol=1e-9, atol=1e-10')
