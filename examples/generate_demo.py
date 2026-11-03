"""Reproduce fictional sequences and hand-authored annotations; no real organism."""
import json
import random
from pathlib import Path

rng=random.Random(42)
def dna(n,weights):
    return ''.join(rng.choices('ACGT',weights=weights,k=n))
a=dna(30000,[2,2,4,2])+dna(30000,[2,4,2,2])+dna(30000,[3,2,2,3])
a=a[:45000]+'N'*300+a[45300:]
b=dna(14000,[4,1,1,4])
sequences=[{'id':'demo_chromosome','nt':a,'length':len(a),'topology':'linear'},
           {'id':'demo_contig_2','nt':b,'length':len(b),'topology':'linear'}]
features=[]
for sequence in sequences:
    for i,start in enumerate(range(500,sequence['length']-1000,1800)):
        kind='tRNA' if i%9==0 else 'CDS'
        features.append({'sequence':sequence['id'],'type':kind,'start':start,'stop':start+(75 if kind=='tRNA' else 850),
                         'strand':'+' if i%2 else '-', 'locus':f"DEMO_{sequence['id']}_{i:04d}",
                         'gene':'', 'product':'illustrative tRNA' if kind=='tRNA' else ('hypothetical protein' if i%3==0 else 'illustrative protein feature')})
payload={'version':{'baklens_fixture':'synthetic; not produced by Bakta'},'run':{'description':'Deterministic teaching fixture, seed 42'},
         'genome':{'name':'Fictional demonstration'},'sequences':sequences,'features':features}
root=Path(__file__).parent
(root/'demo.bakta.json').write_text(json.dumps(payload,indent=2))
(root/'demo.fna').write_text(''.join(f">{s['id']} synthetic demonstration\n"+'\n'.join(s['nt'][i:i+80] for i in range(0,len(s['nt']),80))+'\n' for s in sequences))
