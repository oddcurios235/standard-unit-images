#!/usr/bin/env python3
"""Verify twelve standard-polynomial border cores over two finite fields.
Construct the endpoint and normal map by signed subset and insertion formulas.
A nonzero modular maximal minor certifies nonzero integer determinant.
"""
import json
from pathlib import Path

import numpy as np
BASE = Path(__file__).resolve().parent.parent

def subset_alt(B,p):
 k=len(B);n=B[0].shape[0];F=[None]*(1<<k);F[0]=np.eye(n,dtype=np.int64)
 for mask in range(1,1<<k):
  a=np.zeros((n,n),dtype=np.int64); t=mask
  while t:
   bit=t & -t;i=bit.bit_length()-1;t-=bit
   sign=1 if (mask & (bit-1)).bit_count()%2==0 else -1
   a+=sign*(B[i]@F[mask^bit]);a%=p
  F[mask]=a
 return F

def op_insert(F,k,i,p,u=None):
 """Operator of alternating product with P at index i, B_i excluded if i<k.
 Returns N*N × N*N or N × N*N for result times u.
 """
 n=F[0].shape[0]; active=((1<<k)-1) ^ (1<<i) if i<k else (1<<k)-1
 O=np.zeros((n if u is not None else n*n,n*n),dtype=np.int64)
 sub=active
 while True:
  other=active ^ sub
  inversions=0
  t=sub
  while t:
   bit=t&-t;a=bit.bit_length()-1;t-=bit
   inversions+=(other & ((1<<a)-1)).bit_count()
   if a>i:inversions+=1
  inversions+=(other & ((1<<i)-1)).bit_count()
  sign=1 if inversions%2==0 else -1
  A=F[sub];C=F[other]
  term=np.kron(A,(C@u).reshape(1,n)) if u is not None else np.kron(A,C.T)
  O+=sign*term
  if sub==0:break
  sub=(sub-1)&active
 return O%p

def rank(A,p):
 M=np.array(A,copy=True,dtype=np.int64)%p;m,n=M.shape;r=0
 for j in range(n):
  if r==m:break
  q=np.flatnonzero(M[r:,j]);
  if not len(q):continue
  k=r+int(q[0]);M[[r,k]]=M[[k,r]]
  M[r]=M[r]*pow(int(M[r,j]),-1,p)%p
  nz=np.flatnonzero(M[r+1:,j])+r+1
  if len(nz): M[nz]=(M[nz]-M[nz,j,None]*M[r])%p
  r+=1
 return r

def det(A,p):
 M=np.array(A,dtype=np.int64,copy=True)%p;n=M.shape[0];v=1
 for j in range(n):
  q=np.flatnonzero(M[j:,j]);
  if not len(q):return 0
  k=j+int(q[0]);
  if k!=j:M[[j,k]]=M[[k,j]];v=-v
  d=int(M[j,j]);v=v*d%p;inv=pow(d,-1,p)
  nz=np.flatnonzero(M[j+1:,j])+j+1
  if len(nz):M[nz]=(M[nz]-M[nz,j,None]*inv*M[j])%p
 return int(v%p)

def pivot_rows(A,p):
 """Independent row indices giving a nonsingular maximal minor of A mod p."""
 M=np.array(A.T,copy=True,dtype=np.int64)%p
 row_count,column_count=M.shape
 pivots=[];r=0
 for j in range(column_count):
  if r==row_count: break
  candidates=np.flatnonzero(M[r:,j])
  if not len(candidates): continue
  t=r+int(candidates[0]);M[[r,t]]=M[[t,r]]
  M[r]=M[r]*pow(int(M[r,j]),-1,p)%p
  nz=np.flatnonzero(M[r+1:,j])+r+1
  if len(nz): M[nz]=(M[nz]-M[nz,j,None]*M[r])%p
  pivots.append(j);r+=1
 return pivots

def check(case,p):
 k=2*case['r']-1;N=case['N'];B=[np.array(x,dtype=np.int64)%p for x in case['B']];u=np.array(case['u'],dtype=np.int64)%p
 F=subset_alt(B,p);S=F[-1]
 G=np.concatenate([op_insert(F,k,k,p),* [op_insert(F,k,i,p,u) for i in range(k)]],axis=0)%p
 rows=pivot_rows(G,p)
 minor=det(G[rows,:],p) if len(rows)==N*N else 0
 v={'r':case['r'],'N':N,'p':p,'endpoint_det':det(S,p),'normal_rank':rank(G,p),'normal_target':N*N,'minor_rows':rows,'minor_det':minor,'seed_col':(S@u%p).tolist()}
 return v

REQUIRED_CASES = {(r, N) for r in (3, 4, 5) for N in range(r, 2*r)}


def main():
    inputs = json.loads((BASE / 'certificates/standard_cores.json').read_text())['cases']
    pairs = [(c['r'], c['N']) for c in inputs]
    if len(pairs) != len(REQUIRED_CASES) or set(pairs) != REQUIRED_CASES:
        raise RuntimeError('Core inputs must cover each of the twelve prescribed (r,N) pairs exactly once')

    out = []
    for c in inputs:
        for prime in (101, 103):
            result = check(c, prime)
            out.append(result)
            print(result, flush=True)
    if len(out) != 24 or not all(
        v['endpoint_det'] and v['minor_det']
        and v['normal_rank'] == v['normal_target']
        and any(v['seed_col']) for v in out
    ):
        raise RuntimeError('One or more standard-core certificates failed')

    output = {'method': 'Signed-subset alternating products and signed insertion normal maps; nonzero modular maximal minors', 'cases': out}
    (BASE / 'certificates' / 'standard_core_checks.json').write_text(json.dumps(output, indent=2) + '\n')
    print('PASS: 12 core pairs; 24 full-rank normal maps and nonzero minors (p=101,103).')


if __name__ == '__main__':
    main()
