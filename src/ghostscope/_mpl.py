#!/usr/bin/env python3
"""Concave multinomial profiling on the exact closure of the L2 probability wedge.
Integer active-face candidates avoid treating a numerical optimizer's arbitrary
(g,x,y) as an identified history or mixture proportion.
"""
import json,math,csv
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.special import xlogy

def profile(counts):
    n=np.asarray(counts,dtype=np.int64);N=int(n.sum())
    if N==0:return dict(status='UNESTIMABLE',known_topologies=0)
    den=6*N;candidates=[]
    def add(num,label):
        num=np.array(num,dtype=np.int64)
        if np.all(num>=0) and num.sum()==den and num[0]>=num[2] and num[1]>=num[2]:
            p=num/den;candidates.append((float(xlogy(n,p).sum()),num,label))
    add(6*n,'empirical_interior_or_boundary')
    for i in [0,1]:
        num=6*n.copy();num[i]=num[2]=3*(int(n[i])+int(n[2]));add(num,'pooled_'+['AB','BC'][i]+'_AC')
    add([2*N]*3,'uniform_vertex')
    best,num,face=max(candidates,key=lambda q:q[0]);p=num/den;A=int(num[0]-num[2]);B=int(num[1]-num[2]);lo=B;hi=den-A
    finite=bool(A>0 and B>0 and num[2]>0);g=Fraction(lo+hi,2*den)
    representative=None
    if finite:
        x=1-Fraction(A,den)/(1-g);y=1-Fraction(B,den)/g
        assert 0<g<1 and 0<x<1 and 0<y<1
        representative=dict(gamma=str(g),x=str(x),y=str(y))
        check=[(1-g)*(1-2*x/3)+g*y/3,(1-g)*x/3+g*(1-2*y/3),(1-g)*x/3+g*y/3]
        assert check==[Fraction(int(v),den) for v in num]
    # Pure tree edge P_BC=P_AC with 0<=x<=1; same zero-count convention.
    k=min(3*(int(n[1])+int(n[2])),2*N);pure=np.array([den-2*k,k,k])/den;pure_ll=float(xlogy(n,pure).sum())
    tied=['Inflow','Outflow','Ghost']
    if np.array_equal(num,np.array([den-2*k,k,k])):tied=['Pure Tree']+tied
    attained=['Inflow','Outflow','Ghost'] if finite else []
    if 'Pure Tree' in tied and 0<k<den//3:attained=['Pure Tree']+attained
    return dict(status='FINITE_INTERIOR_MAXIMUM' if finite else 'CLOSURE_SUPREMUM_ONLY',status_domain='strictly positive finite network parameters',known_topologies=N,active_face=face,probabilities=p.tolist(),probability_numerators=num.tolist(),probability_denominator=den,log_pseudolikelihood=best,pure_tree_log_pseudolikelihood=pure_ll,gamma_lower=float(Fraction(lo,den)),gamma_upper=float(Fraction(hi,den)),gamma_lower_numerator=lo,gamma_upper_numerator=hi,gamma_denominator=den,gamma_endpoint_type='open' if finite else 'closure limits; no maximizer in finite strictly positive network domain',tied_histories=tied,tie_domain='equal suprema on the closed probability image',finite_maximizer_histories=attained,one_finite_representative=representative)
