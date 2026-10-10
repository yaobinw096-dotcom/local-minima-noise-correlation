"""Exact independent checks via permutation-prefix dynamic programming."""
from fractions import Fraction as F
from math import comb, factorial
from itertools import combinations
import json
from pathlib import Path

def graph(n, edges):
    adj = [0]*n
    for u,v in edges:
        adj[u] |= 1<<v
        adj[v] |= 1<<u
    return adj

def cycle(n):
    return graph(n, [(i,(i+1)%n) for i in range(n)])

def complete(n):
    return graph(n, combinations(range(n),2))

def bipartite(d):
    return graph(2*d, [(i,d+j) for i in range(d) for j in range(d)])

def union(a,b):
    return a + [x<<len(a) for x in b]

def poly_add(a,b,scale=F(1)):
    out = a[:] + [F(0)]*max(0,len(b)-len(a))
    for i,x in enumerate(b): out[i] += scale*x
    return out

def poly_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j] += x*y
    return out

def poly_equal(a,b):
    n=max(len(a),len(b))
    return a+[F(0)]*(n-len(a)) == b+[F(0)]*(n-len(b))

def prefix_polynomials(adj):
    """Count all orders of every subset, independent of codegree formulas."""
    n=len(adj)
    total=1<<n
    first=[0]*total
    second=[0]*total
    first_by_size=[0]*(n+1)
    second_by_size=[0]*(n+1)
    for mask in range(total):
        m=mask.bit_count()
        first_by_size[m] += first[mask]
        second_by_size[m] += second[mask]
        for v in range(n):
            if mask>>v&1: continue
            inc=int(not(adj[v]&mask))
            nxt=mask|(1<<v)
            first[nxt] += first[mask]+inc*factorial(m)
            second[nxt] += second[mask]+2*inc*first[mask]+inc*factorial(m)
    # t=1-s: P(exactly m labels <=t)=C(n,m)(1-s)^m s^(n-m).
    mean=[F(0)]*(n+1)
    moment2=[F(0)]*(n+1)
    for m in range(n+1):
        for j in range(m+1):
            coeff=F((-1)**j*comb(m,j),factorial(m))
            mean[n-m+j] += coeff*first_by_size[m]
            moment2[n-m+j] += coeff*second_by_size[m]
    variance=poly_add(moment2,poly_mul(mean,mean),F(-1))
    return mean,variance

def threshold_formula(adj):
    n=len(adj)
    degrees={x.bit_count() for x in adj}
    assert len(degrees)==1
    d=degrees.pop()
    p=[F(0)]*(d+2)
    p[0]=F(1,d+1)
    p[d+1]=-p[0]
    ms=[0]*(d+1)
    for u,v in combinations(range(n),2):
        if not(adj[u]>>v&1): ms[(adj[u]&adj[v]).bit_count()] += 1
    variance=poly_add([x*n for x in p],poly_mul(p,p),F(-n*n))
    for c,m in enumerate(ms):
        k=2*d+2-c
        q=poly_add(p,[F(1,k)]+[F(0)]*(k-1)+[F(-1,k)],F(-1))
        variance=poly_add(variance,q,F(4*m,d-c+1))
    return p,variance,ms,d

def cases():
    out={f'C{n}':cycle(n) for n in range(3,13)}
    out.update({f'K{n}':complete(n) for n in range(2,9)})
    out.update({f'K{d},{d}':bipartite(d) for d in range(1,7)})
    out['2C3']=union(cycle(3),cycle(3))
    out['2C4']=union(cycle(4),cycle(4))
    out['C3+C5']=union(cycle(3),cycle(5))
    out['cube']=graph(8,[(u,u^(1<<b)) for u in range(8) for b in range(3) if u<(u^(1<<b))])
    out['prism']=graph(6,[(i,(i+1)%3) for i in range(3)]+[(i+3,3+(i+1)%3) for i in range(3)]+[(i,i+3) for i in range(3)])
    pairs=list(combinations(range(5),2))
    out['Petersen']=graph(10,[(i,j) for i,j in combinations(range(10),2) if not(set(pairs[i])&set(pairs[j]))])
    out['Heawood']=graph(14,[(i,7+(i+s)%7) for i in range(7) for s in (0,1,3)])
    return out

if __name__=='__main__':
    results=[]
    for name,adj in cases().items():
        mean,var=prefix_polynomials(adj)
        p,pred,ms,d=threshold_formula(adj)
        assert poly_equal(mean,[x*len(adj) for x in p]),name+' mean'
        assert poly_equal(var,pred),name+' variance'
        for c in range(1,d+1):
            k=2*d+2-c
            actual=var[k] if k<len(var) else F(0)
            assert actual==F(4*ms[c],(d-c+1)*k),(name,c,actual)
        result={'graph':name,'n':len(adj),'d':d,'m_c':ms,'variance_at_t1':str(var[0]),'exact_polynomial_match':True}
        results.append(result)
        print(name,'OK',result['variance_at_t1'])
    (Path(__file__).resolve().parent/'threshold_results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    print('ALL PASSED',len(results))
