"""Exact noise covariance from order constraints, not pair kernel formulas."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import json
from verify_threshold import cases, poly_equal

@lru_cache(maxsize=None)
def order_probability(prerequisites):
    n=len(prerequisites)
    full=(1<<n)-1
    @lru_cache(maxsize=None)
    def count(mask):
        if mask==full: return 1
        return sum(count(mask|1<<v) for v in range(n)
                   if not(mask>>v&1) and prerequisites[v]&mask==prerequisites[v])
    return F(count(0),factorial(n))

def joint_probability(adj,u,v,retained):
    n=len(adj)
    def variable(world,vertex):
        return vertex if world==0 or retained>>vertex&1 else n+vertex
    constraints=set()
    variables=set()
    for world,center in [(0,u),(1,v)]:
        a=variable(world,center)
        variables.add(a)
        for w in range(n):
            if adj[center]>>w&1:
                b=variable(world,w)
                variables.add(b)
                constraints.add((a,b))
    index={x:i for i,x in enumerate(sorted(variables))}
    pre=[0]*len(index)
    for a,b in constraints: pre[index[b]] |= 1<<index[a]
    return order_probability(tuple(pre))

def exact_noise(adj):
    n=len(adj)
    by_size=[F(0)]*(n+1)
    for retained in range(1<<n):
        product_moment=sum((joint_probability(adj,u,v,retained)
                            for u in range(n) for v in range(n)),F(0))
        by_size[retained.bit_count()] += product_moment
    result=[F(0)]*(n+1)
    for k,moment in enumerate(by_size):
        for j in range(n-k+1):
            result[k+j] += moment*((-1)**j)*comb(n-k,j)
    mean=sum((F(1,a.bit_count()+1) for a in adj),F(0))
    result[0] -= mean*mean
    return result

def noise_formula(adj):
    n=len(adj)
    degrees={a.bit_count() for a in adj}
    assert len(degrees)==1
    d=degrees.pop()
    assert d>=1
    assert all(not(adj[u]&adj[v]) for u,v in combinations(range(n),2) if adj[u]>>v&1),'triangle'
    moments=[0]*(d+1)
    for u,v in combinations(range(n),2):
        if adj[u]>>v&1: continue
        c=(adj[u]&adj[v]).bit_count()
        for k in range(d+1): moments[k] += comb(c,k)
    result=[F(0)]*(d+2)
    result[2]=F(n*d*(d-1),2*(d+1)**2*(2*d+1))
    for k in range(3,d+2):
        result[k] += F(n*d*comb(d+1,k),(d+1)**2*comb(2*d+1,k))
    for k in range(2,d+1):
        result[k] += F(2*moments[k],(d+1)**2*comb(2*d+1,k))
    return result,moments,d

def extract_histogram(coeff,n,d):
    # Independent algebraic inverse from the measured covariance polynomial.
    moments=[0]*(d+1)
    moments[0]=n*(n-1-d)//2
    moments[1]=n*d*(d-1)//2
    for k in range(2,d+1):
        base=F(n*d*(d-1),2*(d+1)**2*(2*d+1)) if k==2 else F(n*d*comb(d+1,k),(d+1)**2*comb(2*d+1,k))
        moments[k]=F((d+1)**2*comb(2*d+1,k),2)*(coeff[k]-base)
    hist=[F(0)]*(d+1)
    for c in range(d,-1,-1):
        hist[c]=moments[c]-sum(comb(j,c)*hist[j] for j in range(c+1,d+1))
    return hist

if __name__=='__main__':
    selected=['K1,1','C4','C5','C6','C7','C8','K3,3','K4,4','2C4','cube']
    all_cases=cases()
    results=[]
    for name in selected:
        adj=all_cases[name]
        actual=exact_noise(adj)
        predicted,moments,d=noise_formula(adj)
        assert poly_equal(actual,predicted),(name,actual,predicted)
        recovered=extract_histogram(actual,len(adj),d)
        hist=[0]*(d+1)
        for u,v in combinations(range(len(adj)),2):
            if not(adj[u]>>v&1): hist[(adj[u]&adj[v]).bit_count()] += 1
        assert recovered==hist,(name,recovered,hist)
        result={'graph':name,'n':len(adj),'d':d,'covariance_coefficients':[str(x) for x in predicted],
                'S_k':moments,'recovered_m_c':[str(x) for x in recovered],
                'exact_polynomial_match':True,'exact_histogram_recovery':True}
        results.append(result)
        print(name,'OK',result['covariance_coefficients'],flush=True)
    (Path(__file__).resolve().parent/'noise_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print('ALL PASSED',len(results),'unique posets',order_probability.cache_info(),flush=True)
