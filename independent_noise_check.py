"""Independent exact check using linear extensions of two minimum-event posets."""
from fractions import Fraction as F
from math import comb, factorial
from itertools import combinations
from functools import lru_cache
import json

def edges_graph(n, edges):
    g = [set() for _ in range(n)]
    for u,v in edges:
        g[u].add(v); g[v].add(u)
    return g

def event_probability(g,u,v,R):
    # First-star labels are distinct latent variables. The second-star label is
    # shared with the first iff its coordinate belongs to R.
    a = sorted(g[u] | {u}); b = sorted(g[v] | {v})
    ids = {(0,w):i for i,w in enumerate(a)}
    size=len(a)
    for w in b:
        if w in R: ids[1,w]=ids[0,w]
        else: ids[1,w]=size; size+=1
    prereq=[0]*size
    for w in g[u]: prereq[ids[0,w]] |= 1<<ids[0,u]
    for w in g[v]: prereq[ids[1,w]] |= 1<<ids[1,v]
    @lru_cache(None)
    def dp(mask):
        if mask==(1<<size)-1: return 1
        return sum(dp(mask | (1<<i)) for i in range(size)
                   if not mask>>i&1 and prereq[i]&mask==prereq[i])
    return F(dp(0),factorial(size))

def direct_polynomial(g):
    n=len(g); d=len(g[0]); result=[F(0) for _ in range(d+2)]
    for u in range(n):
        for v in range(n):
            overlap=sorted((g[u]|{u}) & (g[v]|{v})); t=len(overlap)
            for j in range(t+1):
                for retained in combinations(overlap,j):
                    q=event_probability(g,u,v,set(retained))
                    for k in range(t-j+1):
                        result[j+k] += q*comb(t-j,k)*(-1)**k
    result[0] -= F(n*n,(d+1)**2)
    return result

def proposed(g):
    n=len(g); d=len(g[0]); S=[0]*(d+1)
    for u in range(n):
        for v in range(u+1,n):
            if v not in g[u]:
                c=len(g[u]&g[v])
                for k in range(d+1): S[k]+=comb(c,k)
    p=[F(0)]*(d+2)
    if d>=1: p[2]=F(n*d*(d-1),2*(d+1)**2*(2*d+1))
    for k in range(3,d+2):
        p[k]+=F(n*d*comb(d+1,k),(d+1)**2*comb(2*d+1,k))
    for k in range(2,d+1): p[k]+=F(2*S[k],(d+1)**2*comb(2*d+1,k))
    return p

graphs={}
for n in [4,5,6,7]:
    graphs['C'+str(n)]=edges_graph(n,[(i,(i+1)%n) for i in range(n)])
graphs['K2']=edges_graph(2,[(0,1)])
for d in [3,4]:
    graphs['K'+str(d)+','+str(d)]=edges_graph(2*d,[(i,j) for i in range(d) for j in range(d,2*d)])
graphs['Q3']=edges_graph(8,[(i,i^(1<<j)) for i in range(8) for j in range(3) if i<(i^(1<<j))])
graphs['Petersen']=edges_graph(10,[(i,(i+1)%5) for i in range(5)]+[(i,5+i) for i in range(5)]+[(5+i,5+(i+2)%5) for i in range(5)])
records=[]
for name,g in graphs.items():
    exact=direct_polynomial(g); pred=proposed(g)
    record={'graph':name,'pass':exact==pred,'exact_coefficients':[str(x) for x in exact]}
    records.append(record); print(json.dumps(record),flush=True)
assert all(x['pass'] for x in records)
