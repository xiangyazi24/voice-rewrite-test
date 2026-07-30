from fractions import Fraction as F
from math import comb

N=200; INF=10**9
P=lambda n:34*n**3+51*n**2+27*n+5

def prime(n):
    if n<2:return False
    if n%2==0:return n==2
    d=3
    while d*d<=n:
        if n%d==0:return False
        d+=2
    return True

def vi(n,q):
    if n==0:return INF
    n=abs(n); e=0
    while n%q==0:n//=q;e+=1
    return e

def vq(x,q):
    x=F(x)
    return INF if x==0 else vi(x.numerator,q)-vi(x.denominator,q)

def mq(x,q):
    x=F(x); assert x.denominator%q
    return x.numerator%q*pow(x.denominator%q,-1,q)%q

def seq(M):
    a=[F(0)]*(M+1);b=[0]*(M+1);a[0]=0;b[0]=1
    if M:a[1]=6;b[1]=5
    for m in range(1,M):
        z=P(m)*b[m]-m**3*b[m-1];assert z%(m+1)**3==0
        b[m+1]=z//(m+1)**3
        a[m+1]=(P(m)*a[m]-m**3*a[m-1])/(m+1)**3
    return a,b

def harms(M):
    h=[F(0)]*(M+1)
    for i in range(1,M+1):h[i]=h[i-1]+F(1,i)
    return h

def dh(r,h):
    return sum((F(2)*comb(r,k)**2*comb(r+k,k)**2*(h[r+k]-h[r-k]) for k in range(r+1)),F(0))

def mul(u,w):
    return [u[0]*w[0],u[0]*w[1]+u[1]*w[0],u[0]*w[2]+u[1]*w[1]+u[2]*w[0]]

def jets(M):
    J=[[F(0)]*3 for _ in range(M+1)];J[0]=[F(1),F(0),F(0)]
    if M:J[1]=[F(5),F(12),F(0)]
    for r in range(1,M):
        u=mul([F(P(r)),F(102*r*r+102*r+27),F(102*r+51)],J[r])
        w=mul([F(r**3),F(3*r*r),F(3*r)],J[r-1]);z=[u[i]-w[i] for i in range(3)]
        d0=F((r+1)**3);d1=F(3*(r+1)**2);d2=F(3*(r+1))
        y0=z[0]/d0;y1=(z[1]-d1*y0)/d0;y2=(z[2]-d1*y1-d2*y0)/d0
        J[r+1]=[y0,y1,y2]
    return J

def determinant(x,y):return x[0]*y[1]-x[1]*y[0]

def snf(C,q):
    e1=min(vq(x,q) for c in C for x in c if x)
    mm=[vq(determinant(C[i],C[j]),q) for i in range(len(C)) for j in range(i+1,len(C)) if determinant(C[i],C[j])]
    return (e1,None) if not mm else (e1,min(mm)-e1)

def rank(C,q):
    C=[x for x in C if x!=(0,0)]
    if not C:return 0
    return 2 if any((C[i][0]*C[j][1]-C[i][1]*C[j][0])%q for i in range(len(C)) for j in range(i+1,len(C))) else 1

def tail(a,b,u,w):
    return F(b[u]*b[w])*sum((F(6,t**3*b[t-1]*b[t]) for t in range(u+1,w+1)),F(0))

def family(q,r,am,a,b,J,H):
    assert b[r]%q==0
    beta=b[r]//q%q; db=J[r][1]; hh=J[r][2]
    assert db==dh(r,H) and vq(db,q)>=0 and vq(hh,q)>=0
    delta=mq(db,q); out=[]
    for A in range(1,am+1):
        assert A*q+r<=N and 2*A<q and q*q>2*(A*q+r)
        C=[];L=[];ss=[]
        for k in range(A+1):
            m=k*q+r;Y=(F(q**3)*a[m],F(b[m]));assert vq(Y[0],q)>=1 and vq(Y[1],q)>=1
            z=(mq(Y[0]/q,q),mq(Y[1]/q,q));s=(beta+k*delta)%q
            assert z==(s*mq(a[k],q)%q,s*b[k]%q)
            C.append(Y);L.append(z);ss.append(s)
            if k:
                ep=F(q**3)*a[k*q]-a[k];assert vq(ep,q)>=2;tau=ep/q**2;assert vq(tau,q)>=0
                reg=F(b[r])+k*q*db+k*k*q*q*hh
                assert vq(F(b[m])-b[k]*reg,q)>=3
                assert vq(F(q**3)*a[m]-(a[k]*reg+q*q*tau*b[r]),q)>=3
        R=rank(L,q);S=snf(C,q)
        if R==2:assert S==(1,1)
        elif R==1:assert S[0]==1 and S[1]>=2
        else:assert S[0]>=2
        w=[]
        for i in range(A+1):
            for j in range(i+1,A+1):
                u=i*q+r;z=j*q+r;D=determinant(C[i],C[j])
                assert D==q**3*(a[u]*b[z]-a[z]*b[u])
                assert -D/q**3==tail(a,b,u,z);w.append(vq(D,q))
        out.append((A,A*q+r,tuple(ss),R,S,min(w)))
    return (q,r,beta,str(db),delta,str(hh),out)

a,b=seq(N);H=harms(2*N);J=jets(N)
for r in range(N+1):assert J[r][0]==b[r] and J[r][1]==dh(r,H)
rat=[(r,str(J[r][1])) for r in range(1,N+1) if J[r][1].denominator!=1][:12];assert rat
f11=family(11,5,5,a,b,J,H);f17=family(17,3,8,a,b,J,H)
assert all(x[4]==(1,1) for x in f11[-1])
assert f17[-1][0][4]==(1,2) and all(x[4]==(1,1) for x in f17[-1][1:])
fams=[];exc=[];rows=0;rc={0:0,1:0,2:0}
for q in range(7,N+1):
    if not prime(q):continue
    for r in range(1,q-1):
        if b[r]%q:continue
        am=min((q-1)//2,(q*q-1-2*r)//(2*q),(N-r)//q)
        if am<1:continue
        f=family(q,r,am,a,b,J,H);fams.append((q,r,am,f[2],f[4],f[3]))
        for x in f[-1]:rows+=1;rc[x[3]]+=1;exc.append((q,r,f[2],f[4],x)) if x[4]!=(1,1) else None
full=snf([(F(0),F(1)),(F(1),F(0))],11);assert full==(0,0)
print('maximum_index =',N)
print('rational_jet_examples =',rat)
print('q11_r5 =',f11)
print('q17_r3 =',f17)
print('clean_base_families =',fams)
print('clean_base_family_count =',len(fams))
print('simultaneous_family_rows_checked =',rows)
print('first_layer_rank_counts =',rc)
print('non_diag_q_q_rows =',exc)
print('full_moment_map_snf_exponents =',full)
print('formula_failures = endpoint_failures = wedge_failures = 0')
