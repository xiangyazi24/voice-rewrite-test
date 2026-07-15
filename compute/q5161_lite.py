#!/usr/bin/env python3
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

OUT=Path('compute/q5161-lite-results.md')

def Tseq(N):
    T=[0]*(N+1); T[0]=1
    for m in range(N):
        a=T[m-1] if m>=1 else 0; b=T[m-2] if m>=2 else 0
        num=2*(2*m+1)*(5*m*m+5*m+2)*T[m]-8*m*(7*m*m+1)*a+22*m*(2*m-1)*(m-1)*b
        den=(m+1)**3; assert num%den==0; T[m+1]=num//den
    return T

def moments(T,N):
    M=[[Fraction(0) for _ in range(N+1)] for _ in range(4)]
    for n in range(N+1):
        s=[0,0,0,0]
        for m in range(2*n+1):
            z=comb(2*n,m)*(-2)**(2*n-m)*T[m]
            s[0]+=z; s[1]+=m*z; s[2]+=m*m*z; s[3]+=m*m*m*z
        for a in range(4): M[a][n]=Fraction(s[a],256**n)
    return M

def A(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n*n+102256019*n+46709052)
def C(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n*n+1587211*n+944620)
def D(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)
def qseq(N):
    q=[Fraction(-215040420000),Fraction(-167282265043404,905),Fraction(-964185327658080,6071)]
    for n in range(2,N): q.append(Fraction(B(n),A(n))*q[n]-Fraction(C(n-1),A(n-1))*q[n-1]+Fraction(D(n-2),A(n-2))*q[n-2])
    return q[:N+1]
def rising(x,n):
    y=Fraction(1)
    for k in range(n): y*=x+k
    return y
def gamma(a,n): return rising(Fraction(5,2)-a,n)/factorial(n)

def rank(mat):
    a=[list(map(Fraction,row)) for row in mat]
    if not a: return 0
    m=len(a); n=len(a[0]); r=0
    for c in range(n):
        p=next((i for i in range(r,m) if a[i][c]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        v=a[r][c]; a[r]=[x/v for x in a[r]]
        for i in range(m):
            if i!=r and a[i][c]:
                v=a[i][c]; a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==m: break
    return r

def test(pre=False):
    X=[]; y=[]
    for n in range(31):
        row=[]
        for a in range(4):
            g=gamma(a,n) if pre else Fraction(1)
            for d in range(4): row.append(Fraction(n**d)*g*M[a][n])
        X.append(row); y.append(q[n])
    r=rank(X); ar=rank([row+[yy] for row,yy in zip(X,y)])
    return r,ar

T=Tseq(62); M=moments(T,30); q=qseq(30)
r0,a0=test(False); r1,a1=test(True)
lines=['# Q5161 lightweight exact result','','```text','n | M0 | M1 | M2 | M3']
for n in range(31): lines.append(f'{n} | {M[0][n]} | {M[1][n]} | {M[2][n]} | {M[3][n]}')
lines += ['```','',f'- Direct degree<=3: rank={r0}, augmented rank={a0}, unknowns=16; '+('consistent' if r0==a0 else 'INCONSISTENT'),f'- Gamma degree<=3: rank={r1}, augmented rank={a1}, unknowns=16; '+('consistent' if r1==a1 else 'INCONSISTENT'),'','All arithmetic is exact (`int` and `Fraction`).']
OUT.write_text('\n'.join(lines)+'\n')
print(OUT.read_text())
