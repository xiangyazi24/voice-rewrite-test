#!/usr/bin/env python3
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import sympy as sp

OUT = Path('compute/q5161-fast-results.md')

def cooper_terms(N):
    T=[0]*(N+1); T[0]=1
    for m in range(N):
        tm1=T[m-1] if m>=1 else 0
        tm2=T[m-2] if m>=2 else 0
        rhs=2*(2*m+1)*(5*m*m+5*m+2)*T[m]-8*m*(7*m*m+1)*tm1+22*m*(2*m-1)*(m-1)*tm2
        den=(m+1)**3
        assert rhs%den==0
        T[m+1]=rhs//den
    return T

def moments(T,N):
    M=[[Fraction(0) for _ in range(N+1)] for _ in range(4)]
    for n in range(N+1):
        totals=[0,0,0,0]
        for m in range(2*n+1):
            b=comb(2*n,m)*((-2)**(2*n-m))*T[m]
            totals[0]+=b; totals[1]+=m*b; totals[2]+=m*m*b; totals[3]+=m*m*m*b
        den=256**n
        for a in range(4): M[a][n]=Fraction(totals[a],den)
    return M

def A(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n*n+102256019*n+46709052)
def C(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n*n+1587211*n+944620)
def D(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

def q_terms(N):
    q=[Fraction(-215040420000),Fraction(-167282265043404,905),Fraction(-964185327658080,6071)]
    for n in range(2,N):
        q.append(Fraction(B(n),A(n))*q[n]-Fraction(C(n-1),A(n-1))*q[n-1]+Fraction(D(n-2),A(n-2))*q[n-2])
    return q[:N+1]

def rising(x,n):
    y=Fraction(1)
    for k in range(n): y*=x+k
    return y

def gamma_a(a,n): return rising(Fraction(5,2)-a,n)/factorial(n)
def SR(x): return sp.Rational(x.numerator,x.denominator) if isinstance(x,Fraction) else sp.Integer(x)

def fit(Amax=3,J=0,deg=3,pre=False,gshift=False,nrows=31):
    labels=[(a,j,d) for a in range(Amax+1) for j in range(J+1) for d in range(deg+1)]
    rows=[]; rhs=[]
    for n in range(nrows):
        row=[]
        for a,j,d in labels:
            g=gamma_a(a,n+j if gshift else n) if pre else Fraction(1)
            row.append(SR(Fraction(n**d)*g*M[a][n+j]))
        rows.append(row); rhs.append(SR(q[n]))
    X=sp.Matrix(rows); y=sp.Matrix(rhs)
    r=X.rank(); ar=X.row_join(y).rank()
    return r,ar,len(labels)

N=110
T=cooper_terms(2*N+2)
M=moments(T,N)
q=q_terms(N)

cases=[]
cases.append(('direct A=3 J=0 degree=3',)+fit(3,0,3,False,False,31))
cases.append(('Gamma(n) A=3 J=0 degree=3',)+fit(3,0,3,True,False,31))
for J in (1,2):
    cases.append((f'direct shifts A=3 J={J} degree=3',)+fit(3,J,3,False,False,70))
    cases.append((f'Gamma(n) shifts A=3 J={J} degree=3',)+fit(3,J,3,True,False,70))
    cases.append((f'Gamma(n+j) shifts A=3 J={J} degree=3',)+fit(3,J,3,True,True,70))

lines=['# Q5161 fast exact results','','## Exact moments n=0,...,30','','```text','n | M0 | M1 | M2 | M3']
for n in range(31): lines.append(f'{n} | {M[0][n]} | {M[1][n]} | {M[2][n]} | {M[3][n]}')
lines += ['```','','## Rank tests','']
for name,r,ar,u in cases:
    lines.append(f'- **{name}:** rank={r}, augmented rank={ar}, unknowns={u}; '+('CONSISTENT.' if r==ar else 'INCONSISTENT.'))
lines += ['','All values and ranks were computed exactly with Python integers/Fraction and SymPy rational linear algebra.']
OUT.write_text('\n'.join(lines)+'\n')
print(OUT.read_text())
