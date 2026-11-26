"""Профілі: NACA 4-digit і modified 4-digit (00tt-IT), та профіль тіла BODY-150."""
import numpy as np
D1={0.2:0.200,0.3:0.234,0.4:0.315,0.5:0.465,0.6:0.700}
def mod4_half(x,t,I,m):
    """напівтовщина modified NACA 4-digit (Ladson, NASA TM-4741), хорда 1"""
    d0=0.002; d1=D1[m]; sm=1-m
    A=np.array([[sm**2,sm**3],[2*sm,3*sm**2]]); d2,d3=np.linalg.solve(A,[0.1-d0-d1*sm,-d1])
    fpp=2*d2+6*d3*sm                       # f'' по x у точці m
    a0=0.296904*I/6
    M=np.array([[m,m*m,m**3],[1,2*m,3*m*m],[0,2,6*m]])
    rhs=[0.1-a0*np.sqrt(m),-0.5*a0/np.sqrt(m),fpp+0.25*a0*m**-1.5]
    a1,a2,a3=np.linalg.solve(M,rhs)
    x=np.asarray(x,float); y=np.where(x<m,a0*np.sqrt(x)+a1*x+a2*x**2+a3*x**3,
        d0+d1*(1-x)+d2*(1-x)**2+d3*(1-x)**3)
    return y*t/0.2
def naca4_half(x,t):
    x=np.asarray(x,float)
    return 5*t*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
def body_half(x,t=0.15,I=6,m=0.5,hold=(0.50,0.65)):
    """BODY-150: mod4 0.15 I=6 T=0.5, плоска ділянка максимальної товщини від 0.50 до 0.65,
    задня частина стиснута на [0.65,1]"""
    x=np.asarray(x,float); h0,h1=hold
    xa=np.where(x<=h0,x,np.where(x<=h1,m,m+(x-h1)*(1-m)/(1-h1)))
    return mod4_half(xa,t,I,m)
def half(name,x):
    if name=='BODY': return body_half(x)
    if '-' in name:                       # '0016-64'
        a,b=name.split('-'); return mod4_half(x,int(a[2:])/100,int(b[0]),int(b[1])/10)
    return naca4_half(x,int(name[2:])/100)
def coords(name,n=81):
    """замкнений контур: верх від TE до LE, низ від LE до TE (cos-розподіл)"""
    b=np.linspace(0,np.pi,n); x=0.5*(1-np.cos(b)); y=half(name,x)
    xs=np.r_[x[::-1],x[1:]]; ys=np.r_[y[::-1],-y[1:]]
    return xs,ys
if __name__=='__main__':
    print('BODY: товщина (мм на хорді 1 м): x=0.5 %.1f  0.65 %.1f  0.723 %.1f  0.778 %.1f  0.80 %.1f'%tuple(2000*body_half(v) for v in (0.5,0.65,0.723,0.778,0.80)))
    for n in ('0016-64','0012-64','0009','0006','0005'):
        x=np.linspace(0,1,2001); print(n,'tmax %.4f at %.3f'%(2*half(n,x).max(),x[np.argmax(half(n,x))]))
