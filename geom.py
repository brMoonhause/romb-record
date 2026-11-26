"""Геометрія РОМБ D5 як набори точок (без CAD): спільна для build123d і FreeCAD."""
import numpy as np
from params import *
from foils import coords
def smooth(t): t=np.clip(t,0,1); return t*t*(3-2*t)
# ---------- перерізи крила ----------
def full_sections():
    return [(-y,f,c,o) for y,f,c,o in SECTIONS[::-1] if y>0]+SECTIONS
def section_pts(y,foil,c,off,n=None):
    xs,zs=coords(foil,n or N_PTS); return np.c_[off+c*xs, np.full_like(xs,y), c*zs]
_ST=[(y,section_pts(y,f,c,o)) for y,f,c,o in SECTIONS]
def z_top(x,y):
    """верхня поверхня (лінійчата між станціями, як і loft); 0 поза планформою"""
    y=abs(y)
    for (ya,A),(yb,B) in zip(_ST[:-1],_ST[1:]):
        if ya<=y<=yb:
            t=(y-ya)/(yb-ya); n=(len(A)+1)//2
            P=A[:n]*(1-t)+B[:n]*t; xs,zs=P[::-1,0],P[::-1,2]
            if x<xs[0] or x>xs[-1]: return 0.0
            return float(np.interp(x,xs,zs))
    return 0.0
def le_x(y):
    y=abs(y); ys=[s[0] for s in SECTIONS]; xs=[s[3] for s in SECTIONS]; return float(np.interp(y,ys,xs))
def te_x(y):
    y=abs(y); ys=[s[0] for s in SECTIONS]; xs=[s[3]+s[2] for s in SECTIONS]; return float(np.interp(y,ys,xs))
def body_le_y(x):
    """край планформи (передня кромка) як функція x — для розміщення губи"""
    pts=[(s[3],s[0]) for s in SECTIONS]; return float(np.interp(x,[p[0] for p in pts],[p[1] for p in pts]))
# ---------- NACA 1-series (NASA TM X-3324, табл. I; пропущені рядки інтерпольовані) ----------
_N1=np.array([(0,0),(0.2,4.80),(0.4,6.63),(0.6,8.12),(0.8,9.33),(1,10.38),(1.5,12.72),(2,14.72),
 (4.5,22.96),(5,24.36),(6,27.01),(7,29.47),(8,31.81),(9,34.03),(10,36.13),(11,38.15),(12,40.09),(13,41.94),
 (20,52.70),(21,54.05),(22,55.37),(23,56.66),(24,57.92),(25,59.15),(26,60.35),(27,61.52),(28,62.67),(29,63.79),
 (30,64.89),(31,65.97),(32,67.03),(33,68.07),(34,69.08),(35,70.08),(36,71.05),(37,72.00),(38,72.94),(39,73.85),
 (40,74.75),(41,75.63),(48,81.25),(49,81.99),(50,82.69),(52,84.10),(54,85.45),(56,86.73),(58,87.95),(60,89.11),
 (62,90.20),(64,91.23),(66,92.20),(68,93.11),(70,93.95),(72,94.75),(74,95.48),(76,96.16),(78,96.79),(80,97.35),
 (82,97.87),(84,98.33),(86,98.74),(88,99.09),(90,99.40),(100,100.0)])/100
def naca1(u):
    """y/Y як функція x/X (0..1); інтерполяція в координатах sqrt для гладкості біля носика"""
    u=np.clip(u,0,1); return np.interp(np.sqrt(u),np.sqrt(_N1[:,0]),_N1[:,1])
# ---------- контури ----------
def _resample(P,n):
    Q=np.r_[P,P[:1]]; d=np.r_[0,np.cumsum(np.hypot(*np.diff(Q,axis=0).T))]
    s=np.linspace(0,d[-1],n,endpoint=False)
    return np.c_[np.interp(s,d,Q[:,0]),np.interp(s,d,Q[:,1])]
def rrect(w,h,r,n=72):
    """закруглений прямокутник з центром (0,0), старт (w/2,0), проти год. стрілки"""
    r=min(r,w/2-1e-3,h/2-1e-3); hw,hh=w/2-r,h/2-r; P=[]
    for cx,cz,a0 in ((hw,hh,0),(-hw,hh,np.pi/2),(-hw,-hh,np.pi),(hw,-hh,1.5*np.pi)):
        a=np.linspace(a0,a0+np.pi/2,12); P.append(np.c_[cx+r*np.cos(a),cz+r*np.sin(a)])
    P=np.vstack(P); P=_resample(P,4*n)
    i=np.argmin(np.abs(P[:,1])+(P[:,0]<0)*1e3); return _resample(np.roll(P,-i,axis=0),n)
def dshape(r,split,n=72):
    """половина кола (опукла в +y), плоска сторона на y=split; старт (split+r,0)"""
    a=np.linspace(-np.pi/2,np.pi/2,60); arc=np.c_[split+r*np.cos(a),r*np.sin(a)]
    flat=np.linspace([split,r],[split,-r],30)[1:-1]; P=np.r_[arc,flat]
    i=np.argmin(np.abs(P[:,1])+(P[:,0]<split+r/2)*1e3); return _resample(np.roll(P,-i,axis=0),n)
def area(P): y,z=P[:,0],P[:,1]; return 0.5*abs(np.dot(y,np.roll(z,1))-np.dot(z,np.roll(y,1)))
def rr_dims(A,ar,rf=0.45):
    """w,h,r закругленого прямокутника площею A, h/w=ar, r=rf*w"""
    w=np.sqrt(A/(ar-(4-np.pi)*rf*rf)); return w,ar*w,rf*w
# ---------- повітрозабірник і канал ----------
W_HL,H_HL,R_HL=rr_dims(A_THROAT*INLET_CR,INLET_AR)
W_T,H_T,R_T=rr_dims(A_THROAT,INLET_AR)
Y_INLET_IN=body_le_y(X_INLET)+LIP_T+LIP_GAP          # внутрішня стінка каналу біля губи
YC0=Y_INLET_IN+W_HL/2
def duct_center(x):
    t=smooth((x-X_INLET-L_INLIP)/(X_DUCT2-X_INLET-L_INLIP)); return (1-t)*YC0+t*DUCT_Y_END, t
def duct_area_at(x):
    if x<=X_INLET: return A_THROAT*INLET_CR
    if x<=X_INLET+L_INLIP:                            # еліптична внутрішня губа
        u=(x-X_INLET)/L_INLIP; k=1-np.sqrt(1-(1-u)**2)
        return A_THROAT*(1+(INLET_CR-1)*k)
    u=(x-X_INLET-L_INLIP)/(X_DUCT2-X_INLET-L_INLIP)     # дифузор: площа росте лінійно
    return A_THROAT+(A_DUCT_END-A_THROAT)*u
def duct_shape(x):
    yc,t=duct_center(x); A=duct_area_at(x); ar=INLET_AR*(1-t)+t*1.0
    w,h,r=rr_dims(A,ar); return yc,t,w,h,r
def duct_section(x,side=1):
    yc,t,w,h,r=duct_shape(x); P=rrect(w,h,r)+[yc,0]
    return np.c_[np.full(len(P),x),side*P[:,0],P[:,1]]
def duct_sections(side=1):
    xs=np.r_[X_INLET-8, np.linspace(X_INLET,X_DUCT2,N_DUCT)]
    return [duct_section(x,side) for x in xs]
def cheek_sections(side=1):
    """зовнішня щока: губа NACA 1-series (товщина LIP_T на довжині LIP_L), далі — стінка каналу"""
    out=[]
    for x in np.linspace(X_INLET,X_CHEEK_END,N_CHEEK):
        yc,t,w,h,r=duct_shape(x)
        if x<=X_INLET+LIP_L:
            d=LIP_T*naca1((x-X_INLET)/LIP_L); w,h,r=W_HL+2*d,H_HL+2*d,R_HL+d
        else:
            w,h,r=max(w+2*LIP_T,W_HL+2*LIP_T),max(h+2*LIP_T,H_HL+2*LIP_T),r+LIP_T
        P=rrect(w,h,r)+[yc,0]
        out.append(np.c_[np.full(len(P),x),side*P[:,0],P[:,1]])
    return out
# ---------- двигун і відсік ----------
def engine_profile():
    return [(X_NOZZLE+(xl-9.95)*10,r) for xl,r in ENGINE]
def bay_profile():
    """камера перед двигуном (обидва канали) + огинаюча двигуна з зазором + канал вихлопу"""
    E=engine_profile(); rmax=max(r for x,r in E); gp=BAY_GAP
    xs=[x for x,r in E if r>=rmax-1]; x_case0,x_case1=xs[0],xs[-1]
    x_noz,r_noz=E[-2]
    dome=[(X_CH0,0.5),(X_CH0,R_CH*0.45),(X_CH0+6,R_CH*0.75),(X_CH0+16,R_CH*0.93),(X_CH0+30,R_CH)]
    return dome+[(x_case0,max(R_CH,rmax+gp)),(x_case1+5,rmax+gp),(x_noz,r_noz+gp),(X_EXHAUST,r_noz+gp)]
if __name__=='__main__':
    print(f'губа: {W_HL:.1f}×{H_HL:.1f} мм (r{R_HL:.1f}), горло {W_T:.1f}×{H_T:.1f}; y губи {Y_INLET_IN:.1f}–{Y_INLET_IN+W_HL:.1f}')
    for x in np.r_[X_INLET,X_INLET+L_INLIP,np.linspace(X_INLET+40,X_DUCT2,10)]:
        s=duct_section(x); A=area(s[:,1:])
        gaps=[z_top(x,y)-abs(z) for y,z in s[:,1:] if z_top(x,y)>0]
        print(f'x={x:6.1f} A={A/100:5.1f} см² y {s[:,1].min():5.1f}..{s[:,1].max():5.1f} мін.стінка {min(gaps) if gaps else float("nan"):6.1f}')
    print('відсік:',[(round(a),round(b,1)) for a,b in bay_profile()])
