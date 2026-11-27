"""Геометрія РОМБ D5 як набори точок (без CAD): спільна для build123d і FreeCAD."""
import numpy as np
from params import *
from foils import coords, half
def smooth(t): t=np.clip(t,0,1); return t*t*(3-2*t)
# ---------- гладка поверхня крила ----------
# передня кромка — ламана через станції з заокругленими зламами; задня — пряма; товщина між станціями —
# монотонний кубічний сплайн (PCHIP) по y для кожної точки хорди. Поверхня — гладкий loft через N_SPAN перерізів.
_XC=0.5*(1-np.cos(np.linspace(0,np.pi,N_PTS)))
_SY=np.array([s[0] for s in SECTIONS])
_SH=np.array([half(s[1],_XC)*1.0 for s in SECTIONS])      # напівтовщина / хорда в точках _XC
def _pchip_slopes(x,y):
    h=np.diff(x); d=np.diff(y,axis=0)/h[:,None]; m=np.zeros_like(y)
    for k in range(1,len(x)-1):
        w1=2*h[k]+h[k-1]; w2=h[k]+2*h[k-1]
        same=(d[k-1]*d[k])>0
        m[k]=np.where(same,(w1+w2)/(w1/np.where(same,d[k-1],1)+w2/np.where(same,d[k],1)),0)
    m[0]=d[0]; m[-1]=d[-1]
    return m
_SM=_pchip_slopes(_SY,_SH)
def _hcol(y):
    """напівтовщина/хорда в точках _XC на розмаху y (PCHIP)"""
    y=min(abs(y),_SY[-1]); k=int(np.clip(np.searchsorted(_SY,y)-1,0,len(_SY)-2))
    h=_SY[k+1]-_SY[k]; t=(y-_SY[k])/h
    h00=2*t**3-3*t**2+1; h10=t**3-2*t**2+t; h01=-2*t**3+3*t**2; h11=t**3-t**2
    return h00*_SH[k]+h10*h*_SM[k]+h01*_SH[k+1]+h11*h*_SM[k+1]
def le_x(y):
    """передня кромка: ламана через станції, злами заокруглені квадратичною кривою Безьє"""
    y=abs(y); P=[(s[0],s[3]) for s in SECTIONS]
    for i in range(1,len(P)-1):
        yc,xc=P[i]; b=LE_BLEND.get(yc,0)
        if b>0 and abs(y-yc)<b:
            (y0,x0),(y2,x2)=P[i-1],P[i+1]
            ya,xa=yc-b,xc+(x0-xc)*b/(yc-y0); yb,xb=yc+b,xc+(x2-xc)*b/(y2-yc)
            # параметр Безьє t по y (квадратне рівняння)
            A=ya-2*yc+yb; B=2*(yc-ya); C=ya-y
            t=(-C/B) if abs(A)<1e-9 else (-B+np.sqrt(max(B*B-4*A*C,0)))/(2*A)
            return float((1-t)**2*xa+2*(1-t)*t*xc+t*t*xb)
    return float(np.interp(y,[p[0] for p in P],[p[1] for p in P]))
def te_x(y): return TE_ROOT-abs(y)*np.tan(np.radians(TE_FWD_SWEEP))
def chord(y): return te_x(y)-le_x(y)
def span_stations():
    u=np.linspace(0,1,N_SPAN); ys=SPAN_HALF*(1-np.cos(u*np.pi/2))**1.0   # густіше біля кореня
    return ys
def full_sections():
    ys=span_stations(); return [-y for y in ys[::-1] if y>0]+list(ys)
def section_pts(y,*_):
    c=chord(y); le=le_x(y); h=_hcol(y)
    xs=np.r_[_XC[::-1],_XC[1:]]; zs=np.r_[h[::-1],-h[1:]]
    return np.c_[le+c*xs, np.full_like(xs,y), c*zs]
def z_top(x,y):
    y=abs(y)
    if y>SPAN_HALF: return 0.0
    le=le_x(y); c=chord(y); u=(x-le)/c
    if u<0 or u>1: return 0.0
    return float(c*np.interp(u,_XC,_hcol(y)))
def body_le_y(x):
    ys=np.linspace(0,SPAN_HALF,4001); xs=np.array([le_x(y) for y in ys])
    return float(np.interp(x,xs,ys))
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
# ---------- повітрозабірник у товщині крила і канал ----------
def _lens(x,y0,y1,d,n=72):
    """переріз крила площиною x між y0 і y1, зміщений усередину на d (d<0 — назовні); старт — зовнішній кінчик, проти год. стрілки"""
    ys=np.linspace(y0,y1,240); h=np.array([z_top(x,y)-d for y in ys])
    ok=h>0.3; ys,h=ys[ok],h[ok]
    top=np.c_[ys[::-1],h[::-1]]; bot=np.c_[ys,-h]
    P=np.r_[bot,top]                       # від внутр-низу → зовн. кінчик → верх → внутр-верх
    P=_resample(P,4*n); i=np.argmax(P[:,0]-1e-3*np.abs(P[:,1]))
    return _resample(np.roll(P,-i,axis=0),n)
def _inlet_y():
    """внутрішній край входу підбирається так, щоб площа входу = A_THROAT*INLET_CR"""
    yout=body_le_y(X_INLET)-INLET_EDGE; lo,hi=0.0,yout-5
    for _ in range(40):
        yin=(lo+hi)/2; A=area(_lens(X_INLET,yin+LIP_T,yout,LIP_T))
        if A<A_THROAT*INLET_CR: hi=yin
        else: lo=yin
    return (lo+hi)/2,yout
Y_INLET_IN,Y_INLET_OUT=_inlet_y()
MOUTH=_lens(X_INLET,Y_INLET_IN+LIP_T,Y_INLET_OUT,LIP_T)
YC0=float(np.mean(MOUTH[:,0])); A_MOUTH=area(MOUTH)
def notch_box():
    """усе, що перед заокругленням губи в смузі входу, прибирається (виїмка в кромці)"""
    return (le_x(Y_INLET_IN)-30, X_INLET-LIP_R+2, Y_INLET_IN-1.0, Y_INLET_OUT+40)
def flare_sections(side=1,n=7):
    out=[]
    for u in np.linspace(0,1,n):
        x=X_INLET-LIP_R*(1-u); d=-1.5+(LIP_T+1.5)*np.sqrt(1-(1-u)**2)
        P=_lens(x,Y_INLET_IN+d,Y_INLET_OUT+(6 if d<0 else 0),d)   # внутрішня губа теж заокруглена
        out.append(np.c_[np.full(len(P),x),side*P[:,0],P[:,1]])
    return out
def duct_center(x):
    t=smooth((x-X_INLET-L_INLIP)/(X_DUCT2-X_INLET-L_INLIP)); return (1-t)*YC0+t*DUCT_Y_END, t
def duct_area_at(x):
    if x<=X_INLET: return A_MOUTH
    if x<=X_INLET+L_INLIP:
        u=(x-X_INLET)/L_INLIP; k=1-np.sqrt(1-(1-u)**2)
        return A_THROAT+(A_MOUTH-A_THROAT)*k
    u=(x-X_INLET-L_INLIP)/(X_DUCT2-X_INLET-L_INLIP)
    return A_THROAT+(A_DUCT_END-A_THROAT)*smooth(u)      # плавне розширення: нульовий нахил на кінцях
_M0=MOUTH-[YC0,0]
def duct_section(x,side=1):
    yc,t=duct_center(x); A=duct_area_at(x)
    w,h,r=rr_dims(A_MOUTH,1.0); C=rrect(w,h,r)
    P=(1-t)*_M0+t*C; P=P*np.sqrt(A/area(P))+[yc,0]
    return np.c_[np.full(len(P),x),side*P[:,0],P[:,1]]
def duct_sections(side=1):
    # перший переріз трохи менший за кінець заокруглення губи (сходинка ≈0,7 мм) — щоб поверхні не збігались
    first=_M0*0.975+[YC0,0]; first=np.c_[np.full(len(first),X_INLET-1.0),side*first[:,0],first[:,1]]
    xs=np.linspace(X_INLET+4,X_DUCT2,N_DUCT)
    return [first]+[duct_section(x,side) for x in xs]
def duct_shape(x):
    s=duct_section(x); return s[:,1].mean(),0,s[:,1].max()-s[:,1].min(),2*s[:,2].max(),0
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
    print(f'вхід: y {Y_INLET_IN:.1f}–{Y_INLET_OUT:.1f}, площа {A_MOUTH/100:.1f} см², висота {2*MOUTH[:,1].max():.0f} мм; виїмка від x={notch_box()[0]:.0f}')
    for x in np.r_[X_INLET,X_INLET+L_INLIP,np.linspace(X_INLET+40,X_DUCT2,10)]:
        s=duct_section(x); A=area(s[:,1:])
        gaps=[z_top(x,y)-abs(z) for y,z in s[:,1:] if z_top(x,y)>0]
        print(f'x={x:6.1f} A={A/100:5.1f} см² y {s[:,1].min():5.1f}..{s[:,1].max():5.1f} мін.стінка {min(gaps) if gaps else float("nan"):6.1f}')
    print('відсік:',[(round(a),round(b,1)) for a,b in bay_profile()])
