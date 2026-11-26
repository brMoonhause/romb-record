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
# ---------- верхня поверхня крила z(x,y) (лінійчата між станціями, як і loft) ----------
_ST=[(y,section_pts(y,f,c,o)) for y,f,c,o in SECTIONS]
def z_top(x,y):
    y=abs(y)
    for (ya,A),(yb,B) in zip(_ST[:-1],_ST[1:]):
        if ya<=y<=yb:
            t=(y-ya)/(yb-ya); n=(len(A)+1)//2
            P=A[:n]*(1-t)+B[:n]*t          # верхня половина (від TE до LE)
            xs,zs=P[::-1,0],P[::-1,2]      # LE -> TE
            if x<xs[0] or x>xs[-1]: return 0.0
            return float(np.interp(x,xs,zs))
    return 0.0
# ---------- обтічник ----------
def _ring(pts_yz,n=72):
    """перевибірка замкненого контуру за кутом навколо (0,0)"""
    a=np.arctan2(pts_yz[:,1],pts_yz[:,0]); o=np.argsort(a); a,P=a[o],pts_yz[o]
    A=np.linspace(-np.pi,np.pi,n,endpoint=False)
    r=np.hypot(P[:,0],P[:,1]); rr=np.interp(A,a,r,period=2*np.pi)
    return np.c_[rr*np.cos(A),rr*np.sin(A)]
def body_section_yz(x,ymax=Y_SLOT+0.5,n=81,rc=6.0):
    """переріз тіла площиною x=const для |y|<=ymax, кути біля стінок ніші заокруглені rc"""
    ys=ymax*np.sin(np.linspace(-np.pi/2,np.pi/2,n))          # густіше біля боків
    zt=np.array([z_top(x,y) for y in ys])
    # заокруглення кутів: біля |y|=ymax верхня лінія плавно опускається на rc
    d=ymax-np.abs(ys); zt=np.where(d<rc, zt-rc+np.sqrt(np.clip(rc*rc-(rc-d)**2,0,None)), zt)
    top=np.c_[ys,zt][::-1]                                   # від +y до -y по верху
    bot=np.c_[ys,-zt]                                        # від -y до +y по низу
    P=np.r_[top,bot[1:-1]]
    return P
def fairing_sections(n=160):
    A=_resample(body_section_yz(X_FAIR0),n)
    i0=np.argmin(np.abs(np.arctan2(A[:,1],A[:,0]))); A=np.roll(A,-i0,axis=0)   # старт на +y, z=0
    ang=np.unwrap(np.arctan2(A[:,1],A[:,0]))
    C=np.c_[R_FAIR1*np.cos(ang),R_FAIR1*np.sin(ang)]
    out=[]
    for i in range(N_FAIR+2):
        t=i/(N_FAIR+1); x=X_FAIR0+t*(X_FAIR1-X_FAIR0); w=smooth(t)
        P=(1-w)*A+w*C; out.append(np.c_[np.full(len(P),x),P[:,0],P[:,1]])
    return out
# ---------- канали ----------
def _rrect(yc,w,h,r,n=64):
    """закруглений прямокутник, точки рівномірно по периметру, старт з (yc+w/2, 0) проти год. стрілки"""
    from math import pi
    hw,hh=w/2-r,h/2-r
    segs=[('L',(hw+r, -hh),(hw+r,hh)),('A',(hw,hh),0),('L',(hw,hh+r),(-hw,hh+r)),('A',(-hw,hh),pi/2),
          ('L',(-hw-r,hh),(-hw-r,-hh)),('A',(-hw,-hh),pi),('L',(-hw,-hh-r),(hw,-hh-r)),('A',(hw,-hh),1.5*pi)]
    pts=[]
    for s in segs:
        if s[0]=='L':
            p=np.linspace(s[1],s[2],40,endpoint=False)
        else:
            a=np.linspace(s[2],s[2]+pi/2,20,endpoint=False); p=np.c_[s[1][0]+r*np.cos(a),s[1][1]+r*np.sin(a)]
        pts.append(p)
    P=np.vstack(pts); P=P[np.r_[len(P)//2*0:len(P)]]
    # старт з середини правої сторони
    P=np.roll(P,-20,axis=0)
    return _resample(P+[yc,0],n)
def _dshape(r,split,n=64):
    """половина кола, плоска сторона на y=split, опукла назовні (+y)"""
    a=np.linspace(-np.pi/2,np.pi/2,60); arc=np.c_[split+r*np.cos(a),r*np.sin(a)]
    flat=np.linspace([split,r],[split,-r],30)[1:-1]
    P=np.r_[arc,flat]; i=np.argmin(np.abs(P[:,1])+ (P[:,0]<split+r/2)*1e3)
    return _resample(np.roll(P,-i,axis=0),n)
def _resample(P,n):
    Q=np.r_[P,P[:1]]; d=np.r_[0,np.cumsum(np.hypot(*np.diff(Q,axis=0).T))]
    s=np.linspace(0,d[-1],n,endpoint=False)
    return np.c_[np.interp(s,d,Q[:,0]),np.interp(s,d,Q[:,1])]
def duct_sections(side=1):
    R=_rrect(0,INLET_Y_OUT-INLET_Y_IN,INLET_H,INLET_R); D=_dshape(DUCT_R_END,DUCT_SPLIT)
    yc0=(INLET_Y_IN+INLET_Y_OUT)/2
    xs=np.linspace(X_DUCT0,X_DUCT2,N_DUCT)      # рівномірно — інакше loft «роздуває» канал
    out=[]
    for x in xs:
        t=smooth((x-X_DUCT1)/(X_DUCT2-X_DUCT1))
        P=(1-t)*(R+[yc0,0])+t*D
        out.append(np.c_[np.full(len(P),x),side*P[:,0],P[:,1]])
    return out
def duct_area(sec): y,z=sec[:,1],sec[:,2]; return 0.5*abs(np.dot(y,np.roll(z,1))-np.dot(z,np.roll(y,1)))
def engine_profile():
    return [(X_NOZZLE+(xl-9.95)*10,r) for xl,r in ENGINE]
if __name__=='__main__':
    print('обтічник: висота на X_FAIR0 =',2*z_top(X_FAIR0,0))
    for s in duct_sections():
        x=s[0,0]; A=duct_area(s)
        # мінімальний запас стінки до поверхні (лише там, де канал уже всередині корпусу)
        gaps=[z_top(x,y)-abs(z) for y,z in s[:,1:]]
        print(f'канал x={x:6.1f}  площа {A/100:5.1f} см²  y {s[:,1].min():5.1f}..{s[:,1].max():5.1f}  мін. стінка {min(gaps):6.1f} мм')
