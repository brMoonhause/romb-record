"""Перебір: губа за шпангоутом 450 біля переходу до стріловидності, канали в камеру перед двигуном."""
import sys; sys.path.insert(0,'.')
from geom import z_top, body_le_y, rr_dims, smooth
import numpy as np
A_T=1300.0; CR=1.10; LIP=4.0; GAP=1.0; WALL=4.0; LIN=25.0
X_DEND=700.0; Y_DEND=36.0; A_DEND=2000.0
print(' x губи  h/w   горло w×h     y губи      макс.горб  об\'єм горба  макс.кут вигину  K≈')
best=[]
for xm in (465,480,500,520):
  for ar in (1.0,1.5,2.0,2.5):
    w,h,r=rr_dims(A_T*CR,ar); yin=body_le_y(xm)+LIP+GAP; yc0=yin+w/2
    xs=np.linspace(xm,X_DEND,120); hmax=0; hv=0; ang=0; fr=0
    for i,x in enumerate(xs):
        u=np.clip((x-xm-LIN)/(X_DEND-xm-LIN),0,1); s=smooth(u)
        yc=yc0+(Y_DEND-yc0)*s
        A=A_T if x>xm+LIN else A_T*CR; A=A+(A_DEND-A_T)*u
        ww,hh,rr=rr_dims(A,ar*(1-s)+1.0*s)
        hm=0
        for yy in np.linspace(yc-ww/2,yc+ww/2,9):
            zt=z_top(x,yy); hm=max(hm,hh/2+WALL-zt)
        if x>xm+30: hmax=max(hmax,hm); hv+=max(hm,0)*ww*(xs[1]-xs[0])
        if i: ang=max(ang,np.degrees(np.arctan(abs(yc-yprev)/(xs[1]-xs[0]))))
        yprev=yc; Dh=4*A/(2*(ww+hh)); fr+=0.02*(xs[1]-xs[0])/Dh*(A_T/A)**2
    K=fr+0.0075*ang*1.0+0.18*(1-A_T/A_DEND)**2+0.03+ (A_DEND/A_T)**-2*0.9   # тертя + вигин + дифузор + губа + скид у камеру
    best.append((K,xm,ar))
    print(f'{xm:6.0f} {ar:4.1f}  {w/np.sqrt(CR):5.1f}×{h/np.sqrt(CR):5.1f}  {yin:6.1f}–{yin+w:5.1f}   {hmax:6.1f} мм   {hv/1e3:6.0f} см³     {ang:5.1f}°       {K:.2f}')
