"""Розрахунок повітрозабірників: витрата, горло, губи NACA 1-series, перебір форм каналу."""
import numpy as np, sys
sys.path.insert(0,'.')
from params import *
from geom import z_top, smooth
g=1.4; R=287.0; p0=101325; T0=288.15
MDOT=0.43                       # кг/с, K-210 (паспорт K-210G4+)
def flux(M):                    # кг/с/м² для ізоентропічного потоку
    return p0/np.sqrt(R*T0)*np.sqrt(g)*M*(1+0.2*M*M)**-3
def mach_from_flux(f):
    Ms=np.linspace(0.01,0.99,2000); return np.interp(f,flux(Ms),Ms)
A_face=np.pi*(0.0421**2-0.0125**2)
print(f'вхід двигуна (кільце губи мінус конус): {A_face*1e4:.1f} см²,  M на вході статично {mach_from_flux(MDOT/A_face):.3f}')
for M in (0.35,0.40,0.45,0.50):
    print(f'  горло (обидва) для M_горла={M}: {MDOT/flux(M)*1e4:.1f} см²')
A_T=2*13.0e-4                   # вибір: 13 см² горло на канал
CR=1.10                         # контракція губи A_highlight/A_throat
print(f'вибрано: горло {A_T*1e4/2:.1f} см² × 2  → M_горла статично {mach_from_flux(MDOT/A_T):.3f}; highlight {CR*A_T*1e4/2:.1f} см² × 2')
V=785/3.6; A0=MDOT/(1.225*V)
print(f'на 785 км/год: потрібна площа захоплення {A0*1e4:.1f} см², MFR = A0/A_hl = {A0/(CR*A_T):.2f}')
# ---- перебір: співвідношення сторін горла і положення входу ----
def body_le_y(x):               # край планформи тіла (передня кромка) між станціями 30..60..80
    pts=[(0,0),(53,30),(282.3,60),(417,80),(600,127.5)]
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]; return float(np.interp(x,xs,ys))
def rrect_area(w,h,r): return w*h-(4-np.pi)*r*r
def dims(A,ar):                 # w,h з площі і h/w, r=0.25*w
    w=np.sqrt(A/(ar-(4-np.pi)*0.0625)); return w,ar*w,0.25*w
X_END=X_DUCT2; A_END=np.pi*38.8**2/2   # вихід каналу: половина кола
rows=[]
for xm in (250,280,310,340):
  for ar in (1.0,1.5,2.0,2.5,3.0):
    At=A_T/2*1e6; w,h,r=dims(At*CR,ar)     # highlight
    lip=0.088*np.sqrt(4*At/np.pi)          # товщина губи NACA 1-85: (D-d)/2 = d*(1/0.85-1)/2
    yin=body_le_y(xm)+lip+1.0; yc0=yin+w/2
    hump=0; minwall=99; L=X_END-xm
    xs=np.linspace(xm,X_END,60); fr=0
    for x in xs:
        t=smooth((x-xm-40)/(X_END-xm-40))
        A=At+(A_END-At)*max(0,(x-xm-25)/(L-25))      # площа росте лінійно після горла
        yc=(1-t)*yc0+t*(DUCT_SPLIT+4*38.8/(3*np.pi))
        ww,hh,rr=dims(A,ar*(1-t)+t*1.0)
        for yy in np.linspace(yc-ww/2,yc+ww/2,7):
            zt=z_top(x,yy)
            if zt>0 and x>xm+60:
                hump=max(hump,hh/2+DUCT_WALL-zt); minwall=min(minwall,zt-hh/2)
        Dh=4*A/(2*(ww+hh)); fr+=0.02*(xs[1]-xs[0])/Dh*(At/A)**2
    dy=yc0-(DUCT_SPLIT+4*38.8/(3*np.pi))
    k_bend=0.25*(dy/L)**2*100          # груба оцінка втрат S-вигину (∝ квадрату відносного зміщення)
    k_diff=0.18*(1-At/A_END)**2
    K=fr+k_bend+k_diff+0.03            # +губа
    rows.append((xm,ar,w,h,yin,yin+w,max(hump,0),dy,K))
print('\n x_вх  h/w   w×h горла,мм   y губ    горб над крилом  зміщ.  K втрат (частки q_горла)')
for xm,ar,w,h,yi,yo,hp,dy,K in rows:
    print(f'{xm:5.0f} {ar:4.1f} {w/np.sqrt(1.1):6.1f}×{h/np.sqrt(1.1):5.1f}  {yi:5.1f}–{yo:5.1f}   {hp:6.1f} мм       {dy:5.1f}  {K:.3f}')
