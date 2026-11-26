import sys, math, numpy as np
from build123d import *
sys.path.insert(0,'.'); from params import *; from foils import coords, half
def section_wire(y, foil, c, off, n=N_PTS):
    xs,zs=coords(foil,n)
    pts=[Vector(off+c*x, y, c*z) for x,z in zip(xs,zs)]
    sp=Spline(*pts); te=Line(pts[-1],pts[0])
    return Wire([sp,te])
def oml():
    secs=[(-y,f,c,o) for y,f,c,o in SECTIONS[::-1] if y>0]+SECTIONS
    wires=[section_wire(*s) for s in secs]
    return Solid.make_loft(wires, ruled=True)
def engine():
    pts=[(X_NOZZLE+(xl-9.95)*10, r) for xl,r in ENGINE]
    prof=[Vector(pts[0][0],0,0)]+[Vector(x,0,r) for x,r in pts]+[Vector(pts[-1][0],0,0)]
    face=Face(Wire.make_polygon(prof,close=True))
    return Solid.revolve(face,360,Axis.X)
def fin(side=1):
    xs,zs=coords(FIN_FOIL,41)
    def w(z,le,c):
        p=[Vector(le+c*x, FIN_Y*side + c*t*side*0+ c*t, z) for x,t in zip(xs,zs)]
        return Wire([Spline(*p),Line(p[-1],p[0])])
    s=Solid.make_loft([w(0,FIN_ROOT_LE,FIN_ROOT_C), w(FIN_SPAN,FIN_ROOT_LE+FIN_TIP_DX,FIN_TIP_C)],ruled=True)
    return s.rotate(Axis((0,FIN_Y*side,0),(1,0,0)), -FIN_CANT*side)
def build():
    W=oml()
    slot=Box(400,2*Y_SLOT,400,align=(Align.MIN,Align.CENTER,Align.CENTER)).moved(Location((X_SLOT,0,0)))
    W=W-slot
    parts={}
    for nm,y0,y1 in (('elevon_in_R',Y_ELEV_IN,Y_ELEV_SPLIT),('elevon_out_R',Y_ELEV_SPLIT,SPAN_HALF+5),
                     ('elevon_in_L',-Y_ELEV_SPLIT,-Y_ELEV_IN),('elevon_out_L',-SPAN_HALF-5,-Y_ELEV_SPLIT)):
        b=Box(200,y1-y0,400,align=(Align.MIN,Align.MIN,Align.CENTER)).moved(Location((X_HINGE,y0,0)))
        parts[nm]=W&b; W=W-b
    parts['wing']=W
    parts['fin_R']=fin(1); parts['fin_L']=fin(-1)
    parts['engine_jet200']=engine()
    return parts
if __name__=='__main__':
    P=build()
    for k,v in P.items(): print(f'{k:16s} V={v.volume/1e6:8.3f} л  bbox={v.bounding_box().size}')
    asm=Compound(label='ROMB_D5_jet200',children=[v.solid() if hasattr(v,'solid') and False else v for v in P.values()])
    for (k,v),ch in zip(P.items(),asm.children): ch.label=k
    export_step(asm,'romb_d5_jet200.step')
    for k,v in P.items(): export_stl(v,f'stl_{k}.stl',tolerance=0.2,angular_tolerance=0.2)
