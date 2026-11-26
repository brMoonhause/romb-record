"""РОМБ D5 / jet 200 — збирання моделі в build123d з таблиці params.py → STEP."""
import sys, numpy as np
from build123d import *
sys.path.insert(0,'.')
from params import *
from geom import *
def vol(s):
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp
    p=GProp_GProps(); BRepGProp.VolumeProperties_s(s.wrapped,p,1e-6,False); return p.Mass()
def wire_open_te(P):
    pts=[Vector(*p) for p in P]; return Wire([Spline(*pts),Line(pts[-1],pts[0])])
def wire_closed(P):
    return Wire([Spline(*[Vector(*p) for p in P],periodic=True)])
def revolve_x(prof):
    poly=[Vector(prof[0][0],0,0)]+[Vector(x,0,r) for x,r in prof]+[Vector(prof[-1][0],0,0)]
    return Solid.revolve(Face(Wire.make_polygon(poly,close=True)),360,Axis.X)
def build():
    W=Solid.make_loft([wire_open_te(section_pts(*s)) for s in full_sections()],ruled=True)
    # центральна частина за X_FAIR0 замінюється обтічником, далі — ніша
    W=W-Box(600,2*Y_SLOT,400,align=(Align.MIN,Align.CENTER,Align.CENTER)).moved(Location((X_FAIR0,0,0)))
    F=Solid.make_loft([wire_closed(s) for s in fairing_sections()])
    cowl=Cylinder(R_FAIR1,X_COWL_END-X_FAIR1,rotation=(0,90,0),align=(Align.CENTER,Align.CENTER,Align.MIN)).moved(Location((X_FAIR1,0,0)))
    B=W+F+cowl
    for side in (1,-1):
        D=Solid.make_loft([wire_closed(s) for s in duct_sections(side)])
        B=B-D
    B=B-revolve_x(PLENUM)
    parts={}
    for nm,y0,y1 in (('elevon_in_R',Y_ELEV_IN,Y_ELEV_SPLIT),('elevon_out_R',Y_ELEV_SPLIT,SPAN_HALF+5),
                     ('elevon_in_L',-Y_ELEV_SPLIT,-Y_ELEV_IN),('elevon_out_L',-SPAN_HALF-5,-Y_ELEV_SPLIT)):
        b=Box(200,y1-y0,400,align=(Align.MIN,Align.MIN,Align.CENTER)).moved(Location((X_HINGE,y0,0)))
        parts[nm]=B&b; B=B-b
    parts['body']=B
    for side,nm in ((1,'fin_R'),(-1,'fin_L')):
        xs,zs=coords(FIN_FOIL,41)
        w=lambda z,le,c:wire_open_te([(le+c*x,FIN_Y*side+c*t,z) for x,t in zip(xs,zs)])
        s=Solid.make_loft([w(0,FIN_ROOT_LE,FIN_ROOT_C),w(FIN_SPAN,FIN_ROOT_LE+FIN_TIP_DX,FIN_TIP_C)],ruled=True)
        parts[nm]=s.rotate(Axis((0,FIN_Y*side,0),(1,0,0)),-FIN_CANT*side)
    parts['engine_jet200']=revolve_x(engine_profile())
    return parts
if __name__=='__main__':
    P=build()
    for k,v in P.items(): print(f'{k:16s} V={vol(v)/1e6:8.3f} л valid={v.is_valid} bbox={tuple(round(a,1) for a in v.bounding_box().size)}')
    asm=Compound(label='ROMB_D5_jet200',children=list(P.values()))
    for k,ch in zip(P,asm.children): ch.label=k
    export_step(asm,'romb_d5_jet200.step')
    for k,v in P.items(): export_stl(v,f'stl_{k}.stl',tolerance=0.2,angular_tolerance=0.2)
