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
def box(x0,x1,y0,y1,z0=-200,z1=200):
    return Box(x1-x0,y1-y0,z1-z0,align=(Align.MIN,Align.MIN,Align.MIN)).moved(Location((x0,y0,z0)))
def build():
    W=Solid.make_loft([wire_open_te(section_pts(s)) for s in full_sections()],ruled=False)
    ducts=[Solid.make_loft([wire_closed(s) for s in duct_sections(sd)]) for sd in (1,-1)]
    flares=[Solid.make_loft([wire_closed(s) for s in flare_sections(sd)]) for sd in (1,-1)]
    bay=revolve_x(bay_profile())
    x0,x1,y0,y1=notch_box()
    W=W-box(x0,x1,y0,y1)-box(x0,x1,-y1,-y0)
    for f in flares: W=W-f
    for d in ducts: W=W-d
    W=W-bay
    B=W
    parts={}
    yi,ys,ss=Y_ELEV_IN,Y_ELEV_SPLIT,FIN_STRIP
    B=B-box(X_HINGE,X_HINGE+HINGE_GAP,yi,SPAN_HALF+5)-box(X_HINGE,X_HINGE+HINGE_GAP,-SPAN_HALF-5,-yi)
    for nm,y0,y1 in (('elevon_in_R',yi,ys-ss),('elevon_out_R',ys+ss,SPAN_HALF+5),
                     ('elevon_in_L',-ys+ss,-yi),('elevon_out_L',-SPAN_HALF-5,-ys-ss)):
        b=box(X_HINGE+HINGE_GAP,X_HINGE+200,y0,y1); parts[nm]=B&b; B=B-b
    parts['body']=B
    # --- внутрішній набір ---
    t=T_PLATE/2; S={}
    S['keel']=B&box(KEEL[0],KEEL[1],-t,t)
    for sd,nm in ((1,'R'),(-1,'L')):
        y=sd*SIDE_WEBS[0]; S['web_'+nm]=B&box(SIDE_WEBS[1],SIDE_WEBS[2],y-t,y+t)
    for x,hw in FRAMES: S[f'frame_{int(x)}']=B&box(x-t,x+t,-hw,hw)
    (ya,xa),(yb,xb)=FRONT_SPAR; ang=np.degrees(np.arctan2(xb-xa,yb-ya)); L=np.hypot(xb-xa,yb-ya)
    for sd,nm in ((1,'R'),(-1,'L')):
        sp=Box(T_PLATE,L,400,align=(Align.CENTER,Align.MIN,Align.CENTER)).rotate(Axis.Z,-sd*ang).moved(Location((xa,sd*ya,0)))
        if sd<0: sp=Box(T_PLATE,L,400,align=(Align.CENTER,Align.MAX,Align.CENTER)).rotate(Axis.Z,ang).moved(Location((xa,-ya,0)))
        S['spar_front_'+nm]=B&sp
        S['spar_rear_'+nm]=B&box(REAR_SPAR[0]-t,REAR_SPAR[0]+t,*(sorted((sd*REAR_SPAR[1],sd*REAR_SPAR[2]))))
        for yr in RIBS:
            xf=xa+(yr-ya)*(xb-xa)/(yb-ya) if yr>=ya else le_x(yr)+35
            S[f'rib_{int(yr)}_{nm}']=B&box(xf-2,REAR_SPAR[0],sd*yr-t,sd*yr+t)
        for k,(ysv,xsv) in enumerate(SERVOS):
            lx,ly,lz=SERVO_BOX; S[f'servo_{k}_{nm}']=box(xsv-lx/2,xsv+lx/2,sd*ysv-ly/2,sd*ysv+ly/2,-lz/2,lz/2)
        for y0,y1 in ((yi,ys-ss),(ys+ss,SPAN_HALF)):
            a,b2=sorted((sd*y0,sd*y1))
            S[f'hinge_{int(y0)}_{nm}']=Cylinder(1.0,b2-a,rotation=(90,0,0)).moved(Location((X_HINGE+0.5,(a+b2)/2,0)))
    for side,nm in ((1,'fin_R'),(-1,'fin_L')):
        xs,zs=coords(FIN_FOIL,41)
        w=lambda z,le,c:wire_open_te([(le+c*x,FIN_Y*side+c*tt,z) for x,tt in zip(xs,zs)])
        s=Solid.make_loft([w(0,FIN_ROOT_LE,FIN_ROOT_C),w(FIN_SPAN,FIN_ROOT_LE+FIN_TIP_DX,FIN_TIP_C)],ruled=True)
        parts[nm]=s.rotate(Axis((0,FIN_Y*side,0),(1,0,0)),-FIN_CANT*side)
    parts['engine_jet200']=revolve_x(engine_profile())
    return parts,S
if __name__=='__main__':
    P,S=build()
    allp={**P,**S}
    for k,v in allp.items(): print(f'{k:18s} valid={v.is_valid}',flush=True)
    import os; os.makedirs('stl',exist_ok=True)
    for k,v in allp.items(): export_stl(v,f'stl/{k}.stl',tolerance=0.15,angular_tolerance=0.15)
    print('body V',vol(P['body'])/1e6, flush=True)
    asm=Compound(label='ROMB_D5_jet200',children=list(allp.values()))
    for k,ch in zip(allp,asm.children): ch.label=k
    export_step(asm,'romb_d5_jet200.step')
