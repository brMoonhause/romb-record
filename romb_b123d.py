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
    trunk=Solid.make_loft([wire_closed(s) for s in trunk_sections()])
    bay=revolve_x(bay_profile())
    x0,x1,y0,y1=notch_box()
    W=W-box(x0,x1,y0,y1)-box(x0,x1,-y1,-y0)
    for f in flares: W=W-f
    for d in ducts: W=W-d
    W=W-trunk
    W=W-bay
    B=W
    parts={}
    g=HINGE_GAP
    for sd,sn in ((1,'R'),(-1,'L')):
        for nm,y0,y1 in ELEVONS:
            a,b=sorted((sd*y0,sd*y1))
            el=(B&box(X_HINGE,X_HINGE+200,a,b))+Solid.make_loft([wire_closed(q) for q in nose_sections(a,b)])
            for yh,kind in HINGES:
                if kind=='finger' and y0<yh<y1:
                    h=finger_h(yh); el=el-box(X_HINGE-r_nose(yh)-1,X_HINGE+3.5,sd*yh-1.6,sd*yh+1.6,-h-0.6,h+0.6)
            for yh,rs,lh in HORNS:
                if y0<yh<y1: el=el-box(X_HINGE-r_nose(yh)-1,X_HINGE+1,sd*yh-2,sd*yh+2,-r_nose(yh)-1,0)
            parts[f'elevon_{nm}_{sn}']=el
            B=B-box(X_HINGE,X_HINGE+200,a-g,b+g)
            B=B-Solid.make_loft([wire_closed(q) for q in nose_sections(a-g-0.3,b+g+0.3,grow=g)])
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
        (yra,xra),(yrb,xrb)=REAR_SPAR; ar=np.degrees(np.arctan2(xrb-xra,yrb-yra)); Lr=np.hypot(xrb-xra,yrb-yra)
        if sd>0: rs_=Box(T_PLATE,Lr,400,align=(Align.CENTER,Align.MIN,Align.CENTER)).rotate(Axis.Z,-ar).moved(Location((xra,yra,0)))
        else: rs_=Box(T_PLATE,Lr,400,align=(Align.CENTER,Align.MAX,Align.CENTER)).rotate(Axis.Z,ar).moved(Location((xra,-yra,0)))
        S['spar_rear_'+nm]=B&rs_
        for yr in RIBS:
            xf=xa+(yr-ya)*(xb-xa)/(yb-ya) if yr>=ya else le_x(yr)+35
            S[f'rib_{int(yr)}_{nm}']=B&box(xf-2,x_rear(yr),sd*yr-t,sd*yr+t)
        for k,(ysv,xsv) in enumerate(SERVOS):
            lx,ly,lz=SERVO_BOX; S[f'servo_{k}_{nm}']=box(xsv-lx/2,xsv+lx/2,sd*ysv-ly/2,sd*ysv+ly/2,-lz/2,lz/2)
        for yh,kind in HINGES:
            y=sd*yh; rp=pin_r(yh)
            if kind=='end':
                S[f'pin_{int(yh)}_{nm}']=Cylinder(rp,10,rotation=(90,0,0)).moved(Location((X_HINGE,y,0)))
            else:
                h=finger_h(yh)
                S[f'finger_{int(yh)}_{nm}']=box(x_rear(yh)+1,X_HINGE+2.5,y-1,y+1,-h,h)-Cylinder(rp+0.05,4,rotation=(90,0,0)).moved(Location((X_HINGE,y,0)))
                S[f'pin_{int(yh)}_{nm}']=Cylinder(rp,8,rotation=(90,0,0)).moved(Location((X_HINGE,y,0)))
        for (yh,rs,lh),(ysv,xsv) in zip(HORNS,SERVOS):
            y=sd*yh; p0=Vector(xsv+10,y,-rs); p1=Vector(X_HINGE-2,y,-lh); d=p1-p0
            rr=1.5 if lh>6 else 1.25
            S[f'pushrod_{int(yh)}_{nm}']=Cylinder(rr,d.length,align=(Align.CENTER,Align.CENTER,Align.MIN)).moved(Location(p0,(0,np.degrees(np.arctan2(d.X,d.Z)),0)))
            S[f'hornpin_{int(yh)}_{nm}']=Cylinder(1.0,6,rotation=(90,0,0)).moved(Location((X_HINGE-2,y,-lh)))
            S[f'servoarm_{int(yh)}_{nm}']=box(xsv+10-1.5,xsv+10+1.5,y+(3 if sd>0 else -5),y+(5 if sd>0 else -3),-rs-1.5,1.5)
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
    for fn,dd in (('romb_d5_jet200.step',P),('romb_d5_structure.step',S)):
        asm=Compound(label=fn[:-5],children=list(dd.values()))
        for k,ch in zip(dd,asm.children): ch.label=k
        export_step(asm,fn)
