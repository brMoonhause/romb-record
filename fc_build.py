# ---------------- побудова у FreeCAD ----------------
import FreeCAD as App, Part, math
V=App.Vector
def _w_open(P):
    pts=[V(*map(float,p)) for p in P]; bs=Part.BSplineCurve(); bs.interpolate(pts)
    return Part.Wire([bs.toShape(),Part.LineSegment(pts[-1],pts[0]).toShape()])
def _w_closed(P):
    pts=[V(*map(float,p)) for p in P]; bs=Part.BSplineCurve(); bs.interpolate(pts,PeriodicFlag=True)
    return Part.Wire([bs.toShape()])
def _rev_x(prof):
    poly=[V(prof[0][0],0,0)]+[V(x,0,r) for x,r in prof]+[V(prof[-1][0],0,0),V(prof[0][0],0,0)]
    return Part.Face(Part.makePolygon(poly)).revolve(V(0,0,0),V(1,0,0),360)
def _box(x0,x1,y0,y1,z0=-200,z1=200): return Part.makeBox(x1-x0,y1-y0,z1-z0,V(x0,y0,z0))
def geometry(structure=True):
    """лише геометрія (без документа) — можна рахувати у фоновому потоці"""
    W=Part.makeLoft([_w_open(section_pts(s)) for s in full_sections()],True,False)
    ducts=[Part.makeLoft([_w_closed(s) for s in duct_sections(sd)],True,False) for sd in (1,-1)]
    flares=[Part.makeLoft([_w_closed(s) for s in flare_sections(sd)],True,False) for sd in (1,-1)]
    trunk=Part.makeLoft([_w_closed(s) for s in trunk_sections()],True,False)
    bay=_rev_x(bay_profile())
    x0,x1,y0,y1=notch_box()
    W=W.cut(_box(x0,x1,y0,y1)).cut(_box(x0,x1,-y1,-y0))
    for c in flares+ducts+[trunk]: W=W.cut(c)
    W=W.cut(bay)
    B=W
    g=HINGE_GAP; parts=[]
    loft=lambda S:Part.makeLoft([_w_closed(q) for q in S],True,False)
    for sd,sn in ((1,'R'),(-1,'L')):
        for nm,y0,y1 in ELEVONS:
            a,b=sorted((sd*y0,sd*y1))
            el=B.common(_box(X_HINGE,X_HINGE+200,a,b)).fuse(loft(nose_sections(a,b)))
            for yh,kind in HINGES:
                if kind=='finger' and y0<yh<y1:
                    h=finger_h(yh); el=el.cut(_box(X_HINGE-r_nose(yh)-1,X_HINGE+3.5,sd*yh-1.6,sd*yh+1.6,-h-0.6,h+0.6))
            for yh,rs,lh in HORNS:
                if y0<yh<y1: el=el.cut(_box(X_HINGE-r_nose(yh)-1,X_HINGE+1,sd*yh-2,sd*yh+2,-r_nose(yh)-1,0))
            parts.append((f'Elevon_{nm}_{sn}',el.removeSplitter()))
            B=B.cut(_box(X_HINGE,X_HINGE+200,a-g,b+g)).cut(loft(nose_sections(a-g-0.3,b+g+0.3,grow=g)))
    parts=[('Body',B)]+parts
    for side,nm in ((1,'Fin_R'),(-1,'Fin_L')):
        xs,zs=coords(FIN_FOIL,41)
        w=lambda z,le,c:_w_open([(le+c*x,FIN_Y*side+c*t,z) for x,t in zip(xs,zs)])
        s=Part.makeLoft([w(0,FIN_ROOT_LE,FIN_ROOT_C),w(FIN_SPAN,FIN_ROOT_LE+FIN_TIP_DX,FIN_TIP_C)],True,True)
        s.rotate(V(0,FIN_Y*side,0),V(1,0,0),-FIN_CANT*side); parts.append((nm,s))
    parts.append(('Engine_jet200',_rev_x(engine_profile())))
    S=[]
    if structure:
        t=T_PLATE/2
        S.append(('Keel',B.common(_box(KEEL[0],KEEL[1],-t,t))))
        for sd,n in ((1,'R'),(-1,'L')):
            y=sd*SIDE_WEBS[0]; S.append(('Web_'+n,B.common(_box(SIDE_WEBS[1],SIDE_WEBS[2],y-t,y+t))))
        for x,hw in FRAMES: S.append((f'Frame_{int(x)}',B.common(_box(x-t,x+t,-hw,hw))))
        (ya,xa),(yb,xb)=FRONT_SPAR; ang=math.degrees(math.atan2(xb-xa,yb-ya)); L=math.hypot(xb-xa,yb-ya)
        for sd,n in ((1,'R'),(-1,'L')):
            sp=Part.makeBox(T_PLATE,L,400,V(-t,0,-200)); sp.rotate(V(0,0,0),V(0,0,1),-sd*ang if sd>0 else 0)
            if sd<0: sp=Part.makeBox(T_PLATE,L,400,V(-t,-L,-200)); sp.rotate(V(0,0,0),V(0,0,1),ang)
            sp.translate(V(xa,sd*ya,0)); S.append(('Spar_front_'+n,B.common(sp)))
            (yra,xra),(yrb,xrb)=REAR_SPAR; ar=math.degrees(math.atan2(xrb-xra,yrb-yra)); Lr=math.hypot(xrb-xra,yrb-yra)
            if sd>0: rsp=Part.makeBox(T_PLATE,Lr,400,V(-t,0,-200)); rsp.rotate(V(0,0,0),V(0,0,1),-ar)
            else: rsp=Part.makeBox(T_PLATE,Lr,400,V(-t,-Lr,-200)); rsp.rotate(V(0,0,0),V(0,0,1),ar)
            rsp.translate(V(xra,sd*yra,0)); S.append(('Spar_rear_'+n,B.common(rsp)))
            for yr in RIBS:
                xf=xa+(yr-ya)*(xb-xa)/(yb-ya) if yr>=ya else le_x(yr)+35
                S.append((f'Rib_{int(yr)}_{n}',B.common(_box(xf-2,x_rear(yr),sd*yr-t,sd*yr+t))))
            for k,(ysv,xsv) in enumerate(SERVOS):
                lx,ly,lz=SERVO_BOX; S.append((f'Servo_{k}_{n}',_box(xsv-lx/2,xsv+lx/2,sd*ysv-ly/2,sd*ysv+ly/2,-lz/2,lz/2)))
            for yh,kind in HINGES:
                y=sd*yh; rp=pin_r(yh)
                if kind=='end':
                    S.append((f'Pin_{int(yh)}_{n}',Part.makeCylinder(rp,10,V(X_HINGE,y-5,0),V(0,1,0))))
                else:
                    h=finger_h(yh)
                    S.append((f'Finger_{int(yh)}_{n}',_box(x_rear(yh)+1,X_HINGE+2.5,y-1,y+1,-h,h).cut(Part.makeCylinder(rp+0.05,4,V(X_HINGE,y-2,0),V(0,1,0)))))
                    S.append((f'Pin_{int(yh)}_{n}',Part.makeCylinder(rp,8,V(X_HINGE,y-4,0),V(0,1,0))))
            for (yh,rs,lh),(ysv,xsv) in zip(HORNS,SERVOS):
                y=sd*yh; p0=V(xsv+10,y,-rs); p1=V(X_HINGE-2,y,-lh); d=p1-p0
                S.append((f'Pushrod_{int(yh)}_{n}',Part.makeCylinder(1.5 if lh>6 else 1.25,d.Length,p0,d)))
                S.append((f'Hornpin_{int(yh)}_{n}',Part.makeCylinder(1.0,6,V(X_HINGE-2,y-3,-lh),V(0,1,0))))
    return parts,S
def make_doc(parts,S,docname='ROMB_D5_jet200'):
    if docname in App.listDocuments(): App.closeDocument(docname)
    doc=App.newDocument(docname)
    col={'Body':(0.72,0.77,0.82),'Elev':(0.89,0.63,0.24),'Fin_':(0.35,0.56,0.82),'Engi':(0.75,0.22,0.17),
         'Serv':(0.48,0.36,0.84),'Pin_':(0.6,0.6,0.62),'Fing':(0.6,0.6,0.62),'Push':(0.1,0.1,0.12),'Horn':(0.85,0.34,0.17)}
    grp=doc.addObject('App::DocumentObjectGroup','Structure') if S else None; sn={n for n,_ in S}
    for nm,sh in parts+S:
        o=doc.addObject('Part::Feature',nm); o.Shape=sh
        if grp and nm in sn: grp.addObject(o)
        if App.GuiUp: o.ViewObject.ShapeColor=col.get(nm[:4],(0.15,0.68,0.38))
    doc.recompute()
    if App.GuiUp:
        import FreeCADGui as Gui; Gui.ActiveDocument.ActiveView.viewIsometric(); Gui.SendMsgToActiveView('ViewFit')
    for nm,sh in parts: print(f'{nm:14s} V={sh.Volume/1e6:7.3f} L valid={sh.isValid()}')
    return doc
def build(docname='ROMB_D5_jet200', structure=True):
    parts,S=geometry(structure); return make_doc(parts,S,docname)
if __name__=='__main__': build()
