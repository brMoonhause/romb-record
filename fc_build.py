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
def build(docname='ROMB_D5_jet200', structure=True):
    if docname in App.listDocuments(): App.closeDocument(docname)
    doc=App.newDocument(docname)
    W=Part.makeLoft([_w_open(section_pts(s)) for s in full_sections()],True,False)
    ducts=[Part.makeLoft([_w_closed(s) for s in duct_sections(sd)],True,False) for sd in (1,-1)]
    flares=[Part.makeLoft([_w_closed(s) for s in flare_sections(sd)],True,False) for sd in (1,-1)]
    bay=_rev_x(bay_profile())
    x0,x1,y0,y1=notch_box()
    W=W.cut(_box(x0,x1,y0,y1)).cut(_box(x0,x1,-y1,-y0))
    for c in flares+ducts: W=W.cut(c)
    W=W.cut(bay)
    B=W
    yi,ys,ss=Y_ELEV_IN,Y_ELEV_SPLIT,FIN_STRIP
    B=B.cut(_box(X_HINGE,X_HINGE+HINGE_GAP,yi,SPAN_HALF+5)).cut(_box(X_HINGE,X_HINGE+HINGE_GAP,-SPAN_HALF-5,-yi))
    parts=[]
    for nm,y0,y1 in (('Elevon_in_R',yi,ys-ss),('Elevon_out_R',ys+ss,SPAN_HALF+5),
                     ('Elevon_in_L',-ys+ss,-yi),('Elevon_out_L',-SPAN_HALF-5,-ys-ss)):
        b=_box(X_HINGE+HINGE_GAP,X_HINGE+200,y0,y1); parts.append((nm,B.common(b))); B=B.cut(b)
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
            S.append(('Spar_rear_'+n,B.common(_box(REAR_SPAR[0]-t,REAR_SPAR[0]+t,*sorted((sd*REAR_SPAR[1],sd*REAR_SPAR[2]))))))
            for yr in RIBS:
                xf=xa+(yr-ya)*(xb-xa)/(yb-ya) if yr>=ya else le_x(yr)+35
                S.append((f'Rib_{int(yr)}_{n}',B.common(_box(xf-2,REAR_SPAR[0],sd*yr-t,sd*yr+t))))
            for k,(ysv,xsv) in enumerate(SERVOS):
                lx,ly,lz=SERVO_BOX; S.append((f'Servo_{k}_{n}',_box(xsv-lx/2,xsv+lx/2,sd*ysv-ly/2,sd*ysv+ly/2,-lz/2,lz/2)))
            for y0,y1 in ((yi,ys-ss),(ys+ss,SPAN_HALF)):
                a,b2=sorted((sd*y0,sd*y1))
                S.append((f'Hinge_{int(y0)}_{n}',Part.makeCylinder(1.0,b2-a,V(X_HINGE+0.5,a,0),V(0,1,0))))
    col={'Body':(0.72,0.77,0.82),'Elev':(0.89,0.63,0.24),'Fin_':(0.35,0.56,0.82),'Engi':(0.75,0.22,0.17),
         'Serv':(0.17,0.24,0.31),'Hing':(0.5,0.55,0.55)}
    grp=doc.addObject('App::DocumentObjectGroup','Structure') if S else None
    for nm,sh in parts+S:
        o=doc.addObject('Part::Feature',nm); o.Shape=sh
        if grp and (nm,sh) in S: grp.addObject(o)
        if App.GuiUp: o.ViewObject.ShapeColor=col.get(nm[:4],(0.15,0.68,0.38))
    doc.recompute()
    if App.GuiUp:
        import FreeCADGui as Gui; Gui.ActiveDocument.ActiveView.viewIsometric(); Gui.SendMsgToActiveView('ViewFit')
    for nm,sh in parts: print(f'{nm:14s} V={sh.Volume/1e6:7.3f} L valid={sh.isValid()}')
    return doc
if __name__=='__main__': build()
