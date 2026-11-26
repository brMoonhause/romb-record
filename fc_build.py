# ---------------- побудова у FreeCAD ----------------
import FreeCAD as App, Part
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
def build(docname='ROMB_D5_jet200'):
    if docname in App.listDocuments(): App.closeDocument(docname)
    doc=App.newDocument(docname)
    W=Part.makeLoft([_w_open(section_pts(*s)) for s in full_sections()],True,True)
    W=W.cut(Part.makeBox(600,2*Y_SLOT,400,V(X_FAIR0,-Y_SLOT,-200)))
    F=Part.makeLoft([_w_closed(s) for s in fairing_sections()],True,False)
    cowl=Part.makeCylinder(R_FAIR1,X_COWL_END-X_FAIR1,V(X_FAIR1,0,0),V(1,0,0))
    ducts=[Part.makeLoft([_w_closed(s) for s in duct_sections(side)],True,False) for side in (1,-1)]
    PL=_rev_x(PLENUM)
    # канали й камеру вирізаємо з кожної частини окремо, потім об'єднуємо (так стабільніше в OCC FreeCAD)
    for d in ducts: W=W.cut(d); F=F.cut(d)
    F=F.cut(PL); cowl=cowl.cut(PL)
    B=W.fuse(F).fuse(cowl)
    parts=[]
    for nm,y0,y1 in (('Elevon_in_R',Y_ELEV_IN,Y_ELEV_SPLIT),('Elevon_out_R',Y_ELEV_SPLIT,SPAN_HALF+5),
                     ('Elevon_in_L',-Y_ELEV_SPLIT,-Y_ELEV_IN),('Elevon_out_L',-SPAN_HALF-5,-Y_ELEV_SPLIT)):
        b=Part.makeBox(200,y1-y0,400,V(X_HINGE,y0,-200)); parts.append((nm,B.common(b))); B=B.cut(b)
    parts=[('Body',B)]+parts
    for side,nm in ((1,'Fin_R'),(-1,'Fin_L')):
        xs,zs=coords(FIN_FOIL,41)
        w=lambda z,le,c:_w_open([(le+c*x,FIN_Y*side+c*t,z) for x,t in zip(xs,zs)])
        s=Part.makeLoft([w(0,FIN_ROOT_LE,FIN_ROOT_C),w(FIN_SPAN,FIN_ROOT_LE+FIN_TIP_DX,FIN_TIP_C)],True,True)
        s.rotate(V(0,FIN_Y*side,0),V(1,0,0),-FIN_CANT*side); parts.append((nm,s))
    parts.append(('Engine_jet200',_rev_x(engine_profile())))
    col={'Body':(0.72,0.77,0.82),'Elev':(0.89,0.63,0.24),'Fin_':(0.35,0.56,0.82),'Engi':(0.75,0.22,0.17)}
    for nm,sh in parts:
        o=doc.addObject('Part::Feature',nm); o.Shape=sh
        if App.GuiUp: o.ViewObject.ShapeColor=col[nm[:4]]
    doc.recompute()
    if App.GuiUp:
        import FreeCADGui as Gui; Gui.ActiveDocument.ActiveView.viewIsometric(); Gui.SendMsgToActiveView('ViewFit')
    for nm,sh in parts: print(f'{nm:14s} V={sh.Volume/1e6:7.3f} L valid={sh.isValid()}')
    return doc
build()
