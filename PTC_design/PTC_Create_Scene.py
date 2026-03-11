# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 11:13:09 2023
@author: SAINATH
"""
# -*- coding: utf-8 -*-
try:
    import freecad
except:
    pass
import FreeCAD
import pathlib as Path
import Part
from FreeCAD import Base
import Spreadsheet
import numpy as np
import math
from itertools import product
import numpy as np
import os
cwd = Path.Path().resolve()
path = str(cwd.parents[0])

# define focus
F = 1742.0
# define receiver diameters
D = np.arange(50.0, 120.1, 10.0)

# define receiver locations
H = np.arange(1672.0, 1812.1, 17.5)

# define rim angles for aperture width 
Rim = np.arange(10.0, 150.1, 10.0)

# If aperture width is constant for LS-III
# W = 5760

parameters_set = product(D,H,Rim)
def single_design(arg):
    (
        D,
        H,
        Rim
    ) = arg
    # Alias for active Document
    doc = FreeCAD.newDocument()    
    print(arg)

    # aperture width
    W = 4*F*np.tan(np.deg2rad(Rim)/2)
    
    # define the aperture length
    L = 10000.0

    # define glass tube diameters
    Dgi = D + 38
    Dgo = Dgi + 7
    
    # Alias for active Document
    doc = App.activeDocument()
    # Alias for gui active Document
    #gui = Gui.activeDocument()
    
    # create a Part object that is a Parabola in the XY plane (the parabola is infinite).
    p = Part.Parabola()
    # define de Focal distance in the X axe
    p.Focal = F
    # create a Part object defined by the edge of the parabola with the limits in the Y axe
    edge_p = Part.Edge(p, -W/2, W/2)
    # adds a Part object type to the document and assigns the shape representation of the edge_p
    MyShape = doc.addObject("Part::Feature","Parabola")
    MyShape.Shape = edge_p
    
    # transformating the edge_p
    edge_p.rotate(App.Vector(0.0,0.0,0.0),App.Vector(0.0,1.0,0.0),-90.0)
    edge_p.rotate(App.Vector(0.0,0,0.0),App.Vector(0.0,0.0,1.0),+90.0)
    edge_p.translate(App.Base.Vector(0.0,0.0,0.0))
    MyShape.Shape = edge_p
    
    
    #create an Extrusion object to extrude the edge_p
    Extruded_parabola_1 = doc.addObject("Part::Extrusion","Extruded_parabola_1")
    doc.Extruded_parabola_1.Label = "Extruded_parabola_1(thinfilm_reflective_SiO2_TiO2_SiO2_Al)"
    # define the extrusion for the parabola
    doc.Extruded_parabola_1.Base = doc.Parabola
    doc.Extruded_parabola_1.Dir = (0.0,L,0.0)
    doc.Extruded_parabola_1.Solid = (False)
    doc.Extruded_parabola_1.TaperAngle = (0.0)
    
    # to recalculate the whole document
    FreeCAD.ActiveDocument.recompute()
      
    
    # create a circle for the absorber
    circle_abs = Part.makeCircle(D/2.0, Base.Vector(0.0,0.0,H), Base.Vector(0.0,1.0,0.0))
    #create a FreeCAD object with Part attributes
    circle_abs_Part = doc.addObject("Part::Feature","circle_abs_Part")
    # asign the circle to circle_abs_Part
    circle_abs_Part.Shape = circle_abs
    
    #create a Extrusion object for the extrusion
    Abs_circle_abs_Extrude = doc.addObject("Part::Extrusion","Abs_circle_abs_Extrude")
    doc.Abs_circle_abs_Extrude.Label = "Abs_circle_abs_Extrude(thinfilm_absorber_TiAlN_TiAlON_Si3N4)"
    # define the extrusion of the parabola
    doc.Abs_circle_abs_Extrude.Base = FreeCAD.ActiveDocument.circle_abs_Part
    doc.Abs_circle_abs_Extrude.Dir = (0.0,L,0.0)
    doc.Abs_circle_abs_Extrude.Solid = (False)
    doc.Abs_circle_abs_Extrude.TaperAngle = (0.0)
    
    # to recalculate the whole document
    FreeCAD.ActiveDocument.recompute()
      
    #create a FreeCAD object with Cylinder attributes
    cylinder_out = doc.addObject("Part::Cylinder","cylinder_out")
    cylinder_out.Radius = Dgo/2
    cylinder_out.Height = L
    cylinder_out.Angle = 360.0
    cylinder_out.Placement = App.Placement(App.Vector(0.0,0.0,H),App.Rotation(App.Vector(1,0.0,0.0),-90.0))
    
    #create a FreeCAD object with Cylinder attributes
    cylinder_in = doc.addObject("Part::Cylinder","cylinder_in")
    cylinder_in.Radius = Dgi/2
    cylinder_in.Height = L 
    cylinder_in.Angle = 360.0
    cylinder_in.Placement = App.Placement(App.Vector(0.0,0.0-1.0,H),App.Rotation(App.Vector(1.0,0.0,0.0),-90.0))
    
    #create a FreeCAD object with Part Cut attributes
    tube_glass = doc.addObject("Part::Cut","tube_glass")
    doc.tube_glass.Label = "tube_glass(SiO2_Malitson)"
    doc.tube_glass.Base = doc.cylinder_out
    doc.tube_glass.Tool = doc.cylinder_in
    #gui.hide("cylinder_out")
    #gui.hide("cylinder_in")
    FreeCAD.ActiveDocument.recompute()
    
    # to make transparent the glass tube
    #gui.getObject("tube_glass").Transparency = 80
    
    # this is the face for the outter coating
    o_coating_face = doc.cylinder_out.Shape.Faces[0]
    #create a FreeCAD object for the coating face
    doc.addObject('Part::Feature','O_coating_face').Shape = o_coating_face
    doc.O_coating_face.Label = "O_coating_face(thinfilm_antireflective_MgF2_SiO2_95nm_inf)"
    # to make transparent the glass tube
    #gui.getObject("O_coating_face").Transparency = 40
    FreeCAD.ActiveDocument.recompute()
    
    # this is the face for the inner coating
    i_coating_face = doc.cylinder_in.Shape.Faces[0]
    #create a FreeCAD object for the coating face
    doc.addObject('Part::Feature','I_coating_face').Shape = i_coating_face
    doc.I_coating_face.Label = "I_coating_face(thinfilm_antireflective_MgF2_SiO2_95nm_inf)"
    # to make transparent the glass tube
    #gui.getObject("I_coating_face").Transparency = 40
    FreeCAD.ActiveDocument.recompute()
    
    FreeCAD.ActiveDocument.recompute()
    #Gui.activeDocument().activeView().viewAxonometric()
    #Gui.SendMsgToActiveView("ViewFit")
    
    aperture_th = W * L
    doc.addObject('Spreadsheet::Sheet','Spreadsheet')
    doc.Spreadsheet.set('A1', 'aperture_th')
    doc.Spreadsheet.set('B1', str(aperture_th))
    
    a = list(arg) + [round(aperture_th,2)] 
    b = round(a[2], 1)
    a[2] = b
    print(a)
    if 1 < 2:
        Label_drawing = "PTC_Folder/design_{0}".format(a)+".FCStd"
        doc.saveAs(Label_drawing)

case = 1
total_cases = len(D) * len(H) * len(Rim)
# total_cases = len(list(parameters_set))
for arg in parameters_set:
    case = case + 1
    print(case / total_cases * 100, "%")
    actual_design = single_design(arg)
