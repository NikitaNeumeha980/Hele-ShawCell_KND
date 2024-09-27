#!/usr/bin/env python

import sys
import math
import salome

salome.salome_init()

###
### SHAPER component
###

from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

####################################################
##            Begin of variables section          ##
####################################################
drainageLength =\
    model.addParameter(
        Part_1_doc,
        "drainageLength",
        "1000",
        "channel length"
    )
channelLength =\
    model.addParameter(
        Part_1_doc,
        "channelLength",
        "4000",
        "channel length"
    )
channelWidth =\
    model.addParameter(
        Part_1_doc,
        "channelWidth",
        "200",
        "channel width"
    )
channelDepth =\
    model.addParameter(
        Part_1_doc,
        "channelDepth",
        "10000",
        "channel depth"
    )
NbStep =\
    model.addParameter(
        Part_1_doc,
        "NbStep",
        "4",
        "number of steps"
    )
step =\
    model.addParameter(
        Part_1_doc,
        "step",
        "(channelLength - channelWidth)/NbStep",
        "size of step"
    )
####################################################
##             End of variables section           ##
####################################################

### Create Box
Box_1 =\
    model.addBox(
        Part_1_doc,
        "channelLength",
        "channelWidth",
        "channelDepth"
    )

### Create LinearCopy
LinearCopy_1 =\
    model.addMultiTranslation(
        Part_1_doc,
        [model.selection("SOLID", "Box_1_1")],
        model.selection("EDGE", "PartSet/OY"),
        "step",
        "NbStep + 1",
        keepSubResults = True
    )

### Create Box
Box_2 =\
    model.addBox(
        Part_1_doc,
        "channelWidth",
        "channelLength",
        "channelDepth"
    )

### Create LinearCopy
LinearCopy_2 =\
    model.addMultiTranslation(
        Part_1_doc,
        [model.selection("SOLID", "Box_2_1")],
        model.selection("EDGE", "PartSet/OX"),
        "step",
        "NbStep + 1",
        keepSubResults = True
    )

### Create Box
Box_3 =\
    model.addBox(
        Part_1_doc,
        "drainageLength",
        "channelWidth",
        "channelDepth"
    )

### Create LinearCopy
LinearCopy_3 =\
    model.addMultiTranslation(
        Part_1_doc,
        [model.selection("SOLID", "Box_3_1")],
        model.selection("EDGE", "PartSet/OY"),
        "step",
        "NbStep + 1",
        keepSubResults = True
    )

### Create Copy
Copy_1 =\
    model.addCopy(
        Part_1_doc,
        [model.selection("COMPOUND", "LinearCopy_3_1")],
        4
    )

### Create Translation
Translation_1 =\
    model.addTranslation(
        Part_1_doc,
        [model.selection("COMPOUND", "LinearCopy_3_1_1")],
        vector = ["-drainageLength", 0, 0],
        keepSubResults = True
    )

### Create Translation
Translation_2 =\
    model.addTranslation(
        Part_1_doc,
        [model.selection("COMPOUND", "LinearCopy_3_1_2")],
        vector = ["channelLength", 0, 0],
        keepSubResults = True
    )

### Create Rotation
Rotation_1 =\
    model.addRotation(
        Part_1_doc,
        [model.selection("COMPOUND", "LinearCopy_3_1_3")],
        axis = model.selection("EDGE", "PartSet/OZ"),
        angle = 270,
        keepSubResults = True
    )

### Create Rotation
Rotation_2 =\
    model.addRotation(
        Part_1_doc,
        [model.selection("COMPOUND", "LinearCopy_3_1_4")],
        axis = model.selection("EDGE", "PartSet/OZ"),
        angle = 90,
        keepSubResults = True
    )

### Create Translation
Translation_3 =\
    model.addTranslation(
        Part_1_doc,
        [model.selection("COMPOUND", "LinearCopy_3_1_4")],
        vector = ["channelLength", "channelLength", 0],
        keepSubResults = True
    )

### Create Fuse
Fuse_1 =\
    model.addFuse(
        Part_1_doc,
        [
            model.selection("COMPOUND", "LinearCopy_1_1"),
            model.selection("COMPOUND", "LinearCopy_2_1"),
            model.selection("COMPOUND", "LinearCopy_3_1_1"),
            model.selection("COMPOUND", "LinearCopy_3_1_2"),
            model.selection("COMPOUND", "LinearCopy_3_1_3"),
            model.selection("COMPOUND", "LinearCopy_3_1_4")
        ],
        removeEdges = True,
        keepSubResults = True
    )

### Create Translation
Translation_4 =\
    model.addTranslation(
        Part_1_doc,
        [model.selection("SOLID", "Fuse_1_1")],
        vector = ["-0.5*channelLength", "-0.5*channelLength", "-0.5*channelDepth"],
        keepSubResults = True
    )

### Create Group
Group_1 =\
    model.addGroup(
        Part_1_doc,
        "Faces",
        [
            model.selection("SOLID", "Translation_4_1"),
            model.filters(
                Part_1_doc,
                [model.addFilter(name = "ExternalFaces")]
            )
        ]
    )

Group_1.setName("allFaces")
Group_1.result().setName("allFaces")

### Create Plane
Plane_4 = model.addPlane(Part_1_doc, model.selection("FACE", "PartSet/XOY"), "0.5*channelDepth", False)

### Create Group
Group_2 =\
    model.addGroup(
        Part_1_doc,
        "Faces",
        [
            model.selection("SOLID", "Translation_4_1"),
            model.filters(
                Part_1_doc,
                [model.addFilter(name = "OnGeometry", args = [Plane_4.result()])]
            )
        ]
    )
Group_2.setName("top")
Group_2.result().setName("top")

#### Create Group
#Group_2 =\
    #model.addGroup(
        #Part_1_doc,
        #"Faces",
        #[
            #model.selection("SOLID", "Translation_4_1"),
            #model.filters(Part_1_doc, [model.addFilter(name = "ExternalFaces")])
        #]
    #)

#Group_2.setName("top")
#Group_2.result().setName("top")

### Create Plane
Plane_5 = model.addPlane(Part_1_doc, model.selection("FACE", "PartSet/XOY"), "0.5*channelDepth", True)

model.end()

###
### SHAPERSTUDY component
###

#model.publishToShaperStudy()
#import SHAPERSTUDY

#Translation_4_1, = SHAPERSTUDY.shape(model.featureStringId(Translation_4))

####
#### SMESH component
####

#import  SMESH, SALOMEDS
#from salome.smesh import smeshBuilder

#smesh = smeshBuilder.New()
##smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 ## multiples meshes built in parallel, complex and numerous mesh edition (performance)

#lattice = smesh.Mesh(Translation_4_1)
#NETGEN_1D_2D = lattice.Triangle(algo=smeshBuilder.NETGEN_1D2D)
#NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
#NETGEN_2D_Parameters_1.SetMaxSize(
        #min(
            #channelLength.value(),
            #channelWidth.value(),
            #channelDepth.value()
        #)
    #)
#NETGEN_2D_Parameters_1.SetMinSize(
        #min(
            #channelLength.value(),
            #channelWidth.value(),
            #channelDepth.value()
        #)
    #)
#NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
#NETGEN_2D_Parameters_1.SetOptimize( 1 )
#NETGEN_2D_Parameters_1.SetFineness( 5 )
#NETGEN_2D_Parameters_1.SetGrowthRate( 0.1 )
#NETGEN_2D_Parameters_1.SetNbSegPerEdge( 3 )
#NETGEN_2D_Parameters_1.SetNbSegPerRadius( 5 )
#NETGEN_2D_Parameters_1.SetChordalError( -1 )
#NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
#NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
#NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
#NETGEN_2D_Parameters_1.SetWorstElemMeasure( 21874 )
#NETGEN_2D_Parameters_1.SetUseDelauney( 0 )
#NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
#NETGEN_2D_Parameters_1.SetCheckChartBoundary( 184 )
#isDone = lattice.Compute()

#try:
  #lattice.ExportSTL( r'/media/alexshtil/STORAGE/calculations/OpenFOAM/regularPorousMedia/twoPhase/incompressible/lattice/src/lattice.stl', 1 )
  #pass
#except:
  #print('ExportSTL() failed. Invalid file name?')

### Set names of Mesh objects
#smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
#smesh.SetName(lattice.GetMesh(), 'lattice')
#smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')

#salome.myStudy.SaveAs("/media/alexshtil/STORAGE/calculations/OpenFOAM/regularPorousMedia/twoPhase/incompressible/lattice/src/lattice.hdf", False, False)


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
