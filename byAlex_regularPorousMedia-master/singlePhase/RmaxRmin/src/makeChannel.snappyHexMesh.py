import sys
import math
import numpy
import salome

import GEOM
from salome.geom import geomBuilder
import SALOMEDS

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder
from salome.StdMeshers import StdMeshersBuilder

salome.salome_init()
theStudy = salome.myStudy

####################################################
##            Begin of variables section          ##
####################################################
ABSOLUTE_CASE_PATH = sys.argv[1]

Rmin = float(sys.argv[2])
Rmax = float(sys.argv[3])
channelLength = float(sys.argv[4])

numberOfSegmentsAlongChannelLength0 = 20

backgroundMeshBoxLength = channelLength
backgroundMeshBoxWidth = 1.1 * Rmax

delta = min(backgroundMeshBoxWidth, channelLength / numberOfSegmentsAlongChannelLength0)
numberOfSegmentsAlongChannelLength = int(channelLength / delta)

backgroundMeshBoxDepth = channelLength / numberOfSegmentsAlongChannelLength

tolerance = 1e-07

print("\n***********************************")
print("ABSOLUTE CASE PATH: " + ABSOLUTE_CASE_PATH)
print("Rmin: " + str(Rmin))
print("Rmax: " + str(Rmax))
print("channelLength: " + str(channelLength))
####################################################
##             End of variables section           ##
####################################################

###
### GEOM component
###

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

#Make background mesh box
backgroundMeshBoxVertexOnBottom_1 = geompy.MakeVertex(
    0,
    0,
    0
)

backgroundMeshBoxVertexOnBottom_2 = geompy.MakeVertex(
    backgroundMeshBoxWidth,
    0,
    0
)

backgroundMeshBoxVertexOnBottom_3 = geompy.MakeVertex(
    backgroundMeshBoxWidth,
    backgroundMeshBoxDepth,
    0
)

backgroundMeshBoxVertexOnBottom_4 = geompy.MakeVertex(
    0,
    backgroundMeshBoxDepth,
    0
)

backgroundMeshBoxVertexOnTop_1 = geompy.MakeVertex(
    0,
    0,
    backgroundMeshBoxLength
)

backgroundMeshBoxVertexOnTop_2 = geompy.MakeVertex(
    backgroundMeshBoxWidth,
    0,
    backgroundMeshBoxLength
)

backgroundMeshBoxVertexOnTop_3 = geompy.MakeVertex(
    backgroundMeshBoxWidth,
    backgroundMeshBoxDepth,
    backgroundMeshBoxLength
)

backgroundMeshBoxVertexOnTop_4 = geompy.MakeVertex(
    0,
    backgroundMeshBoxDepth,
    backgroundMeshBoxLength
)

#Make faces of porous volume
backgroundMeshBoxLeftFace = geompy.MakeQuad4Vertices(backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_4, backgroundMeshBoxVertexOnTop_4, backgroundMeshBoxVertexOnTop_1)
backgroundMeshBoxRightFace = geompy.MakeQuad4Vertices(backgroundMeshBoxVertexOnBottom_2, backgroundMeshBoxVertexOnBottom_3, backgroundMeshBoxVertexOnTop_3, backgroundMeshBoxVertexOnTop_2)
backgroundMeshBoxOppositeFace = geompy.MakeQuad4Vertices(backgroundMeshBoxVertexOnBottom_4, backgroundMeshBoxVertexOnBottom_3, backgroundMeshBoxVertexOnTop_3, backgroundMeshBoxVertexOnTop_4)
backgroundMeshBoxFrontFace = geompy.MakeQuad4Vertices(backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_2, backgroundMeshBoxVertexOnTop_2, backgroundMeshBoxVertexOnTop_1)
backgroundMeshBoxTopFace = geompy.MakeQuad4Vertices(backgroundMeshBoxVertexOnTop_1, backgroundMeshBoxVertexOnTop_2, backgroundMeshBoxVertexOnTop_3, backgroundMeshBoxVertexOnTop_4)
backgroundMeshBoxBottomFace = geompy.MakeQuad4Vertices(backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_2, backgroundMeshBoxVertexOnBottom_3, backgroundMeshBoxVertexOnBottom_4)

backgroundMeshBox = geompy.MakeHexa(backgroundMeshBoxLeftFace, backgroundMeshBoxRightFace, backgroundMeshBoxOppositeFace, backgroundMeshBoxFrontFace, backgroundMeshBoxTopFace, backgroundMeshBoxBottomFace)

edgeForSubMeshOnDepth = geompy.GetEdge(backgroundMeshBox, backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_4)
edgeForSubMeshOnWidth = geompy.GetEdge(backgroundMeshBox, backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_2)

faceVertexOnBottom_1 = geompy.MakeVertex(
    0,
    0,
    0
)
faceVertexOnBottom_2 = geompy.MakeVertex(
    Rmax,
    0,
    0
)

faceVertexOnTop_1 = geompy.MakeVertex(
    0,
    0,
    channelLength
)
faceVertexOnTop_2 = geompy.MakeVertex(
    Rmax,
    0,
    channelLength
)

#Make channel profile
channelProfile1D = geompy.MakeCurveParametric("0.5 *" + str(Rmax + Rmin) + "+ 0.5 *" + str(Rmax - Rmin) + "* cos(2 * pi * t /" + str(channelLength) + ")", "0", "t", 0, channelLength, 100, GEOM.Interpolation, True)

channelWall = geompy.MakePrismVecH(channelProfile1D, OY, backgroundMeshBoxDepth)

sk = geompy.Sketcher3D()
sk.addPointsAbsolute(Rmax, 0, 0)
sk.addPointsAbsolute(0, 0, 0)
sk.addPointsAbsolute(0, 0, channelLength)
sk.addPointsAbsolute(Rmax, 0, channelLength)
closingContour = sk.wire()

channel2D = geompy.MakeFaceWires([channelProfile1D, closingContour], 1)

planeChannel3D = geompy.MakePrismVecH2Ways(channel2D, OY, 2.0 * backgroundMeshBoxDepth)

geompy.addToStudy( planeChannel3D, 'planeChannel3DRmax' + str(Rmax) + 'channelLength' + str(channelLength) )

geompy.addToStudyInFather( planeChannel3D, channelWall, 'channelWall' )

geompy.addToStudyInFather( planeChannel3D, channel2D, 'channel2D' )

geompy.addToStudyInFather( planeChannel3D, backgroundMeshBox, 'backgroundMeshBox' )

geompy.addToStudyInFather( backgroundMeshBox, edgeForSubMeshOnDepth, 'edgeForSubMeshOnDepth' )
geompy.addToStudyInFather( backgroundMeshBox, edgeForSubMeshOnWidth, 'edgeForSubMeshOnWidth' )

geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxBottomFace, 'backgroundMeshBoxBottomFace' )
geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxTopFace, 'backgroundMeshBoxTopFace' )
geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxOppositeFace, 'backgroundMeshBoxOppositeFace' )
geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxFrontFace, 'backgroundMeshBoxFrontFace' )

##Save geometry to file
geompy.ExportSTL(channelWall, ABSOLUTE_CASE_PATH + "/channelWall.stl", True, 1e-4, True)
geompy.ExportSTL(planeChannel3D, ABSOLUTE_CASE_PATH + "/planeChannel3D.stl", True, 1e-4, True)

print("Save geometry to " + ABSOLUTE_CASE_PATH + "/channelGeometry.hdf\n")

salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/channelGeometry.hdf", salome.myStudy, False)

#locationInMesh output
locationInMeshFile = open(ABSOLUTE_CASE_PATH + "/locationInMesh", "w", 0)
locationInMeshFile.write("locationInMesh (" + str(0.5 * Rmin) + " " + str(0.5 * backgroundMeshBoxDepth) + " " + str(0.5 * channelLength) + ");")

###
### SMESH component
###

smesh = smeshBuilder.New(theStudy)
backgroundMeshBoxMesh = smesh.Mesh(backgroundMeshBox)

globalRegular1D = backgroundMeshBoxMesh.Segment()
globalNumberOfSegments = globalRegular1D.NumberOfSegments(numberOfSegmentsAlongChannelLength)

quadrangle2D = backgroundMeshBoxMesh.Quadrangle(algo=smeshBuilder.QUADRANGLE)
quadrangleParameters = quadrangle2D.QuadrangleParameters(StdMeshersBuilder.QUAD_QUADRANGLE_PREF,-1,[],[])

hexa3D = backgroundMeshBoxMesh.Hexahedron(algo=smeshBuilder.Hexa)

#Sub-mesh on width
regular1DOnWidth = backgroundMeshBoxMesh.Segment(geom=edgeForSubMeshOnWidth)
numberOfSegmentsOnWidth = regular1DOnWidth.NumberOfSegments(
    max(
        [
            int(round(backgroundMeshBoxWidth / backgroundMeshBoxLength * numberOfSegmentsAlongChannelLength)),
            1
        ]
    )
)
propagationOf1DHypOnWidth = regular1DOnWidth.Propagation()

#Sub-mesh on depth
regular1DOnDepth = backgroundMeshBoxMesh.Segment(geom=edgeForSubMeshOnDepth)
numberOfSegmentsOnDepth = regular1DOnDepth.NumberOfSegments(1)
propagationOf1DHypOnDepth = regular1DOnDepth.Propagation()

isDone = backgroundMeshBoxMesh.Compute()

smesh.SetName(backgroundMeshBoxMesh.GetMesh(), 'backgroundMeshBoxMeshRmax' + str(Rmax) + 'channelLength' + str(channelLength))

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,backgroundMeshBoxBottomFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(backgroundMeshBoxMesh.GetMesh())
inlet = backgroundMeshBoxMesh.GroupOnFilter( SMESH.FACE, 'inlet', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,backgroundMeshBoxTopFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(backgroundMeshBoxMesh.GetMesh())
outlet = backgroundMeshBoxMesh.GroupOnFilter( SMESH.FACE, 'outlet', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,backgroundMeshBoxFrontFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(backgroundMeshBoxMesh.GetMesh())
wedge_1 = backgroundMeshBoxMesh.GroupOnFilter( SMESH.FACE, 'wedge_1', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,backgroundMeshBoxOppositeFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(backgroundMeshBoxMesh.GetMesh())
wedge_2 = backgroundMeshBoxMesh.GroupOnFilter( SMESH.FACE, 'wedge_2', aFilter )

backgroundMeshBoxMesh.ExportUNV( ABSOLUTE_CASE_PATH + "/backgroundMesh.unv" )

##Save geometry and mesh to file
print("Save geometry and mesh to " + ABSOLUTE_CASE_PATH + "/channelGeometryAndBackgroundMesh.hdf\n")

salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/channelGeometryAndBackgroundMesh.hdf", salome.myStudy, False)

print("done")
print("***********************************\n")

