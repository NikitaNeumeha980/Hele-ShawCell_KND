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

tolerance = 1e-07

sectorAngle = 2.0

#For debuging
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

middleLineOnBottom = geompy.MakeLineTwoPnt(faceVertexOnBottom_1, faceVertexOnBottom_2)
middleLineOnTop = geompy.MakeLineTwoPnt(faceVertexOnTop_1, faceVertexOnTop_2)

bottomFace = geompy.MakeRevolution2Ways(middleLineOnBottom, OZ, 0.5*sectorAngle*math.pi/180.0)
topFace = geompy.MakeRevolution2Ways(middleLineOnTop, OZ, 0.5*sectorAngle*math.pi/180.0)

middleFace = geompy.MakeQuad4Vertices(faceVertexOnBottom_1, faceVertexOnBottom_2, faceVertexOnTop_2, faceVertexOnTop_1)
leftFace = geompy.MakeRotation(middleFace, OZ, 0.5*sectorAngle*math.pi/180.0)
rightFace = geompy.MakeRotation(middleFace, OZ, -0.5*sectorAngle*math.pi/180.0)

#For debuging
#geompy.addToStudy( bottomFace, 'bottomFace' )
#geompy.addToStudy( topFace, 'topFace' )
#geompy.addToStudy( middleFace, 'middleFace' )
#geompy.addToStudy( leftFace, 'leftFace' )
#geompy.addToStudy( rightFace, 'rightFace' )

#Make channel profile
channelProfile = geompy.MakeCurveParametric("0.5 * " + str(Rmax + Rmin) + "+ 0.5 * " + str(Rmax - Rmin) + "* cos(2 * pi * t / " + str(channelLength) +")", "0", "t", 0, channelLength, 20, GEOM.Interpolation, True)

channelProfileFace = geompy.MakeRevolution2Ways(channelProfile, OZ, 0.5*sectorAngle*math.pi/180.0)

sk = geompy.Sketcher3D()
sk.addPointsAbsolute(Rmax, 0, 0)
sk.addPointsAbsolute(0, 0, 0)
sk.addPointsAbsolute(0, 0, channelLength)
sk.addPointsAbsolute(Rmax, 0, channelLength)
closingContour = sk.wire()

channel2D = geompy.MakeFaceWires([channelProfile, closingContour], 1)

geompy.addToStudy( channel2D, 'channel2DRmax' + str(Rmax) + 'channelLength' + str(channelLength) )

geompy.addToStudyInFather( channel2D, bottomFace, 'bottomFace' )
geompy.addToStudyInFather( channel2D, topFace, 'topFace' )
geompy.addToStudyInFather( channel2D, middleFace, 'middleFace' )
geompy.addToStudyInFather( channel2D, leftFace, 'leftFace' )
geompy.addToStudyInFather( channel2D, rightFace, 'rightFace' )
geompy.addToStudyInFather( channel2D, channelProfileFace, 'channelProfileFace' )

channel2DchannelProfile = geompy.GetEdge(channel2D, faceVertexOnBottom_2, faceVertexOnTop_2)
channel2DEdgeOnBottom = geompy.GetEdge(channel2D, faceVertexOnBottom_1, faceVertexOnBottom_2)
channel2DEdgeOnTop = geompy.GetEdge(channel2D, faceVertexOnTop_1, faceVertexOnTop_2)

geompy.addToStudyInFather( channel2D, channel2DchannelProfile, 'channel2DchannelProfile' )
geompy.addToStudyInFather( channel2D, channel2DEdgeOnBottom, 'channel2DEdgeOnBottom' )
geompy.addToStudyInFather( channel2D, channel2DEdgeOnTop, 'channel2DEdgeOnTop' )

##Save geometry to file
print("Save geometry to " + ABSOLUTE_CASE_PATH + "/channelGeometry.hdf\n")

salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/channelGeometry.hdf", salome.myStudy, False)

###
### SMESH component
###

smesh = smeshBuilder.New(theStudy)
channelMesh = smesh.Mesh(channel2D)

globalRegular1D = channelMesh.Segment()
globalNumberOfSegments = globalRegular1D.NumberOfSegments(200)

if (
    channelLength >= 10.0 * Rmin and
    Rmax / Rmin <= 3.0
):
    quadFromMedialAxis1D2D = channelMesh.Quadrangle(algo=smeshBuilder.QUAD_MA_PROJ)
    numberOfLayersForQuadFromMedialAxis1D2D = quadFromMedialAxis1D2D.NumberOfLayers(50)
else:
    regular1D = channelMesh.Segment(geom=channel2DchannelProfile)
    numberOfSegmentsForRegular1D = regular1D.NumberOfSegments(200,None,[])
    numberOfSegmentsForRegular1D.SetConversionMode( 1 )
    numberOfSegmentsForRegular1D.SetExpressionFunction( '(t*(1-t))^0.2 + 5.0*exp(-150.0*(t-0.5)^2)' )
    
    quadrangle2D = channelMesh.Quadrangle(algo=smeshBuilder.QUADRANGLE)
    quadrangleParameters = quadrangle2D.QuadrangleParameters(StdMeshersBuilder.QUAD_QUADRANGLE_PREF,-1,[],[])

subMeshOnChannel2DEdgeOnBottom = channelMesh.Segment(geom=channel2DEdgeOnBottom)
startAndEndLengthOnChannel2DEdgeOnBottom = subMeshOnChannel2DEdgeOnBottom.StartEndLength(0.01*Rmax,0.1*Rmax,[])

subMeshOnChannel2DEdgeOnTop = channelMesh.Segment(geom=channel2DEdgeOnTop)
startAndEndLengthOnChannel2DEdgeOnTop = subMeshOnChannel2DEdgeOnTop.StartEndLength(0.1*Rmax,0.01*Rmax,[])

isDone = channelMesh.Compute()

channelMesh.RotationSweepObjects( [ channelMesh ], [ channelMesh ], [ channelMesh ], SMESH.AxisStruct( 0, 0, 0, 0, 0, 1 ), sectorAngle*math.pi/180.0, 1, tolerance, 0 )
channelMesh.RotateObject( channelMesh, SMESH.AxisStruct( 0, 0, 0, 0, 0, 1 ), -0.5*sectorAngle*math.pi/180.0, 0 )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,bottomFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(channelMesh.GetMesh())
inlet = channelMesh.GroupOnFilter( SMESH.FACE, 'inlet', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,topFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(channelMesh.GetMesh())
outlet = channelMesh.GroupOnFilter( SMESH.FACE, 'outlet', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,leftFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(channelMesh.GetMesh())
wedge_1 = channelMesh.GroupOnFilter( SMESH.FACE, 'wedge_1', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,rightFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(channelMesh.GetMesh())
wedge_2 = channelMesh.GroupOnFilter( SMESH.FACE, 'wedge_2', aFilter )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,channelProfileFace)
aCriteria.append(aCriterion)
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(channelMesh.GetMesh())
channelWall= channelMesh.GroupOnFilter( SMESH.FACE, 'channelWall', aFilter )

smesh.SetName(channelMesh.GetMesh(), 'channelMesh' + str(Rmax) + 'channelLength' + str(channelLength))

channelMesh.ExportUNV( ABSOLUTE_CASE_PATH + "/channel.unv" )

##Save geometry and mesh to file
print("Save geometry and mesh to " + ABSOLUTE_CASE_PATH + "/channelGeometryAndMesh.hdf\n")

salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/channelGeometryAndMesh.hdf", salome.myStudy, False)

print("done")
print("***********************************\n")

#if salome.sg.hasDesktop():
  #salome.sg.updateObjBrowser(True)
