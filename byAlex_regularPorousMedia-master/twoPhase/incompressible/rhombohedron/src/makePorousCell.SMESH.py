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
#ABSOLUTE_CASE_PATH = sys.argv[0]
#ABSOLUTE_CASE_PATH = ABSOLUTE_CASE_PATH[0:ABSOLUTE_CASE_PATH.rfind("/")]
#ABSOLUTE_CASE_PATH = ABSOLUTE_CASE_PATH[0:ABSOLUTE_CASE_PATH.rfind("/")]

#For debuging
print("\n***********************************")
print("ABSOLUTE CASE PATH: " + ABSOLUTE_CASE_PATH)

theta = float(sys.argv[2])
intersectionParameter = float(sys.argv[3])

theta = math.radians(theta)
alpha = math.acos(math.cos(theta) / math.cos(0.5 * theta))

#scaleFactor = float(sys.argv[4])

#increasingScaleFactor = float(sys.argv[5])

#UTranslationDirectionNbTimes = int(sys.argv[6])
#VTranslationDirectionNbTimes = int(sys.argv[7])
#WTranslationDirectionNbTimes = int(sys.argv[8])

increasingScaleFactor = float(sys.argv[4])

UTranslationDirectionNbTimes = int(sys.argv[5])
VTranslationDirectionNbTimes = int(sys.argv[6])
WTranslationDirectionNbTimes = int(sys.argv[7])

onePoreCellSize = 1.0 * increasingScaleFactor
grainSize0 = 0.5 * increasingScaleFactor
filletFactor = 0.15

tolerance = 1e-07
####################################################
##             End of variables section           ##
####################################################

###
### GEOM component
###

geompy = geomBuilder.New(theStudy)

#For debuging
print("theta=" + str(math.degrees(theta)))
print("alpha=" + str(math.degrees(math.acos(math.cos(theta) / math.cos(0.5 * theta)))))
print("intersectionParameter=" + str(intersectionParameter))

grainSize = grainSize0 / (1.0 - intersectionParameter)

#Building the sceleton
#Make the grain's centeres
grainCentreOnBottom_1 = geompy.MakeVertex(
    0 * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(0)
)
grainCentreOnBottom_2 = geompy.MakeVertex(
    onePoreCellSize * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(0)
)
grainCentreOnBottom_3 = geompy.MakeVertex(
    onePoreCellSize + onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(0)
)
grainCentreOnBottom_4 = geompy.MakeVertex(
    onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(0)
)

grainCentreOnTop_1 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + 0 * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(alpha)
)
grainCentreOnTop_2 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + onePoreCellSize * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(alpha)
)
grainCentreOnTop_3 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + onePoreCellSize + onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(alpha)
)
grainCentreOnTop_4 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(alpha)
)

#Make the grains with its rotations
#Grain #1
grain_1 = geompy.MakeSpherePntR(grainCentreOnBottom_1, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + math.pi), onePoreCellSize * math.sin(0.5 * theta + math.pi), 0)
grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * theta + math.pi)
grain_1 = geompy.MakeRotation(grain_1, rotationVector_2, 0.5 * math.pi)

grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * math.pi)
mirroredGrain_1 = geompy.MakeRotation(grain_1, rotationVector_1, math.pi)

#Grain #2
grain_2 = geompy.MakeSpherePntR(grainCentreOnBottom_2, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 1.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnBottom_2)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnBottom_2)
grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
grain_2 = geompy.MakeRotation(grain_2, rotationVector_2, 0.5 * math.pi)

grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * math.pi)
mirroredGrain_2 = geompy.MakeRotation(grain_2, rotationVector_1, math.pi)

#Grain #3
grain_3 = geompy.MakeSpherePntR(grainCentreOnBottom_3, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta), onePoreCellSize * math.sin(0.5 * theta), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnBottom_3)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnBottom_3)
grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * theta)
grain_3 = geompy.MakeRotation(grain_3, rotationVector_2, 0.5 * math.pi)

grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * math.pi)
mirroredGrain_3 = geompy.MakeRotation(grain_3, rotationVector_1, math.pi)

#Grain #4
grain_4 = geompy.MakeSpherePntR(grainCentreOnBottom_4, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 0.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnBottom_4)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnBottom_4)
grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
grain_4 = geompy.MakeRotation(grain_4, rotationVector_2, 0.5 * math.pi)

grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * math.pi)
mirroredGrain_4 = geompy.MakeRotation(grain_4, rotationVector_1, math.pi)

#Grain #5
grain_5 = geompy.MakeSpherePntR(grainCentreOnTop_1, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + math.pi), onePoreCellSize * math.sin(0.5 * theta + math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnTop_1)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnTop_1)
grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * theta + math.pi)
grain_5 = geompy.MakeRotation(grain_5, rotationVector_2, 0.5 * math.pi)

grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * math.pi)
mirroredGrain_5 = geompy.MakeRotation(grain_5, rotationVector_1, math.pi)

#Grain #6
grain_6 = geompy.MakeSpherePntR(grainCentreOnTop_2, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 1.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnTop_2)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnTop_2)
grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
grain_6 = geompy.MakeRotation(grain_6, rotationVector_2, 0.5 * math.pi)

grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * math.pi)
mirroredGrain_6 = geompy.MakeRotation(grain_6, rotationVector_1, math.pi)

#Grain #7
grain_7 = geompy.MakeSpherePntR(grainCentreOnTop_3, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta), onePoreCellSize * math.sin(0.5 * theta), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnTop_3)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnTop_3)
grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * theta)
grain_7 = geompy.MakeRotation(grain_7, rotationVector_2, 0.5 * math.pi)

grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * math.pi)
mirroredGrain_7 = geompy.MakeRotation(grain_7, rotationVector_1, math.pi)

#Grain #8
grain_8 = geompy.MakeSpherePntR(grainCentreOnTop_4, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 0.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCentreOnBottom_1, grainCentreOnTop_4)
geompy.TranslateTwoPoints(rotationVector_2, grainCentreOnBottom_1, grainCentreOnTop_4)
grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
grain_8 = geompy.MakeRotation(grain_8, rotationVector_2, 0.5 * math.pi)

grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * math.pi)
mirroredGrain_8 = geompy.MakeRotation(grain_8, rotationVector_1, math.pi)

#For debuging
geompy.addToStudy( grain_1, 'grain_1')
geompy.addToStudy( grain_2, 'grain_2')
geompy.addToStudy( grain_3, 'grain_3')
geompy.addToStudy( grain_4, 'grain_4')
geompy.addToStudy( grain_5, 'grain_5')
geompy.addToStudy( grain_6, 'grain_6')
geompy.addToStudy( grain_7, 'grain_7')
geompy.addToStudy( grain_8, 'grain_8')
geompy.addToStudy( mirroredGrain_1, 'mirroredGrain_1')
geompy.addToStudy( mirroredGrain_2, 'mirroredGrain_2')
geompy.addToStudy( mirroredGrain_3, 'mirroredGrain_3')
geompy.addToStudy( mirroredGrain_4, 'mirroredGrain_4')
geompy.addToStudy( mirroredGrain_5, 'mirroredGrain_5')
geompy.addToStudy( mirroredGrain_6, 'mirroredGrain_6')
geompy.addToStudy( mirroredGrain_7, 'mirroredGrain_7')
geompy.addToStudy( mirroredGrain_8, 'mirroredGrain_8')

#building the porous volume
#Make faces of porous volume
onePoreBoxLeftFace = geompy.MakeQuad4Vertices(grainCentreOnBottom_1, grainCentreOnBottom_4, grainCentreOnTop_4, grainCentreOnTop_1)
onePoreBoxRightFace = geompy.MakeQuad4Vertices(grainCentreOnBottom_2, grainCentreOnBottom_3, grainCentreOnTop_3, grainCentreOnTop_2)
onePoreBoxOppositeFace = geompy.MakeQuad4Vertices(grainCentreOnBottom_4, grainCentreOnBottom_3, grainCentreOnTop_3, grainCentreOnTop_4)
onePoreBoxFrontFace = geompy.MakeQuad4Vertices(grainCentreOnBottom_1, grainCentreOnBottom_2, grainCentreOnTop_2, grainCentreOnTop_1)
onePoreBoxTopFace = geompy.MakeQuad4Vertices(grainCentreOnTop_1, grainCentreOnTop_2, grainCentreOnTop_3, grainCentreOnTop_4)
onePoreBoxBottomFace = geompy.MakeQuad4Vertices(grainCentreOnBottom_1, grainCentreOnBottom_2, grainCentreOnBottom_3, grainCentreOnBottom_4)

onePoreBox = geompy.MakeHexa(onePoreBoxLeftFace, onePoreBoxRightFace, onePoreBoxOppositeFace, onePoreBoxFrontFace, onePoreBoxTopFace, onePoreBoxBottomFace)

#For debuging
geompy.addToStudy( onePoreBox, 'onePoreBox')

fusedGrains = geompy.MakeFuseList([grain_1, grain_2, grain_3, grain_4, grain_5, grain_6, grain_7, grain_8, mirroredGrain_1, mirroredGrain_2, mirroredGrain_3, mirroredGrain_4, mirroredGrain_5, mirroredGrain_6, mirroredGrain_7, mirroredGrain_8])

#For debuging
geompy.addToStudy( fusedGrains, 'fusedGrains')

onePore = geompy.MakeCut(onePoreBox, fusedGrains)

onePore = geompy.RemoveExtraEdges(onePore, False)

#Get edges IDs with fillet
checkFilletSwitcher = 0

allEdgesList = []
allEdgesWithFillet = []
allEdgesIDsWithFillet = []
allEdgesWithoutFillet = []

allEdgesList = geompy.SubShapeAll(onePore, geompy.ShapeType["EDGE"])
allEdgesWithFillet = geompy.CreateGroup(onePore, geompy.ShapeType["EDGE"])

geompy.UnionList(allEdgesWithFillet, allEdgesList)

#Make middle points on boundaries for searching boundaries
onePoreBoxMiddlePointOnLeftFace = geompy.MakeVertexOnSurface(onePoreBoxLeftFace, 0.5, 0.5)
onePoreBoxMiddlePointOnRightFace = geompy.MakeVertexOnSurface(onePoreBoxRightFace, 0.5, 0.5)
onePoreBoxMiddlePointOnOppositeFace = geompy.MakeVertexOnSurface(onePoreBoxOppositeFace, 0.5, 0.5)
onePoreBoxMiddlePointOnFrontFace = geompy.MakeVertexOnSurface(onePoreBoxFrontFace, 0.5, 0.5)
onePoreBoxMiddlePointOnTopFace = geompy.MakeVertexOnSurface(onePoreBoxTopFace, 0.5, 0.5)
onePoreBoxMiddlePointOnBottomFace = geompy.MakeVertexOnSurface(onePoreBoxBottomFace, 0.5, 0.5)

#Make normals for searching boundaries
onePoreBoxNormalVectorOnLeftFace = geompy.GetNormal(onePoreBoxLeftFace)
onePoreBoxNormalVectorOnRightFace = geompy.GetNormal(onePoreBoxRightFace)
onePoreBoxNormalVectorOnOppositeFace = geompy.GetNormal(onePoreBoxOppositeFace)
onePoreBoxNormalVectorOnFrontFace = geompy.GetNormal(onePoreBoxFrontFace)
onePoreBoxNormalVectorOnTopFace = geompy.GetNormal(onePoreBoxTopFace)
onePoreBoxNormalVectorOnBottomFace = geompy.GetNormal(onePoreBoxBottomFace)

allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnLeftFace, onePoreBoxMiddlePointOnLeftFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnRightFace, onePoreBoxMiddlePointOnRightFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnOppositeFace, onePoreBoxMiddlePointOnOppositeFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnFrontFace, onePoreBoxMiddlePointOnFrontFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnTopFace, onePoreBoxMiddlePointOnTopFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnBottomFace, onePoreBoxMiddlePointOnBottomFace, GEOM.ST_ON))

geompy.DifferenceList(allEdgesWithFillet, allEdgesWithoutFillet)
allEdgesIDsWithFillet = geompy.GetObjectIDs(allEdgesWithFillet)
allEdgesWithFillet = geompy.ExtractShapes(allEdgesWithFillet, geompy.ShapeType["EDGE"], True)

if (
    filletFactor != 0 and
    intersectionParameter != 0
):
    #Make fillet
    print("filletRadius=" + str(filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])) + "\n")

    #if(
        #filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]) > 0.02
    #):
    fusedFilletedGrains = geompy.MakeFilletAll(fusedGrains, filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

    #fusedFilletedGrains = geompy.MakeChamferAll(fusedGrains, filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

    #For debuging
    geompy.addToStudy( fusedFilletedGrains, 'fusedFilletedGrains')

    onePore = geompy.MakeCut(onePoreBox, fusedFilletedGrains)

    onePore = geompy.RemoveExtraEdges(onePore, False)

    checkFilletSwitcher = 1
    #else:
        #print("Can not make a fillet!\n")

#For debuging
geompy.addToStudy( onePore, 'onePore')

#One pore centre of mass calculation\
onePoreCentreOfMassFile = open(ABSOLUTE_CASE_PATH + "/onePoreCentreOfMass", "w", 0)
onePoreCentreOfMass = geompy.MakeCDG(onePore)
#onePoreCentreOfMassFile.write(
    #"centre ("
    #+ str(scaleFactor / increasingScaleFactor * geompy.PointCoordinates(onePoreCentreOfMass)[0]) + " "
    #+ str(scaleFactor / increasingScaleFactor * geompy.PointCoordinates(onePoreCentreOfMass)[1]) + " "
    #+ str(scaleFactor / increasingScaleFactor * geompy.PointCoordinates(onePoreCentreOfMass)[2]) + ");"
#)
onePoreCentreOfMassFile.write(
    "centre ("
    + str(geompy.PointCoordinates(onePoreCentreOfMass)[0]) + " "
    + str(geompy.PointCoordinates(onePoreCentreOfMass)[1]) + " "
    + str(geompy.PointCoordinates(onePoreCentreOfMass)[2]) + ");"
)
onePoreCentreOfMassFile.close()

#Make translation vectors for building porous cell
UTranslationDirectionVector = geompy.MakeVector(grainCentreOnBottom_1, grainCentreOnBottom_2)
VTranslationDirectionVector = geompy.MakeVector(grainCentreOnBottom_1, grainCentreOnBottom_4)
WTranslationDirectionVector = geompy.MakeVector(grainCentreOnBottom_1, grainCentreOnTop_1)

#Building porous cell
porousCell = geompy.MakeMultiTranslation2D(onePore, UTranslationDirectionVector, onePoreCellSize, UTranslationDirectionNbTimes, VTranslationDirectionVector, onePoreCellSize, VTranslationDirectionNbTimes)

porousCell = geompy.MakeMultiTranslation1D(porousCell, WTranslationDirectionVector, onePoreCellSize, WTranslationDirectionNbTimes)

listOfSolids = []
listOfSolids = geompy.SubShapeAllSortedCentres(porousCell, geompy.ShapeType["SOLID"])

porousCell = geompy.MakeFuseList(listOfSolids)

porousCell = geompy.MakePartition(
        [porousCell],
        [
            geompy.MakeMultiTranslation1D(
                geompy.MakeMultiTranslation2D(
                    geompy.MakeScaleTransform(
                        onePoreBox,
                        onePoreCentreOfMass,
                        0.85
                    ),
                    UTranslationDirectionVector,
                    onePoreCellSize,
                    UTranslationDirectionNbTimes,
                    VTranslationDirectionVector,
                    onePoreCellSize,
                    VTranslationDirectionNbTimes
                ),
                WTranslationDirectionVector,
                onePoreCellSize,
                WTranslationDirectionNbTimes
            )
        ],
        [],
        [],
        geompy.ShapeType["SOLID"],
        0,
        [],
        0
    )

listOfSolids = geompy.SubShapeAllSortedCentres(porousCell, geompy.ShapeType["SOLID"])

throats = geompy.CreateGroup(porousCell, geompy.ShapeType["SOLID"])
geompy.UnionList(throats, listOfSolids)

geompy.DifferenceList(
    throats,
    [
        geompy.GetShapesNearPoint(porousCell, vertexIndex, geompy.ShapeType["SOLID"])
        for vertexIndex in
        geompy.SubShapeAllSortedCentres(
            geompy.MakeMultiTranslation1D(
                geompy.MakeMultiTranslation2D(
                    onePoreCentreOfMass,
                    UTranslationDirectionVector,
                    onePoreCellSize,
                    UTranslationDirectionNbTimes,
                    VTranslationDirectionVector,
                    onePoreCellSize,
                    VTranslationDirectionNbTimes
                ),
                WTranslationDirectionVector,
                onePoreCellSize,
                WTranslationDirectionNbTimes
            ),
            geompy.ShapeType["VERTEX"]
        )
    ]
)

#For debuging
geompy.addToStudy( porousCell, 'porousCell')
geompy.addToStudyInFather( porousCell, throats, 'throats' )

#Make the porous cell box's vertexes
porousCellBoxVertexOnBottom_1 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * 0 * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(0)
)
porousCellBoxVertexOnBottom_2 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(0)
)
porousCellBoxVertexOnBottom_3 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize + UTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.cos(theta),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(theta),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(0)
)
porousCellBoxVertexOnBottom_4 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.cos(theta),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(theta),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(0)
)

porousCellBoxVertexOnTop_1 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * 0 * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(alpha)
)
porousCellBoxVertexOnTop_2 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(alpha)
)
porousCellBoxVertexOnTop_3 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * onePoreCellSize + UTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.cos(theta),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(theta),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(alpha)
)
porousCellBoxVertexOnTop_4 = geompy.MakeVertex(
    UTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.cos(theta),
    VTranslationDirectionNbTimes * onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * onePoreCellSize * math.cos(0) * math.sin(theta),
    WTranslationDirectionNbTimes * onePoreCellSize * math.sin(alpha)
)

#Make faces of porous volume
porousCellBoxLeftFace = geompy.MakeQuad4Vertices(porousCellBoxVertexOnBottom_1, porousCellBoxVertexOnBottom_4, porousCellBoxVertexOnTop_4, porousCellBoxVertexOnTop_1)
porousCellBoxRightFace = geompy.MakeQuad4Vertices(porousCellBoxVertexOnBottom_2, porousCellBoxVertexOnBottom_3, porousCellBoxVertexOnTop_3, porousCellBoxVertexOnTop_2)
porousCellBoxOppositeFace = geompy.MakeQuad4Vertices(porousCellBoxVertexOnBottom_4, porousCellBoxVertexOnBottom_3, porousCellBoxVertexOnTop_3, porousCellBoxVertexOnTop_4)
porousCellBoxFrontFace = geompy.MakeQuad4Vertices(porousCellBoxVertexOnBottom_1, porousCellBoxVertexOnBottom_2, porousCellBoxVertexOnTop_2, porousCellBoxVertexOnTop_1)
porousCellBoxTopFace = geompy.MakeQuad4Vertices(porousCellBoxVertexOnTop_1, porousCellBoxVertexOnTop_2, porousCellBoxVertexOnTop_3, porousCellBoxVertexOnTop_4)
porousCellBoxBottomFace = geompy.MakeQuad4Vertices(porousCellBoxVertexOnBottom_1, porousCellBoxVertexOnBottom_2, porousCellBoxVertexOnBottom_3, porousCellBoxVertexOnBottom_4)

#Make middle points on boundaries for searching boundaries
porousCellBoxMiddlePointOnLeftFace = geompy.MakeVertexOnSurface(porousCellBoxLeftFace, 0.5, 0.5)
porousCellBoxMiddlePointOnRightFace = geompy.MakeVertexOnSurface(porousCellBoxRightFace, 0.5, 0.5)
porousCellBoxMiddlePointOnOppositeFace = geompy.MakeVertexOnSurface(porousCellBoxOppositeFace, 0.5, 0.5)
porousCellBoxMiddlePointOnFrontFace = geompy.MakeVertexOnSurface(porousCellBoxFrontFace, 0.5, 0.5)
porousCellBoxMiddlePointOnTopFace = geompy.MakeVertexOnSurface(porousCellBoxTopFace, 0.5, 0.5)
porousCellBoxMiddlePointOnBottomFace = geompy.MakeVertexOnSurface(porousCellBoxBottomFace, 0.5, 0.5)

#Make normals for searching boundaries
porousCellBoxNormalVectorOnLeftFace = geompy.GetNormal(porousCellBoxLeftFace)
porousCellBoxNormalVectorOnRightFace = geompy.GetNormal(porousCellBoxRightFace)
porousCellBoxNormalVectorOnOppositeFace = geompy.GetNormal(porousCellBoxOppositeFace)
porousCellBoxNormalVectorOnFrontFace = geompy.GetNormal(porousCellBoxFrontFace)
porousCellBoxNormalVectorOnTopFace = geompy.GetNormal(porousCellBoxTopFace)
porousCellBoxNormalVectorOnBottomFace = geompy.GetNormal(porousCellBoxBottomFace)

porousCellBox = geompy.MakeHexa(porousCellBoxLeftFace, porousCellBoxRightFace, porousCellBoxOppositeFace, porousCellBoxFrontFace, porousCellBoxTopFace, porousCellBoxBottomFace)

#For debuging
geompy.addToStudy( porousCellBox, 'porousCellBox')

facesListOnLeftSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnLeftFace, porousCellBoxMiddlePointOnLeftFace, GEOM.ST_ON)
facesListOnRightSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnRightFace, porousCellBoxMiddlePointOnRightFace, GEOM.ST_ON)
facesListOnOppositeSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnOppositeFace, porousCellBoxMiddlePointOnOppositeFace, GEOM.ST_ON)
facesListOnFrontSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnFrontFace, porousCellBoxMiddlePointOnFrontFace, GEOM.ST_ON)
facesListOnTopSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnTopFace, porousCellBoxMiddlePointOnTopFace, GEOM.ST_ON)
facesListOnBottomSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnBottomFace, porousCellBoxMiddlePointOnBottomFace, GEOM.ST_ON)

##Get fases IDs without viscous layer
allFacesIDsWithoutViscousLayers = []

allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnLeftFace, porousCellBoxMiddlePointOnLeftFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnRightFace, porousCellBoxMiddlePointOnRightFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnOppositeFace, porousCellBoxMiddlePointOnOppositeFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnFrontFace, porousCellBoxMiddlePointOnFrontFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnTopFace, porousCellBoxMiddlePointOnTopFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnBottomFace, porousCellBoxMiddlePointOnBottomFace, GEOM.ST_ON))

#Create groups
#Get faces of non-skeleton wall
nonSkeletonFaces = []

nonSkeletonFaces.extend(facesListOnLeftSide)
nonSkeletonFaces.extend(facesListOnRightSide)
nonSkeletonFaces.extend(facesListOnOppositeSide)
nonSkeletonFaces.extend(facesListOnFrontSide)
nonSkeletonFaces.extend(facesListOnTopSide)
nonSkeletonFaces.extend(facesListOnBottomSide)

listOfSolids = geompy.SubShapeAllSortedCentres(
    geompy.MakeMultiTranslation1D(
        geompy.MakeMultiTranslation2D(
            geompy.MakeScaleTransform(
                onePoreBox,
                onePoreCentreOfMass,
                0.85
            ),
            UTranslationDirectionVector,
            onePoreCellSize,
            UTranslationDirectionNbTimes,
            VTranslationDirectionVector,
            onePoreCellSize,
            VTranslationDirectionNbTimes
        ),
        WTranslationDirectionVector,
        onePoreCellSize,
        WTranslationDirectionNbTimes
    ),
    geompy.ShapeType["SOLID"]
)

nonSkeletonFaces.extend(
    sum(
        [
            geompy.GetShapesOnShape(
                solidIndex,
                porousCell,
                geompy.ShapeType["FACE"],
                GEOM.ST_ON
            )
            for solidIndex in listOfSolids
        ],
        []
    )
)

##Get faces of skeleton wall
allFacesList = []

allFacesList = geompy.SubShapeAll(porousCell, geompy.ShapeType["FACE"])
skeletonWall = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
geompy.UnionList(skeletonWall, allFacesList)
geompy.DifferenceList(skeletonWall, nonSkeletonFaces)

leftSide = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
rightSide = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
oppositeSide = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
frontSide = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
topSide = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
bottomSide = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])

geompy.UnionList(leftSide, facesListOnLeftSide)
geompy.UnionList(rightSide, facesListOnRightSide)
geompy.UnionList(oppositeSide, facesListOnOppositeSide)
geompy.UnionList(frontSide, facesListOnFrontSide)
geompy.UnionList(topSide, facesListOnTopSide)
geompy.UnionList(bottomSide, facesListOnBottomSide)

geompy.addToStudyInFather( porousCell, skeletonWall, 'skeletonWall' )
geompy.addToStudyInFather( porousCell, leftSide, 'leftSide' )
geompy.addToStudyInFather( porousCell, rightSide, 'rightSide' )
geompy.addToStudyInFather( porousCell, oppositeSide, 'oppositeSide' )
geompy.addToStudyInFather( porousCell, frontSide, 'frontSide' )
geompy.addToStudyInFather( porousCell, topSide, 'topSide' )
geompy.addToStudyInFather( porousCell, bottomSide, 'bottomSide' )

##Save geometry to file
print("Save geometry to " + ABSOLUTE_CASE_PATH + "/porousCellGeometry.hdf\n")

salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/porousCellGeometry.hdf", salome.myStudy, False)

###
### SMESH component
###

####################################################
##         Begin of mesh variables section        ##
####################################################
#globalNETGEN_2D_MaxSize = 0.05 * geompy.BasicProperties(onePore)[1] ** (1.0 / 2.0)
#globalNETGEN_2D_MinSize = 0.001 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
#patchesNETGEN_2D_MaxSize = 0.05 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
#patchesNETGEN_2D_MinSize = 0.001 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
#NETGEN_2D_QuadAllowed = 0
#growthRate2D = 0.1

#globalLocalLengthSize = 0.005 * max([geompy.BasicProperties(allEdgesWithFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithFillet))])
#localLengthSize = 0.01 * max([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])

globalLocalLengthSize = 0.02 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)
localLengthSize = 0.06 * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])

globalNETGEN_2D_MaxSize = 2.5 * localLengthSize
globalNETGEN_2D_MinSize = 1.0 * localLengthSize
throatsNETGEN_2D_MaxSize = 1.0 * localLengthSize
throatsNETGEN_2D_MinSize = 1.0 * localLengthSize
patchesNETGEN_2D_MaxSize = 1.0 * localLengthSize
patchesNETGEN_2D_MinSize = 1.0 * localLengthSize
NETGEN_2D_QuadAllowed = 0
globalGrowthRate2D = 0.02
throatsGrowthRate2D = 1e-4
patchesGrowthRate2D = 1e-4

globalNETGEN_3D_MaxSize = 4.0 * globalLocalLengthSize
globalNETGEN_3D_MinSize = 1.0 * localLengthSize
throatsNETGEN_3D_MaxSize = 2.0 * localLengthSize
throatsNETGEN_3D_MinSize = 1.0 * localLengthSize
globalGrowthRate3D = 0.05
throatsGrowthRate3D = 0.02

#viscousLayersSize3D = 10.0 * 1e-3 * onePoreCellSize
globalViscousLayersSize3D = min(onePoreCellSize / increasingScaleFactor, 0.025 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0))
globalViscousLayersNumber = 3
globalViscousLayersGrowth = 1.0

throatsViscousLayersSize3D = min(onePoreCellSize / increasingScaleFactor, 0.025 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0))
throatsViscousLayersNumber = 3
throatsViscousLayersGrowth = 1.0

####################################################
##          End of mesh variables section         ##
####################################################

smesh = smeshBuilder.New(theStudy)
porousCellMesh = smesh.Mesh(porousCell)

if(
    checkFilletSwitcher != 0
):
    globalNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)

    globalNETGEN1D2DParameters = globalNETGEN1D2D.Parameters()
    globalNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    globalNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    globalNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    globalNETGEN1D2DParameters.SetOptimize( 1 )
    globalNETGEN1D2DParameters.SetFineness( 5 )
    globalNETGEN1D2DParameters.SetGrowthRate( globalGrowthRate2D )
    globalNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    globalNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    globalNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    globalNETGEN1D2DParameters.SetSecondOrder( 0 )
    globalNETGEN1D2DParameters.SetFuseEdges( 1 )

    ##Sub-mesh for leftSide boundary
    #leftSideNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=leftSide)

    #leftSideNETGEN1D2DParameters = leftSideNETGEN1D2D.Parameters()
    #leftSideNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    #leftSideNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    #leftSideNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    #leftSideNETGEN1D2DParameters.SetOptimize( 1 )
    #leftSideNETGEN1D2DParameters.SetFineness( 5 )
    #leftSideNETGEN1D2DParameters.SetGrowthRate( growthRate2D )
    ##leftSideNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    ##leftSideNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    #leftSideNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    #leftSideNETGEN1D2DParameters.SetSecondOrder( 0 )
    #leftSideNETGEN1D2DParameters.SetFuseEdges( 1 )

    ##Sub-mesh for oppositeSide boundary
    #oppositeSideNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=oppositeSide)

    #oppositeSideNETGEN1D2DParameters = oppositeSideNETGEN1D2D.Parameters()
    #oppositeSideNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    #oppositeSideNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    #oppositeSideNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    #oppositeSideNETGEN1D2DParameters.SetOptimize( 1 )
    #oppositeSideNETGEN1D2DParameters.SetFineness( 5 )
    #oppositeSideNETGEN1D2DParameters.SetGrowthRate( growthRate2D )
    ##oppositeSideNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    ##oppositeSideNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    #oppositeSideNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    #oppositeSideNETGEN1D2DParameters.SetSecondOrder( 0 )
    #oppositeSideNETGEN1D2DParameters.SetFuseEdges( 1 )

    ##Sub-mesh for topSide boundary
    #topSideNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=topSide)

    #topSideNETGEN1D2DParameters = topSideNETGEN1D2D.Parameters()
    #topSideNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    #topSideNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    #topSideNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    #topSideNETGEN1D2DParameters.SetOptimize( 1 )
    #topSideNETGEN1D2DParameters.SetFineness( 5 )
    #topSideNETGEN1D2DParameters.SetGrowthRate( growthRate2D )
    ##topSideNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    ##topSideNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    #topSideNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    #topSideNETGEN1D2DParameters.SetSecondOrder( 0 )
    #topSideNETGEN1D2DParameters.SetFuseEdges( 1 )
else:
    globalRegular1D = porousCellMesh.Segment()

    globalLocalLength = globalRegular1D.LocalLength( globalNETGEN_2D_MinSize, None, 1e-6 )

    #globalFixedPoints = globalRegular1D.FixedPoints1D([ 0.1, 0.9 ],[ 10, 50, 10 ],[])
    #globalFixedPoints.SetObjectEntry( "porousCell" )

    #globalDeflection = globalRegular1D.Deflection1D(4e-5)

    globalNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D)

    globalNETGEN2DParameters = globalNETGEN2D.Parameters()
    globalNETGEN2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    globalNETGEN2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    globalNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    globalNETGEN2DParameters.SetOptimize( 1 )
    globalNETGEN2DParameters.SetFineness( 5 )
    globalNETGEN2DParameters.SetGrowthRate( globalGrowthRate2D )
    #globalNETGEN2DParameters.SetNbSegPerEdge( 2 )
    #globalNETGEN2DParameters.SetNbSegPerRadius( 3 )
    globalNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
    globalNETGEN2DParameters.SetSecondOrder( 0 )
    globalNETGEN2DParameters.SetFuseEdges( 1 )

#Sub-mesh for leftSide boundary
leftSideRegular1D = porousCellMesh.Segment(geom=leftSide)

leftSideLocalLength = leftSideRegular1D.LocalLength(patchesNETGEN_2D_MinSize, None, 1e-6)
#leftSideDeflection = Regular_1D_1.Deflection1D(1e-5)

leftSideNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=leftSide)

leftSideNETGEN2DParameters = leftSideNETGEN2D.Parameters()
leftSideNETGEN2DParameters.SetMaxSize( patchesNETGEN_2D_MaxSize )
leftSideNETGEN2DParameters.SetMinSize( patchesNETGEN_2D_MinSize )
leftSideNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
leftSideNETGEN2DParameters.SetOptimize( 1 )
leftSideNETGEN2DParameters.SetFineness( 5 )
leftSideNETGEN2DParameters.SetGrowthRate( patchesGrowthRate2D )
leftSideNETGEN2DParameters.SetNbSegPerEdge( 2 )
leftSideNETGEN2DParameters.SetNbSegPerRadius( 3 )
leftSideNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
leftSideNETGEN2DParameters.SetSecondOrder( 0 )
leftSideNETGEN2DParameters.SetFuseEdges( 1 )

#Sub-mesh for oppositeSide boundary
oppositeSideRegular1D = porousCellMesh.Segment(geom=oppositeSide)

#oppositeSideDeflection = oppositeSideRegular1D.Deflection1D(1e-5)
oppositeSideLocalLength = oppositeSideRegular1D.LocalLength(patchesNETGEN_2D_MinSize, None, 1e-6)

oppositeSideNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=oppositeSide)

oppositeSideNETGEN2DParameters = oppositeSideNETGEN2D.Parameters()
oppositeSideNETGEN2DParameters.SetMaxSize( patchesNETGEN_2D_MaxSize )
oppositeSideNETGEN2DParameters.SetMinSize( patchesNETGEN_2D_MinSize )
oppositeSideNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
oppositeSideNETGEN2DParameters.SetOptimize( 1 )
oppositeSideNETGEN2DParameters.SetFineness( 5 )
oppositeSideNETGEN2DParameters.SetGrowthRate( patchesGrowthRate2D )
oppositeSideNETGEN2DParameters.SetNbSegPerEdge( 2 )
oppositeSideNETGEN2DParameters.SetNbSegPerRadius( 3 )
oppositeSideNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
oppositeSideNETGEN2DParameters.SetSecondOrder( 0 )
oppositeSideNETGEN2DParameters.SetFuseEdges( 1 )

#Sub-mesh for topSide boundary
topSideRegular1D = porousCellMesh.Segment(geom=topSide)

#topSideDeflection = topSideRegular1D.Deflection1D(1e-5)
topSideLocalLength = topSideRegular1D.LocalLength(patchesNETGEN_2D_MinSize, None, 1e-6)

topSideNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=topSide)

topSideNETGEN2DParameters = topSideNETGEN2D.Parameters()
topSideNETGEN2DParameters.SetMaxSize( patchesNETGEN_2D_MaxSize )
topSideNETGEN2DParameters.SetMinSize( patchesNETGEN_2D_MinSize )
topSideNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
topSideNETGEN2DParameters.SetOptimize( 1 )
topSideNETGEN2DParameters.SetFineness( 5 )
topSideNETGEN2DParameters.SetGrowthRate( patchesGrowthRate2D )
topSideNETGEN2DParameters.SetNbSegPerEdge( 2 )
topSideNETGEN2DParameters.SetNbSegPerRadius( 3 )
topSideNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
topSideNETGEN2DParameters.SetSecondOrder( 0 )
topSideNETGEN2DParameters.SetFuseEdges( 1 )

#Projection leftSide's sub-mesh on rightSide boundary
rightSideProjection1D2D = porousCellMesh.Projection1D2D(geom=rightSide)
leftSideSourceFace = rightSideProjection1D2D.SourceFace(leftSide,None,None,None,None,None)

#Projection oppositeSide's sub-mesh on frontSide boundary
frontSideProjection1D2D = porousCellMesh.Projection1D2D(geom=frontSide)
oppositeSideSourceFace = frontSideProjection1D2D.SourceFace(oppositeSide,None,None,None,None,None)

#Projection topSide's sub-mesh on bottomSide boundary
bottomSideProjection1D2D = porousCellMesh.Projection1D2D(geom=bottomSide)
topSideSourceFace = bottomSideProjection1D2D.SourceFace(topSide,None,None,None,None,None)

throatsNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=throats)

throatsNETGEN1D2DParameters = throatsNETGEN1D2D.Parameters()
throatsNETGEN1D2DParameters.SetMaxSize( throatsNETGEN_2D_MaxSize )
throatsNETGEN1D2DParameters.SetMinSize( throatsNETGEN_2D_MinSize )
throatsNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
throatsNETGEN1D2DParameters.SetOptimize( 1 )
throatsNETGEN1D2DParameters.SetFineness( 5 )
throatsNETGEN1D2DParameters.SetGrowthRate( throatsGrowthRate2D )
throatsNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
throatsNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
throatsNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
throatsNETGEN1D2DParameters.SetSecondOrder( 0 )
throatsNETGEN1D2DParameters.SetFuseEdges( 1 )

gloabalNETGEN3D = porousCellMesh.Tetrahedron()

globalNETGEN3DParameters = gloabalNETGEN3D.Parameters()
globalNETGEN3DParameters.SetMaxSize( globalNETGEN_3D_MaxSize )
globalNETGEN3DParameters.SetMinSize( globalNETGEN_3D_MinSize )
globalNETGEN3DParameters.SetOptimize( 1 )
globalNETGEN3DParameters.SetFineness( 5 )
globalNETGEN3DParameters.SetGrowthRate( globalGrowthRate3D )
globalNETGEN3DParameters.SetUseSurfaceCurvature( 0 )
globalNETGEN3DParameters.SetSecondOrder( 0 )
globalNETGEN3DParameters.SetFuseEdges( 1 )

throatsNETGEN3D = porousCellMesh.Tetrahedron(algo=smeshBuilder.NETGEN_3D,geom=throats)

throatsNETGEN3DParameters = throatsNETGEN3D.Parameters()
throatsNETGEN3DParameters.SetMaxSize( throatsNETGEN_3D_MaxSize )
throatsNETGEN3DParameters.SetMinSize( throatsNETGEN_3D_MinSize )
throatsNETGEN3DParameters.SetOptimize( 1 )
throatsNETGEN3DParameters.SetFineness( 5 )
throatsNETGEN3DParameters.SetGrowthRate( throatsGrowthRate3D )
throatsNETGEN3DParameters.SetUseSurfaceCurvature( 0 )
throatsNETGEN3DParameters.SetSecondOrder( 0 )
throatsNETGEN3DParameters.SetFuseEdges( 1 )

if (
    intersectionParameter != 0
):
    #viscousLayersSize2D = 0.025 / (VTranslationDirectionNbTimes * WTranslationDirectionNbTimes) * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
    #viscousLayersSize3D = 0.01 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)
    #viscousLayersSize3D = 0.03 / (VTranslationDirectionNbTimes * WTranslationDirectionNbTimes) * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
    #viscousLayersSize3D = 0.5 * localLengthSize
    #Viscous_Layers_2D_1 = NETGEN_1D_2D_1.ViscousLayers2D(viscousLayersSize2D,viscousLayersNumber,viscousLayersGrowth)
    #Viscous_Layers_2D_2 = NETGEN_1D_2D_2.ViscousLayers2D(viscousLayersSize2D,viscousLayersNumber,viscousLayersGrowth)
    #Viscous_Layers_2D_3 = NETGEN_1D_2D_3.ViscousLayers2D(viscousLayersSize2D,viscousLayersNumber,viscousLayersGrowth)
    globalViscousLayers3D = gloabalNETGEN3D.ViscousLayers(
        globalViscousLayersSize3D,
        globalViscousLayersNumber,
        globalViscousLayersGrowth,
        allFacesIDsWithoutViscousLayers,
        1,
        StdMeshersBuilder.FACE_OFFSET
    )
    #throatsViscousLayers3D = throatsNETGEN3D.ViscousLayers(
        #throatsViscousLayersSize3D,
        #throatsViscousLayersNumber,
        #throatsViscousLayersGrowth,
        #allFacesIDsWithoutViscousLayers,
        #1,
        #StdMeshersBuilder.FACE_OFFSET
    #)
    #Viscous_Layers_3D = NETGEN_3D.ViscousLayers(viscousLayersSize3D,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.SURF_OFFSET_SMOOTH)
    #Viscous_Layers_3D = NETGEN_3D.ViscousLayers(viscousLayersSize3D,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.NODE_OFFSET)

#Compute the mesh
isDone = porousCellMesh.Compute()

skeletonWallBoundary = porousCellMesh.GroupOnGeom(skeletonWall,'skeletonWall',SMESH.FACE)
leftSideBoundary = porousCellMesh.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
rightSideBoundary = porousCellMesh.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
oppositeSideBoundary = porousCellMesh.GroupOnGeom(oppositeSide,'oppositeSide',SMESH.FACE)
frontSideBoundary = porousCellMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
topSideBoundary = porousCellMesh.GroupOnGeom(topSide,'topSide',SMESH.FACE)
bottomSideBoundary = porousCellMesh.GroupOnGeom(bottomSide,'bottomSide',SMESH.FACE)

## Get viscous layers
#aCriteria = []
#aCriterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,skeletonWall,SMESH.FT_Undefined,SMESH.FT_Undefined,0.9 * viscousLayersSize3D)
#aCriteria.append(aCriterion)
#aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
#aFilter_1.SetMesh(porousCellMesh.GetMesh())
#viscousLayers = porousCellMesh.GroupOnFilter( SMESH.VOLUME, 'viscousLayers', aFilter_1 )

## Get film
#aCriteria = []
#aCriterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,skeletonWall,SMESH.FT_Undefined,SMESH.FT_Undefined,1e-7)
#aCriteria.append(aCriterion)
#aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
#aFilter_1.SetMesh(porousCellMesh.GetMesh())
#film = porousCellMesh.GroupOnFilter( SMESH.VOLUME, 'film', aFilter_1 )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,skeletonWall,SMESH.FT_Undefined,SMESH.FT_Undefined,1.1 * globalViscousLayersSize3D)
aCriteria.append(aCriterion)
aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_1.SetMesh(porousCellMesh.GetMesh())
film = porousCellMesh.GroupOnFilter( SMESH.VOLUME, 'film', aFilter_1 )

porousCellMesh.ExportUNV( ABSOLUTE_CASE_PATH + "/porousCell.unv" )

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser(1)

##Save geometry and mesh to file
print("Save geometry and mesh to " + ABSOLUTE_CASE_PATH + "/porousCellGeometryAndMesh.hdf\n")

salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/porousCellGeometryAndMesh.hdf", salome.myStudy, False)

print("done")
print("***********************************\n")

