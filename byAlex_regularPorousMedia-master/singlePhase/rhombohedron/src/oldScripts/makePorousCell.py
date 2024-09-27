# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

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

UTranslationDirectionNbTimes = int(sys.argv[4])
VTranslationDirectionNbTimes = int(sys.argv[5])
WTranslationDirectionNbTimes = int(sys.argv[6])

onePoreCellSize = 1.0
grainSize0 = 0.5
filletFactor = 0
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
grainCenterOnBottom_1 = geompy.MakeVertex(
    0 * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(0)
)
grainCenterOnBottom_2 = geompy.MakeVertex(
    onePoreCellSize * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(0)
)
grainCenterOnBottom_3 = geompy.MakeVertex(
    onePoreCellSize + onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(0)
)
grainCenterOnBottom_4 = geompy.MakeVertex(
    onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(0)
)

grainCenterOnTop_1 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + 0 * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(alpha)
)
grainCenterOnTop_2 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + onePoreCellSize * math.cos(0) * math.cos(0),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(0),
    onePoreCellSize * math.sin(alpha)
)
grainCenterOnTop_3 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + onePoreCellSize + onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(alpha)
)
grainCenterOnTop_4 = geompy.MakeVertex(
    onePoreCellSize * math.cos(alpha) * math.cos(0.5 * theta) + onePoreCellSize * math.cos(0) * math.cos(theta),
    onePoreCellSize * math.cos(alpha) * math.sin(0.5 * theta) + onePoreCellSize * math.cos(0) * math.sin(theta),
    onePoreCellSize * math.sin(alpha)
)

#Make the cell's vertex
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

UTranslationDirectionVector = geompy.MakeVector(grainCenterOnBottom_1, grainCenterOnBottom_2)
VTranslationDirectionVector = geompy.MakeVector(grainCenterOnBottom_1, grainCenterOnBottom_4)
WTranslationDirectionVector = geompy.MakeVector(grainCenterOnBottom_1, grainCenterOnTop_1)

#building the porous volume
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

#Make the grains with its rotations
#Grain #1
grain_1 = geompy.MakeSpherePntR(grainCenterOnBottom_1, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + math.pi), onePoreCellSize * math.sin(0.5 * theta + math.pi), 0)
grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * theta + math.pi)
grain_1 = geompy.MakeRotation(grain_1, rotationVector_2, 0.5 * math.pi)

grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * math.pi)
mirroredGrain_1 = geompy.MakeRotation(grain_1, rotationVector_1, math.pi)

#Grain #2
grain_2 = geompy.MakeSpherePntR(grainCenterOnBottom_2, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 1.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_2)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_2)
grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
grain_2 = geompy.MakeRotation(grain_2, rotationVector_2, 0.5 * math.pi)

grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * math.pi)
mirroredGrain_2 = geompy.MakeRotation(grain_2, rotationVector_1, math.pi)

#Grain #3
grain_3 = geompy.MakeSpherePntR(grainCenterOnBottom_3, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta), onePoreCellSize * math.sin(0.5 * theta), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_3)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_3)
grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * theta)
grain_3 = geompy.MakeRotation(grain_3, rotationVector_2, 0.5 * math.pi)

grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * math.pi)
mirroredGrain_3 = geompy.MakeRotation(grain_3, rotationVector_1, math.pi)

#Grain #4
grain_4 = geompy.MakeSpherePntR(grainCenterOnBottom_4, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 0.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_4)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_4)
grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
grain_4 = geompy.MakeRotation(grain_4, rotationVector_2, 0.5 * math.pi)

grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * math.pi)
mirroredGrain_4 = geompy.MakeRotation(grain_4, rotationVector_1, math.pi)

#Grain #5
grain_5 = geompy.MakeSpherePntR(grainCenterOnTop_1, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + math.pi), onePoreCellSize * math.sin(0.5 * theta + math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_1)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_1)
grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * theta + math.pi)
grain_5 = geompy.MakeRotation(grain_5, rotationVector_2, 0.5 * math.pi)

grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * math.pi)
mirroredGrain_5 = geompy.MakeRotation(grain_5, rotationVector_1, math.pi)

#Grain #6
grain_6 = geompy.MakeSpherePntR(grainCenterOnTop_2, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 1.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_2)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_2)
grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
grain_6 = geompy.MakeRotation(grain_6, rotationVector_2, 0.5 * math.pi)

grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * math.pi)
mirroredGrain_6 = geompy.MakeRotation(grain_6, rotationVector_1, math.pi)

#Grain #7
grain_7 = geompy.MakeSpherePntR(grainCenterOnTop_3, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta), onePoreCellSize * math.sin(0.5 * theta), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_3)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_3)
grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * theta)
grain_7 = geompy.MakeRotation(grain_7, rotationVector_2, 0.5 * math.pi)

grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * math.pi)
mirroredGrain_7 = geompy.MakeRotation(grain_7, rotationVector_1, math.pi)

#Grain #8
grain_8 = geompy.MakeSpherePntR(grainCenterOnTop_4, grainSize)
rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, onePoreCellSize)
rotationVector_2 = geompy.MakeVectorDXDYDZ(onePoreCellSize * math.cos(0.5 * theta + 0.5 * math.pi), onePoreCellSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_4)
geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_4)
grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
grain_8 = geompy.MakeRotation(grain_8, rotationVector_2, 0.5 * math.pi)

grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * math.pi)
mirroredGrain_8 = geompy.MakeRotation(grain_8, rotationVector_1, math.pi)

#building the porous volume
#Make faces of porous volume
onePoreBoxLeftFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_4, grainCenterOnTop_4, grainCenterOnTop_1)
onePoreBoxRightFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_2, grainCenterOnBottom_3, grainCenterOnTop_3, grainCenterOnTop_2)
onePoreBoxOppositeFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_4, grainCenterOnBottom_3, grainCenterOnTop_3, grainCenterOnTop_4)
onePoreBoxFrontFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_2, grainCenterOnTop_2, grainCenterOnTop_1)
onePoreBoxTopFace = geompy.MakeQuad4Vertices(grainCenterOnTop_1, grainCenterOnTop_2, grainCenterOnTop_3, grainCenterOnTop_4)
onePoreBoxBottomFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_2, grainCenterOnBottom_3, grainCenterOnBottom_4)

onePoreBox = geompy.MakeHexa(onePoreBoxLeftFace, onePoreBoxRightFace, onePoreBoxOppositeFace, onePoreBoxFrontFace, onePoreBoxTopFace, onePoreBoxBottomFace)

fusedGrains = geompy.MakeFuseList([grain_1, grain_2, grain_3, grain_4, grain_5, grain_6, grain_7, grain_8, mirroredGrain_1, mirroredGrain_2, mirroredGrain_3, mirroredGrain_4, mirroredGrain_5, mirroredGrain_6, mirroredGrain_7, mirroredGrain_8])

onePore = geompy.MakeCut(onePoreBox, fusedGrains)

onePore = geompy.RemoveExtraEdges(onePore, False)

#Get edges IDs with fillet
if (
    filletFactor != 0 and
    intersectionParameter != 0
):
    allEdgesList = []
    allEdgesWithFillet = []
    allEdgesIDsWithFillet = []
    allEdgesWithoutFillet = []

    allEdgesList = geompy.SubShapeAll(onePore, geompy.ShapeType["EDGE"])
    allEdgesWithFillet = geompy.CreateGroup(onePore, geompy.ShapeType["EDGE"])

    geompy.UnionList(allEdgesWithFillet, allEdgesList)

    allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], normalVectorOnLeftFace, middlePointOnLeftFace, GEOM.ST_ON))
    allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], normalVectorOnRightFace, middlePointOnRightFace, GEOM.ST_ON))
    allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], normalVectorOnOppositeFace, middlePointOnOppositeFace, GEOM.ST_ON))
    allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], normalVectorOnFrontFace, middlePointOnFrontFace, GEOM.ST_ON))
    allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], normalVectorOnTopFace, middlePointOnTopFace, GEOM.ST_ON))
    allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], normalVectorOnBottomFace, middlePointOnBottomFace, GEOM.ST_ON))

    geompy.DifferenceList(allEdgesWithFillet, allEdgesWithoutFillet)
    allEdgesIDsWithFillet = geompy.GetObjectIDs(allEdgesWithFillet)

    #Make fillet
    print("filletRadius=" + str(filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])) + "\n")

    fusedFilletedGrains = geompy.MakeFilletAll(fusedGrains, filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

    onePore = geompy.MakeCut(onePoreBox, fusedFilletedGrains)

    onePore = geompy.RemoveExtraEdges(onePore, False)

porousCell = geompy.MakeMultiTranslation2D(onePore, UTranslationDirectionVector, onePoreCellSize, UTranslationDirectionNbTimes, VTranslationDirectionVector, onePoreCellSize, VTranslationDirectionNbTimes)

porousCell = geompy.MakeMultiTranslation1D(porousCell, WTranslationDirectionVector, onePoreCellSize, WTranslationDirectionNbTimes)

listOfSolids = []
listOfSolids = geompy.SubShapeAllSortedCentres(porousCell, geompy.ShapeType["SOLID"])

porousCell = geompy.MakeFuseList(listOfSolids)

facesListOnLeftSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnLeftFace, porousCellBoxMiddlePointOnLeftFace, GEOM.ST_ON)
facesListOnRightSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnRightFace, porousCellBoxMiddlePointOnRightFace, GEOM.ST_ON)
facesListOnOppositeSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnOppositeFace, porousCellBoxMiddlePointOnOppositeFace, GEOM.ST_ON)
facesListOnFrontSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnFrontFace, porousCellBoxMiddlePointOnFrontFace, GEOM.ST_ON)
facesListOnTopSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnTopFace, porousCellBoxMiddlePointOnTopFace, GEOM.ST_ON)
facesListOnBottomSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnBottomFace, porousCellBoxMiddlePointOnBottomFace, GEOM.ST_ON)

## Get fases IDs without viscous layer
allFacesIDsWithoutViscousLayers = []

allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnLeftFace, porousCellBoxMiddlePointOnLeftFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnRightFace, porousCellBoxMiddlePointOnRightFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnOppositeFace, porousCellBoxMiddlePointOnOppositeFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnFrontFace, porousCellBoxMiddlePointOnFrontFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnTopFace, porousCellBoxMiddlePointOnTopFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnBottomFace, porousCellBoxMiddlePointOnBottomFace, GEOM.ST_ON))

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

#geompy.ExportSTL(onePore, ABSOLUTE_CASE_PATH + "/stl/porousCell.stl", True, 1e-4, True)

#geompy.ExportSTL(leftSide, ABSOLUTE_CASE_PATH + "/stl/leftSide.stl", True, 1e-4, True)
#geompy.ExportSTL(rightSide, ABSOLUTE_CASE_PATH + "/stl/rightSide.stl", True, 1e-4, True)
#geompy.ExportSTL(oppositeSide, ABSOLUTE_CASE_PATH + "/stl/oppositeSide.stl", True, 1e-4, True)
#geompy.ExportSTL(frontSide, ABSOLUTE_CASE_PATH + "/stl/frontSide.stl", True, 1e-4, True)
#geompy.ExportSTL(topSide, ABSOLUTE_CASE_PATH + "/stl/topSide.stl", True, 1e-4, True)
#geompy.ExportSTL(bottomSide, ABSOLUTE_CASE_PATH + "/stl/bottomSide.stl", True, 1e-4, True)

###
### SMESH component
###

####################################################
##         Begin of mesh variables section        ##
####################################################
NETGEN_3D_MaxSize = 0.2 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)
NETGEN_3D_MinSize = 0.002 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)

NETGEN_2D_MaxSize = 0.2 * geompy.BasicProperties(onePore)[1] ** (1.0 / 2.0)
NETGEN_2D_MinSize = 0.002 * geompy.BasicProperties(onePore)[1] ** (1.0 / 2.0)
NETGEN_2D_MaxSize_1 = 0.2 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
NETGEN_2D_MinSize_1 = 0.002 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
NETGEN_2D_QuadAllowed = 0

viscousLayersNumber = 2
viscousLayersGrowth = 1.5

Local_Length_Size = 0.1 / (VTranslationDirectionNbTimes * WTranslationDirectionNbTimes) * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
####################################################
##          End of mesh variables section         ##
####################################################

smesh = smeshBuilder.New(theStudy)
porousCellMesh = smesh.Mesh(porousCell)

Regular_1D = porousCellMesh.Segment()
NETGEN_2D_ONLY = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D)
NETGEN_3D = porousCellMesh.Tetrahedron()

#Global parameters for 1D, 2D and 3D mesh
Local_Length = Regular_1D.LocalLength( Local_Length_Size, None, 1e-6 )
#Deflection = Regular_1D.Deflection1D(4e-5)

NETGEN_Parameters_2D_ONLY = NETGEN_2D_ONLY.Parameters()
NETGEN_Parameters_2D_ONLY.SetMaxSize( NETGEN_2D_MaxSize )
NETGEN_Parameters_2D_ONLY.SetMinSize( NETGEN_2D_MinSize )
NETGEN_Parameters_2D_ONLY.SetQuadAllowed( NETGEN_2D_QuadAllowed )
NETGEN_Parameters_2D_ONLY.SetOptimize( 1 )
NETGEN_Parameters_2D_ONLY.SetFineness( 5 )
NETGEN_Parameters_2D_ONLY.SetGrowthRate( 0.1 )
NETGEN_Parameters_2D_ONLY.SetUseSurfaceCurvature( 1 )
NETGEN_Parameters_2D_ONLY.SetSecondOrder( 0 )
NETGEN_Parameters_2D_ONLY.SetFuseEdges( 1 )

NETGEN_Parameters_3D = NETGEN_3D.Parameters()
NETGEN_Parameters_3D.SetMaxSize( NETGEN_3D_MaxSize )
NETGEN_Parameters_3D.SetMinSize( NETGEN_3D_MinSize )
NETGEN_Parameters_3D.SetOptimize( 1 )
NETGEN_Parameters_3D.SetFineness( 5 )
NETGEN_Parameters_3D.SetGrowthRate( 0.1 )
NETGEN_Parameters_3D.SetUseSurfaceCurvature( 1 )
NETGEN_Parameters_3D.SetSecondOrder( 0 )
NETGEN_Parameters_3D.SetFuseEdges( 1 )

if (
    intersectionParameter != 0
):
    viscousLayersSize = 0.025 / (VTranslationDirectionNbTimes * WTranslationDirectionNbTimes) * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
    Viscous_Layers = NETGEN_3D.ViscousLayers(viscousLayersSize,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.FACE_OFFSET)
    #Viscous_Layers = NETGEN_3D.ViscousLayers(viscousLayersSize,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.SURF_OFFSET_SMOOTH)
    #Viscous_Layers = NETGEN_3D.ViscousLayers(viscousLayersSize,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.NODE_OFFSET)

#Sub-mesh for leftSide boundary
Regular_1D_1 = porousCellMesh.Segment(geom=leftSide)
NETGEN_2D_ONLY_1 = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=leftSide)

#Deflection_1 = Regular_1D_1.Deflection1D(1e-5)
Local_Length_1 = Regular_1D_1.LocalLength(Local_Length_Size, None, 1e-6)

NETGEN_Parameters_2D_ONLY_1 = NETGEN_2D_ONLY_1.Parameters()
NETGEN_Parameters_2D_ONLY_1.SetMaxSize( NETGEN_2D_MaxSize_1 )
NETGEN_Parameters_2D_ONLY_1.SetMinSize( NETGEN_2D_MinSize_1 )
NETGEN_Parameters_2D_ONLY_1.SetQuadAllowed( NETGEN_2D_QuadAllowed )
NETGEN_Parameters_2D_ONLY_1.SetOptimize( 1 )
NETGEN_Parameters_2D_ONLY_1.SetFineness( 5 )
NETGEN_Parameters_2D_ONLY_1.SetGrowthRate( 0.1 )
NETGEN_Parameters_2D_ONLY_1.SetUseSurfaceCurvature( 1 )
NETGEN_Parameters_2D_ONLY_1.SetSecondOrder( 0 )
NETGEN_Parameters_2D_ONLY_1.SetFuseEdges( 1 )

#Projection leftSide's sub-mesh on rightSide boundary
Projection_1D2D_1 = porousCellMesh.Projection1D2D(geom=rightSide)
Source_Face_1 = Projection_1D2D_1.SourceFace(leftSide,None,None,None,None,None)

#Sub-mesh for oppositeSide boundary
Regular_1D_2 = porousCellMesh.Segment(geom=oppositeSide)
NETGEN_2D_ONLY_2 = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=oppositeSide)

#Deflection_2 = Regular_1D_2.Deflection1D(1e-5)
Local_Length_2 = Regular_1D_2.LocalLength(Local_Length_Size, None, 1e-6)

NETGEN_Parameters_2D_ONLY_2 = NETGEN_2D_ONLY_2.Parameters()
NETGEN_Parameters_2D_ONLY_2.SetMaxSize( NETGEN_2D_MaxSize_1 )
NETGEN_Parameters_2D_ONLY_2.SetMinSize( NETGEN_2D_MinSize_1 )
NETGEN_Parameters_2D_ONLY_2.SetQuadAllowed( NETGEN_2D_QuadAllowed )
NETGEN_Parameters_2D_ONLY_2.SetOptimize( 1 )
NETGEN_Parameters_2D_ONLY_2.SetFineness( 5 )
NETGEN_Parameters_2D_ONLY_2.SetGrowthRate( 0.1 )
NETGEN_Parameters_2D_ONLY_2.SetUseSurfaceCurvature( 1 )
NETGEN_Parameters_2D_ONLY_2.SetSecondOrder( 0 )
NETGEN_Parameters_2D_ONLY_2.SetFuseEdges( 1 )

#status = porousCellMesh.AddHypothesis(NETGEN_Parameters_2D_ONLY_1,oppositeSide)

#Projection oppositeSide's sub-mesh on frontSide boundary
Projection_1D2D_2 = porousCellMesh.Projection1D2D(geom=frontSide)
Source_Face_2 = Projection_1D2D_2.SourceFace(oppositeSide,None,None,None,None,None)

#Sub-mesh for topSide boundary
Regular_1D_3 = porousCellMesh.Segment(geom=topSide)
NETGEN_2D_ONLY_3 = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=topSide)

#Deflection_3 = Regular_1D_3.Deflection1D(1e-5)
Local_Length_3 = Regular_1D_3.LocalLength(Local_Length_Size, None, 1e-6)

NETGEN_Parameters_2D_ONLY_3 = NETGEN_2D_ONLY_3.Parameters()
NETGEN_Parameters_2D_ONLY_3.SetMaxSize( NETGEN_2D_MaxSize_1 )
NETGEN_Parameters_2D_ONLY_3.SetMinSize( NETGEN_2D_MinSize_1 )
NETGEN_Parameters_2D_ONLY_3.SetQuadAllowed( NETGEN_2D_QuadAllowed )
NETGEN_Parameters_2D_ONLY_3.SetOptimize( 1 )
NETGEN_Parameters_2D_ONLY_3.SetFineness( 5 )
NETGEN_Parameters_2D_ONLY_3.SetGrowthRate( 0.1 )
NETGEN_Parameters_2D_ONLY_3.SetUseSurfaceCurvature( 1 )
NETGEN_Parameters_2D_ONLY_3.SetSecondOrder( 0 )
NETGEN_Parameters_2D_ONLY_3.SetFuseEdges( 1 )

#status = porousCellMesh.AddHypothesis(NETGEN_Parameters_2D_ONLY_1,topSide)

#Projection topSide's sub-mesh on bottomSide boundary
Projection_1D2D_3 = porousCellMesh.Projection1D2D(geom=bottomSide)
Source_Face_3 = Projection_1D2D_3.SourceFace(topSide,None,None,None,None,None)

#Compute the mesh
isDone = porousCellMesh.Compute()

leftSideBoundary = porousCellMesh.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
rightSideBoundary = porousCellMesh.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
oppositeSideBoundary = porousCellMesh.GroupOnGeom(oppositeSide,'oppositeSide',SMESH.FACE)
frontSideBoundary = porousCellMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
topSideBoundary = porousCellMesh.GroupOnGeom(topSide,'topSide',SMESH.FACE)
bottomSideBoundary = porousCellMesh.GroupOnGeom(bottomSide,'bottomSide',SMESH.FACE)

porousCellMesh.ExportUNV( ABSOLUTE_CASE_PATH + "/porousCell.unv" )

print("done")
print("***********************************\n")

#if salome.sg.hasDesktop():
#    salome.sg.updateObjBrowser(True)
