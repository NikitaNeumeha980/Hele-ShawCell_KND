# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import sys
import math
import salome

salome.salome_init()
theStudy = salome.myStudy

####################################################
##            Begin of variables section          ##
####################################################
#For debuging
#CASE_DIR = "/home/alexshtil/PETER"

#CASE_DIR = sys.argv[0]
#print("CASE DIRECTORY: " + CASE_DIR)

## Log file
#outputFile = open(CASE_DIR + "/log.outputFile", "w", 0)

minTheta = 80
maxTheta = 90
thetaStep = 5

minIntersectionParameter = 0
maxIntersectionParameter = 0.1
intersectionParameterStep = 0.05

UTranslationDirectionNbTimes = 2
VTranslationDirectionNbTimes = 2
WTranslationDirectionNbTimes = 2
poreSize = 1.0
grainSize0 = 0.5
####################################################
##             End of variables section           ##
####################################################

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import SALOMEDS

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder
from salome.StdMeshers import StdMeshersBuilder

geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

for theta in [(minTheta + x * thetaStep) for x in range(0, int((maxTheta - minTheta) / thetaStep) + 1)]:

    theta = math.radians(theta)
    alpha = math.acos(math.cos(theta) / math.cos(0.5 * theta))

    #For debuging
    print('theta=' + str(math.degrees(theta)))
    print('alpha=' + str(math.degrees(math.acos(math.cos(theta) / math.cos(0.5 * theta)))))

    for intersectionParameter in [(minIntersectionParameter + x * intersectionParameterStep) for x in range(0, int((maxIntersectionParameter - minIntersectionParameter) / intersectionParameterStep) + 1)]:

        #For debuging
        print('intersectionParameter=' + str(intersectionParameter))

        grainSize = grainSize0 / (1.0 - intersectionParameter)

        #Building the sceleton
        #Make the grain's centeres
        grainCenterOnBottom_1 = geompy.MakeVertex(
                0 * math.cos(0) * math.cos(0),
                poreSize * math.cos(0) * math.sin(0),
                poreSize * math.sin(0)
            )
        grainCenterOnBottom_2 = geompy.MakeVertex(
                poreSize * math.cos(0) * math.cos(0),
                poreSize * math.cos(0) * math.sin(0),
                poreSize * math.sin(0)
            )
        grainCenterOnBottom_3 = geompy.MakeVertex(
                poreSize + poreSize * math.cos(0) * math.cos(theta),
                poreSize * math.cos(0) * math.sin(theta),
                poreSize * math.sin(0)
            )
        grainCenterOnBottom_4 = geompy.MakeVertex(
                poreSize * math.cos(0) * math.cos(theta),
                poreSize * math.cos(0) * math.sin(theta),
                poreSize * math.sin(0)
            )

        grainCenterOnTop_1 = geompy.MakeVertex(
                poreSize * math.cos(alpha) * math.cos(0.5 * theta) + 0 * math.cos(0) * math.cos(0),
                poreSize * math.cos(alpha) * math.sin(0.5 * theta) + poreSize * math.cos(0) * math.sin(0),
                poreSize * math.sin(alpha)
            )
        grainCenterOnTop_2 = geompy.MakeVertex(
                poreSize * math.cos(alpha) * math.cos(0.5 * theta) + poreSize * math.cos(0) * math.cos(0),
                poreSize * math.cos(alpha) * math.sin(0.5 * theta) + poreSize * math.cos(0) * math.sin(0),
                poreSize * math.sin(alpha)
            )
        grainCenterOnTop_3 = geompy.MakeVertex(
                poreSize * math.cos(alpha) * math.cos(0.5 * theta) + poreSize + poreSize * math.cos(0) * math.cos(theta),
                poreSize * math.cos(alpha) * math.sin(0.5 * theta) + poreSize * math.cos(0) * math.sin(theta),
                poreSize * math.sin(alpha)
            )
        grainCenterOnTop_4 = geompy.MakeVertex(
                poreSize * math.cos(alpha) * math.cos(0.5 * theta) + poreSize * math.cos(0) * math.cos(theta),
                poreSize * math.cos(alpha) * math.sin(0.5 * theta) + poreSize * math.cos(0) * math.sin(theta),
                poreSize * math.sin(alpha)
            )

        #Make the cell's vertex
        porousCellBoxVertexOnBottom_1 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * 0 * math.cos(0) * math.cos(0),
                VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(0),
                WTranslationDirectionNbTimes * poreSize * math.sin(0)
            )
        porousCellBoxVertexOnBottom_2 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize * math.cos(0) * math.cos(0),
                VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(0),
                WTranslationDirectionNbTimes * poreSize * math.sin(0)
            )
        porousCellBoxVertexOnBottom_3 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize + UTranslationDirectionNbTimes * poreSize * math.cos(0) * math.cos(theta),
                VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(theta),
                WTranslationDirectionNbTimes * poreSize * math.sin(0)
            )
        porousCellBoxVertexOnBottom_4 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize * math.cos(0) * math.cos(theta),
                VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(theta),
                WTranslationDirectionNbTimes * poreSize * math.sin(0)
            )

        porousCellBoxVertexOnTop_1 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * 0 * math.cos(0) * math.cos(0),
                VTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(0),
                WTranslationDirectionNbTimes * poreSize * math.sin(alpha)
            )
        porousCellBoxVertexOnTop_2 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * poreSize * math.cos(0) * math.cos(0),
                VTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(0),
                WTranslationDirectionNbTimes * poreSize * math.sin(alpha)
            )
        porousCellBoxVertexOnTop_3 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * poreSize + UTranslationDirectionNbTimes * poreSize * math.cos(0) * math.cos(theta),
                VTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(theta),
                WTranslationDirectionNbTimes * poreSize * math.sin(alpha)
            )
        porousCellBoxVertexOnTop_4 = geompy.MakeVertex(
                UTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.cos(0.5 * theta) + UTranslationDirectionNbTimes * poreSize * math.cos(0) * math.cos(theta),
                VTranslationDirectionNbTimes * poreSize * math.cos(alpha) * math.sin(0.5 * theta) + VTranslationDirectionNbTimes * poreSize * math.cos(0) * math.sin(theta),
                WTranslationDirectionNbTimes * poreSize * math.sin(alpha)
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
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta + math.pi), poreSize * math.sin(0.5 * theta + math.pi), 0)
        grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * theta + math.pi)
        grain_1 = geompy.MakeRotation(grain_1, rotationVector_2, 0.5 * math.pi)

        #Grain #2
        grain_2 = geompy.MakeSpherePntR(grainCenterOnBottom_2, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta + 1.5 * math.pi), poreSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_2)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_2)
        grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
        grain_2 = geompy.MakeRotation(grain_2, rotationVector_2, 0.5 * math.pi)

        #Grain #3
        grain_3 = geompy.MakeSpherePntR(grainCenterOnBottom_3, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta), poreSize * math.sin(0.5 * theta), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_3)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_3)
        grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * theta)
        grain_3 = geompy.MakeRotation(grain_3, rotationVector_2, 0.5 * math.pi)

        #Grain #4
        grain_4 = geompy.MakeSpherePntR(grainCenterOnBottom_4, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta + 0.5 * math.pi), poreSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_4)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_4)
        grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
        grain_4 = geompy.MakeRotation(grain_4, rotationVector_2, 0.5 * math.pi)

        #Grain #5
        grain_5 = geompy.MakeSpherePntR(grainCenterOnTop_1, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta + math.pi), poreSize * math.sin(0.5 * theta + math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_1)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_1)
        grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * theta + math.pi)
        grain_5 = geompy.MakeRotation(grain_5, rotationVector_2, 0.5 * math.pi)

        #Grain #6
        grain_6 = geompy.MakeSpherePntR(grainCenterOnTop_2, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta + 1.5 * math.pi), poreSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_2)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_2)
        grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
        grain_6 = geompy.MakeRotation(grain_6, rotationVector_2, 0.5 * math.pi)

        #Grain #7
        grain_7 = geompy.MakeSpherePntR(grainCenterOnTop_3, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta), poreSize * math.sin(0.5 * theta), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_3)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_3)
        grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * theta)
        grain_7 = geompy.MakeRotation(grain_7, rotationVector_2, 0.5 * math.pi)

        #Grain #8
        grain_8 = geompy.MakeSpherePntR(grainCenterOnTop_4, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, poreSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(poreSize * math.cos(0.5 * theta + 0.5 * math.pi), poreSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_4)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_4)
        grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
        grain_8 = geompy.MakeRotation(grain_8, rotationVector_2, 0.5 * math.pi)

        #building the porous volume
        #Make faces of porous volume
        onePoreBoxLeftFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_4, grainCenterOnTop_4, grainCenterOnTop_1)
        onePoreBoxRightFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_2, grainCenterOnBottom_3, grainCenterOnTop_3, grainCenterOnTop_2)
        onePoreBoxOppositeFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_4, grainCenterOnBottom_3, grainCenterOnTop_3, grainCenterOnTop_4)
        onePoreBoxFrontFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_2, grainCenterOnTop_2, grainCenterOnTop_1)
        onePoreBoxTopFace = geompy.MakeQuad4Vertices(grainCenterOnTop_1, grainCenterOnTop_2, grainCenterOnTop_3, grainCenterOnTop_4)
        onePoreBoxBottomFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_2, grainCenterOnBottom_3, grainCenterOnBottom_4)

        onePoreBox = geompy.MakeHexa(onePoreBoxLeftFace, onePoreBoxRightFace, onePoreBoxOppositeFace, onePoreBoxFrontFace, onePoreBoxTopFace, onePoreBoxBottomFace)

        #For debuging
        #geompy.addToStudy( onePoreBox, 'onePoreBox')

        onePore = geompy.MakeCutList(onePoreBox, [grain_1, grain_2, grain_3, grain_4, grain_5, grain_6, grain_7, grain_8], True)

        porousCell = geompy.MakeMultiTranslation2D(onePore, UTranslationDirectionVector, poreSize, UTranslationDirectionNbTimes, VTranslationDirectionVector, poreSize, VTranslationDirectionNbTimes)

        porousCell = geompy.MakeMultiTranslation1D(porousCell, WTranslationDirectionVector, poreSize, WTranslationDirectionNbTimes)

        listOfSolids = []
        listOfSolids = geompy.SubShapeAllSortedCentres(porousCell, geompy.ShapeType["SOLID"])

        porousCell = geompy.MakeFuseList(listOfSolids)

        porousCellSceleton = geompy.MakeCutList(porousCellBox, [porousCell], True)

        geompy.addToStudy( porousCell, 'porousCell_theta' + str(math.degrees(theta)) + 'IP' + str(intersectionParameter))
        geompy.addToStudy( porousCellSceleton, 'porousCellSceleton_theta' + str(math.degrees(theta)) + 'IP' + str(intersectionParameter))

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

        geompy.addToStudyInFather( porousCell, leftSide, 'leftSide' )
        geompy.addToStudyInFather( porousCell, rightSide, 'rightSide' )
        geompy.addToStudyInFather( porousCell, oppositeSide, 'oppositeSide' )
        geompy.addToStudyInFather( porousCell, frontSide, 'frontSide' )
        geompy.addToStudyInFather( porousCell, topSide, 'topSide' )
        geompy.addToStudyInFather( porousCell, bottomSide, 'bottomSide' )

        ###
        ### SMESH component
        ###

        smesh = smeshBuilder.New(theStudy)
        Mesh_1 = smesh.Mesh(porousCell)

        #Global parameters for 2D and 3D mesh
        NETGEN_3D = Mesh_1.Tetrahedron()
        NETGEN_3D_Parameters = NETGEN_3D.Parameters()
        NETGEN_3D_Parameters.SetMaxSize( 0.05 )
        NETGEN_3D_Parameters.SetOptimize( 1 )
        NETGEN_3D_Parameters.SetFineness( 4 )
        NETGEN_3D_Parameters.SetMinSize( 0.001 )
        NETGEN_3D_Parameters.SetUseSurfaceCurvature( 1 )
        NETGEN_3D_Parameters.SetSecondOrder( 0 )
        NETGEN_3D_Parameters.SetFuseEdges( 1 )
        NETGEN_3D_Parameters.SetQuadAllowed( 0 )

        if (
            intersectionParameter != 0
        ):
            Viscous_Layers = NETGEN_3D.ViscousLayers(0.001,5,1.1,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.SURF_OFFSET_SMOOTH)
            Viscous_Layers.SetTotalThickness( 0.01 )
            Viscous_Layers.SetNumberLayers( 5 )
            Viscous_Layers.SetStretchFactor( 1.1 )
            Viscous_Layers.SetMethod( StdMeshersBuilder.SURF_OFFSET_SMOOTH )
            Viscous_Layers.SetFaces( allFacesIDsWithoutViscousLayers, 1 )
            
            ## Set names of Mesh objects
            #smesh.SetName(Viscous_Layers, 'Viscous Layers')

        NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
        NETGEN_2D_Parameters = NETGEN_1D_2D.Parameters()
        NETGEN_2D_Parameters.SetMaxSize( 0.05 )
        NETGEN_2D_Parameters.SetSecondOrder( 0 )
        NETGEN_2D_Parameters.SetOptimize( 1 )
        NETGEN_2D_Parameters.SetFineness( 4 )
        NETGEN_2D_Parameters.SetMinSize( 0.001 )
        NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
        NETGEN_2D_Parameters.SetFuseEdges( 1 )
        NETGEN_2D_Parameters.SetQuadAllowed( 0 )

        isDone = Mesh_1.Compute()

        #Sub-mesh for leftSide boundary
        NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=leftSide)
        NETGEN_2D_Parameters_1 = NETGEN_1D_2D_1.Parameters()
        NETGEN_2D_Parameters_1.SetMaxSize( 0.02 )
        NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
        NETGEN_2D_Parameters_1.SetOptimize( 1 )
        NETGEN_2D_Parameters_1.SetFineness( 4 )
        NETGEN_2D_Parameters_1.SetMinSize( 0.001 )
        NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
        NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
        NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )

        #Projection leftSide's sub-mesh on rightSide boundary
        Projection_1D2D_1 = Mesh_1.Projection1D2D(geom=rightSide)
        Source_Face_1 = Projection_1D2D_1.SourceFace(leftSide,None,None,None,None,None)

        #Sub-mesh for oppositeSide boundary
        NETGEN_1D_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=oppositeSide)
        status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_1,oppositeSide)

        #Projection oppositeSide's sub-mesh on frontSide boundary
        Projection_1D2D_2 = Mesh_1.Projection1D2D(geom=frontSide)
        Source_Face_2 = Projection_1D2D_2.SourceFace(oppositeSide,None,None,None,None,None)

        #Sub-mesh for topSide boundary
        NETGEN_1D_2D_3 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=topSide)
        status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_1,topSide)

        #Projection topSide's sub-mesh on bottomSide boundary
        Projection_1D2D_3 = Mesh_1.Projection1D2D(geom=bottomSide)
        Source_Face_3 = Projection_1D2D_3.SourceFace(topSide,None,None,None,None,None)

        isDone = Mesh_1.Compute()

        leftSideBoundary = Mesh_1.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
        rightSideBoundary = Mesh_1.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
        oppositeSideBoundary = Mesh_1.GroupOnGeom(oppositeSide,'oppositeSide',SMESH.FACE)
        frontSideBoundary = Mesh_1.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
        topSideBoundary = Mesh_1.GroupOnGeom(topSide,'topSide',SMESH.FACE)
        bottomSideBoundary = Mesh_1.GroupOnGeom(bottomSide,'bottomSide',SMESH.FACE)

        #Sub_mesh_1 = NETGEN_1D_2D_1.GetSubMesh()
        #Sub_mesh_2 = Projection_1D2D_1.GetSubMesh()
        #Sub_mesh_3 = NETGEN_1D_2D_2.GetSubMesh()
        #Sub_mesh_4 = Projection_1D2D_2.GetSubMesh()
        #Sub_mesh_5 = NETGEN_1D_2D_3.GetSubMesh()
        #Sub_mesh_6 = Projection_1D2D_3.GetSubMesh()

        ## Set names of Mesh objects
        #smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
        #smesh.SetName(Projection_1D2D_1.GetAlgorithm(), 'Projection_1D2D_1')
        #smesh.SetName(Projection_1D2D_2.GetAlgorithm(), 'Projection_1D2D_2')
        #smesh.SetName(Projection_1D2D_3.GetAlgorithm(), 'Projection_1D2D_3')
        #smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
        #smesh.SetName(NETGEN_1D_2D_1.GetAlgorithm(), 'NETGEN 1D-2D_1')

        #smesh.SetName(NETGEN_3D_Parameters, 'NETGEN 3D Parameters')
        #smesh.SetName(NETGEN_2D_Parameters, 'NETGEN 2D Parameters')
        #smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')

        #smesh.SetName(Source_Face_1, 'Source Face_1')
        #smesh.SetName(Source_Face_2, 'Source Face_2')
        #smesh.SetName(Source_Face_3, 'Source Face_3')
        
        #smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

        #smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
        #smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
        #smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
        #smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
        #smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
        #smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser(True)
