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

minTheta = 60
maxTheta = 60
thetaStep = 5

minIntersectionParameter = 0.01
maxIntersectionParameter = 0.11
intersectionParameterStep = 0.05

cellSize = 1.0
grainSize0 = 0.5
filletFactor = 0
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
                cellSize * math.cos(0) * math.sin(0),
                cellSize * math.sin(0)
            )
        grainCenterOnBottom_2 = geompy.MakeVertex(
                cellSize * math.cos(0) * math.cos(0),
                cellSize * math.cos(0) * math.sin(0),
                cellSize * math.sin(0)
            )
        grainCenterOnBottom_3 = geompy.MakeVertex(
                cellSize + cellSize * math.cos(0) * math.cos(theta),
                cellSize * math.cos(0) * math.sin(theta),
                cellSize * math.sin(0)
            )
        grainCenterOnBottom_4 = geompy.MakeVertex(
                cellSize * math.cos(0) * math.cos(theta),
                cellSize * math.cos(0) * math.sin(theta),
                cellSize * math.sin(0)
            )

        grainCenterOnTop_1 = geompy.MakeVertex(
                cellSize * math.cos(alpha) * math.cos(0.5 * theta) + 0 * math.cos(0) * math.cos(0),
                cellSize * math.cos(alpha) * math.sin(0.5 * theta) + cellSize * math.cos(0) * math.sin(0),
                cellSize * math.sin(alpha)
            )
        grainCenterOnTop_2 = geompy.MakeVertex(
                cellSize * math.cos(alpha) * math.cos(0.5 * theta) + cellSize * math.cos(0) * math.cos(0),
                cellSize * math.cos(alpha) * math.sin(0.5 * theta) + cellSize * math.cos(0) * math.sin(0),
                cellSize * math.sin(alpha)
            )
        grainCenterOnTop_3 = geompy.MakeVertex(
                cellSize * math.cos(alpha) * math.cos(0.5 * theta) + cellSize + cellSize * math.cos(0) * math.cos(theta),
                cellSize * math.cos(alpha) * math.sin(0.5 * theta) + cellSize * math.cos(0) * math.sin(theta),
                cellSize * math.sin(alpha)
            )
        grainCenterOnTop_4 = geompy.MakeVertex(
                cellSize * math.cos(alpha) * math.cos(0.5 * theta) + cellSize * math.cos(0) * math.cos(theta),
                cellSize * math.cos(alpha) * math.sin(0.5 * theta) + cellSize * math.cos(0) * math.sin(theta),
                cellSize * math.sin(alpha)
            )

        #Make the grains with its rotations
        #Grain #1
        grain_1 = geompy.MakeSpherePntR(grainCenterOnBottom_1, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta + math.pi), cellSize * math.sin(0.5 * theta + math.pi), 0)
        grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * theta + math.pi)
        grain_1 = geompy.MakeRotation(grain_1, rotationVector_2, 0.5 * math.pi)

        grain_1 = geompy.MakeRotation(grain_1, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_1 = geompy.MakeRotation(grain_1, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #2
        grain_2 = geompy.MakeSpherePntR(grainCenterOnBottom_2, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta + 1.5 * math.pi), cellSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_2)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_2)
        grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
        grain_2 = geompy.MakeRotation(grain_2, rotationVector_2, 0.5 * math.pi)

        grain_2 = geompy.MakeRotation(grain_2, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_2 = geompy.MakeRotation(grain_2, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #3
        grain_3 = geompy.MakeSpherePntR(grainCenterOnBottom_3, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta), cellSize * math.sin(0.5 * theta), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_3)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_3)
        grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * theta)
        grain_3 = geompy.MakeRotation(grain_3, rotationVector_2, 0.5 * math.pi)

        grain_3 = geompy.MakeRotation(grain_3, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_3 = geompy.MakeRotation(grain_3, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #4
        grain_4 = geompy.MakeSpherePntR(grainCenterOnBottom_4, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta + 0.5 * math.pi), cellSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnBottom_4)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnBottom_4)
        grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
        grain_4 = geompy.MakeRotation(grain_4, rotationVector_2, 0.5 * math.pi)

        grain_4 = geompy.MakeRotation(grain_4, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_4 = geompy.MakeRotation(grain_4, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #5
        grain_5 = geompy.MakeSpherePntR(grainCenterOnTop_1, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta + math.pi), cellSize * math.sin(0.5 * theta + math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_1)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_1)
        grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * theta + math.pi)
        grain_5 = geompy.MakeRotation(grain_5, rotationVector_2, 0.5 * math.pi)

        grain_5 = geompy.MakeRotation(grain_5, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_5 = geompy.MakeRotation(grain_5, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #6
        grain_6 = geompy.MakeSpherePntR(grainCenterOnTop_2, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta + 1.5 * math.pi), cellSize * math.sin(0.5 * theta + 1.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_2)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_2)
        grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * theta + 1.5 * math.pi)
        grain_6 = geompy.MakeRotation(grain_6, rotationVector_2, 0.5 * math.pi)

        grain_6 = geompy.MakeRotation(grain_6, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_6 = geompy.MakeRotation(grain_6, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #7
        grain_7 = geompy.MakeSpherePntR(grainCenterOnTop_3, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta), cellSize * math.sin(0.5 * theta), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_3)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_3)
        grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * theta)
        grain_7 = geompy.MakeRotation(grain_7, rotationVector_2, 0.5 * math.pi)

        grain_7 = geompy.MakeRotation(grain_7, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_7 = geompy.MakeRotation(grain_7, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #Grain #8
        grain_8 = geompy.MakeSpherePntR(grainCenterOnTop_4, grainSize)
        rotationVector_1 = geompy.MakeVectorDXDYDZ(0, 0, cellSize)
        rotationVector_2 = geompy.MakeVectorDXDYDZ(cellSize * math.cos(0.5 * theta + 0.5 * math.pi), cellSize * math.sin(0.5 * theta + 0.5 * math.pi), 0)
        geompy.TranslateTwoPoints(rotationVector_1, grainCenterOnBottom_1, grainCenterOnTop_4)
        geompy.TranslateTwoPoints(rotationVector_2, grainCenterOnBottom_1, grainCenterOnTop_4)
        grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * theta + 0.5 * math.pi)
        grain_8 = geompy.MakeRotation(grain_8, rotationVector_2, 0.5 * math.pi)

        grain_8 = geompy.MakeRotation(grain_8, rotationVector_1, 0.5 * math.pi)
        mirroredGrain_8 = geompy.MakeRotation(grain_8, rotationVector_1, math.pi)

        #For debuging
        #geompy.addToStudy( rotationVector_1, 'rotationVector_1')
        #geompy.addToStudy( rotationVector_2, 'rotationVector_2')

        #For debuging
        #geompy.addToStudy( grain_1, 'grain_1_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_2, 'grain_2_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_3, 'grain_3_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_4, 'grain_4_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_5, 'grain_5_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_6, 'grain_6_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_7, 'grain_7_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))
        #geompy.addToStudy( grain_8, 'grain_8_theta' + str(math.degrees(theta)) + '_intersectionParameter' + str(intersectionParameter))

        #building the porous volume
        #Make faces of porous volume
        leftFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_4, grainCenterOnTop_4, grainCenterOnTop_1)
        rightFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_2, grainCenterOnBottom_3, grainCenterOnTop_3, grainCenterOnTop_2)
        oppositeFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_4, grainCenterOnBottom_3, grainCenterOnTop_3, grainCenterOnTop_4)
        frontFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_2, grainCenterOnTop_2, grainCenterOnTop_1)
        topFace = geompy.MakeQuad4Vertices(grainCenterOnTop_1, grainCenterOnTop_2, grainCenterOnTop_3, grainCenterOnTop_4)
        bottomFace = geompy.MakeQuad4Vertices(grainCenterOnBottom_1, grainCenterOnBottom_2, grainCenterOnBottom_3, grainCenterOnBottom_4)

        #For debuging
        #geompy.addToStudy( leftFace, 'leftFace')
        #geompy.addToStudy( rightFace, 'rightFace')
        #geompy.addToStudy( oppositeFace, 'oppositeFace')
        #geompy.addToStudy( frontFace, 'frontFace')
        #geompy.addToStudy( topFace, 'topFace')
        #geompy.addToStudy( bottomFace, 'bottomFace')

        #Make middle points on boundaries for searching boundaries
        middlePointOnLeftFace = geompy.MakeVertexOnSurface(leftFace, 0.5, 0.5)
        middlePointOnRightFace = geompy.MakeVertexOnSurface(rightFace, 0.5, 0.5)
        middlePointOnOppositeFace = geompy.MakeVertexOnSurface(oppositeFace, 0.5, 0.5)
        middlePointOnFrontFace = geompy.MakeVertexOnSurface(frontFace, 0.5, 0.5)
        middlePointOnTopFace = geompy.MakeVertexOnSurface(topFace, 0.5, 0.5)
        middlePointOnBottomFace = geompy.MakeVertexOnSurface(bottomFace, 0.5, 0.5)

        #For debuging
        #geompy.addToStudy( middlePointOnLeftFace, 'middlePointOnLeftFace')
        #geompy.addToStudy( middlePointOnRightFace, 'middlePointOnRightFace')
        #geompy.addToStudy( middlePointOnOppositeFace, 'middlePointOnOppositeFace')
        #geompy.addToStudy( middlePointOnFrontFace, 'middlePointOnFrontFace')
        #geompy.addToStudy( middlePointOnTopFace, 'middlePointOnTopFace')
        #geompy.addToStudy( middlePointOnBottomFace, 'middlePointOnBottomFace')

        #Make normals for searching boundaries
        normalVectorOnLeftFace = geompy.GetNormal(leftFace)
        normalVectorOnRightFace = geompy.GetNormal(rightFace)
        normalVectorOnOppositeFace = geompy.GetNormal(oppositeFace)
        normalVectorOnFrontFace = geompy.GetNormal(frontFace)
        normalVectorOnTopFace = geompy.GetNormal(topFace)
        normalVectorOnBottomFace = geompy.GetNormal(bottomFace)

        #For debuging
        #geompy.addToStudy( normalVectorOnLeftFace, 'normalVectorOnLeftFace')
        #geompy.addToStudy( normalVectorOnRightFace, 'normalVectorOnRightFace')
        #geompy.addToStudy( normalVectorOnOppositeFace, 'normalVectorOnOppositeFace')
        #geompy.addToStudy( normalVectorOnFrontFace, 'normalVectorOnFrontFace')
        #geompy.addToStudy( normalVectorOnTopFace, 'normalVectorOnTopFace')
        #geompy.addToStudy( normalVectorOnBottomFace, 'normalVectorOnBottomFace')

        onePoreBox = geompy.MakeHexa(leftFace, rightFace, oppositeFace, frontFace, topFace, bottomFace)

        #For debuging
        #geompy.addToStudy( onePoreBox, 'onePoreBox')

        fusedGrains = geompy.MakeFuseList([grain_1, grain_2, grain_3, grain_4, grain_5, grain_6, grain_7, grain_8, mirroredGrain_1, mirroredGrain_2, mirroredGrain_3, mirroredGrain_4, mirroredGrain_5, mirroredGrain_6, mirroredGrain_7, mirroredGrain_8])

        #onePore = geompy.MakeCutList(onePoreBox, [grain_1, grain_2, grain_3, grain_4, grain_5, grain_6, grain_7, grain_8], True)

        #For debuging
        #geompy.addToStudy( fusedGrains, 'fusedGrains_theta' + str(math.degrees(theta)) + 'IP' + str(intersectionParameter))

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
            #onePore = geompy.MakeFillet(onePore, filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]), geompy.ShapeType["EDGE"], allEdgesIDsWithFillet)

            fusedFilletedGrains = geompy.MakeFilletAll(fusedGrains, filletFactor * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

            #For debuging
            #geompy.addToStudy( fusedFilletedGrains, 'fusedFilletedGrains_theta' + str(math.degrees(theta)) + 'IP' + str(intersectionParameter))

            onePore = geompy.MakeCut(onePoreBox, fusedFilletedGrains)

            onePore = geompy.RemoveExtraEdges(onePore, False)

        geompy.addToStudy( onePore, 'onePore_theta' + str(math.degrees(theta)) + 'IP' + str(intersectionParameter))

        facesListOnLeftSide = geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["FACE"], normalVectorOnLeftFace, middlePointOnLeftFace, GEOM.ST_ON)
        facesListOnRightSide = geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["FACE"], normalVectorOnRightFace, middlePointOnRightFace, GEOM.ST_ON)
        facesListOnOppositeSide = geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["FACE"], normalVectorOnOppositeFace, middlePointOnOppositeFace, GEOM.ST_ON)
        facesListOnFrontSide = geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["FACE"], normalVectorOnFrontFace, middlePointOnFrontFace, GEOM.ST_ON)
        facesListOnTopSide = geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["FACE"], normalVectorOnTopFace, middlePointOnTopFace, GEOM.ST_ON)
        facesListOnBottomSide = geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["FACE"], normalVectorOnBottomFace, middlePointOnBottomFace, GEOM.ST_ON)

        ## Get fases IDs without viscous layer
        allFacesIDsWithoutViscousLayers = []

        allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(onePore, geompy.ShapeType["FACE"], normalVectorOnLeftFace, middlePointOnLeftFace, GEOM.ST_ON))
        allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(onePore, geompy.ShapeType["FACE"], normalVectorOnRightFace, middlePointOnRightFace, GEOM.ST_ON))
        allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(onePore, geompy.ShapeType["FACE"], normalVectorOnOppositeFace, middlePointOnOppositeFace, GEOM.ST_ON))
        allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(onePore, geompy.ShapeType["FACE"], normalVectorOnFrontFace, middlePointOnFrontFace, GEOM.ST_ON))
        allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(onePore, geompy.ShapeType["FACE"], normalVectorOnTopFace, middlePointOnTopFace, GEOM.ST_ON))
        allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(onePore, geompy.ShapeType["FACE"], normalVectorOnBottomFace, middlePointOnBottomFace, GEOM.ST_ON))

        leftSide = geompy.CreateGroup(onePore, geompy.ShapeType["FACE"])
        rightSide = geompy.CreateGroup(onePore, geompy.ShapeType["FACE"])
        oppositeSide = geompy.CreateGroup(onePore, geompy.ShapeType["FACE"])
        frontSide = geompy.CreateGroup(onePore, geompy.ShapeType["FACE"])
        topSide = geompy.CreateGroup(onePore, geompy.ShapeType["FACE"])
        bottomSide = geompy.CreateGroup(onePore, geompy.ShapeType["FACE"])

        geompy.UnionList(leftSide, facesListOnLeftSide)
        geompy.UnionList(rightSide, facesListOnRightSide)
        geompy.UnionList(oppositeSide, facesListOnOppositeSide)
        geompy.UnionList(frontSide, facesListOnFrontSide)
        geompy.UnionList(topSide, facesListOnTopSide)
        geompy.UnionList(bottomSide, facesListOnBottomSide)

        geompy.addToStudyInFather( onePore, leftSide, 'leftSide' )
        geompy.addToStudyInFather( onePore, rightSide, 'rightSide' )
        geompy.addToStudyInFather( onePore, oppositeSide, 'oppositeSide' )
        geompy.addToStudyInFather( onePore, frontSide, 'frontSide' )
        geompy.addToStudyInFather( onePore, topSide, 'topSide' )
        geompy.addToStudyInFather( onePore, bottomSide, 'bottomSide' )

        ###
        ### SMESH component
        ###

        ####################################################
        ##         Begin of mesh variables section        ##
        ####################################################
        NETGEN_3D_MaxSize = 0.1 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)
        NETGEN_3D_MinSize = 0.001 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)

        NETGEN_2D_MaxSize = 0.1 * (geompy.BasicProperties(onePore)[1] - 6.0 * geompy.BasicProperties(leftSide)[1]) ** (1.0 / 2.0)
        NETGEN_2D_MinSize = 0.001 * (geompy.BasicProperties(onePore)[1] - 6.0 * geompy.BasicProperties(leftSide)[1]) ** (1.0 / 2.0)
        NETGEN_2D_MaxSize_1 = 0.1 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
        NETGEN_2D_MinSize_1 = 0.001 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
        NETGEN_2D_QuadAllowed = 0

        #viscousLayersSize = 0.005
        viscousLayersNumber = 2
        viscousLayersGrowth = 1.5

        Local_Length_Size = 0.05 *geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
        ####################################################
        ##          End of mesh variables section         ##
        ####################################################

        smesh = smeshBuilder.New(theStudy)
        onePoreMesh = smesh.Mesh(onePore)

        Regular_1D = onePoreMesh.Segment()
        NETGEN_2D_ONLY = onePoreMesh.Triangle(algo=smeshBuilder.NETGEN_2D)
        NETGEN_3D = onePoreMesh.Tetrahedron()

        #Global parameters for 1D, 2D and 3D mesh
        #Local_Length = Regular_1D.LocalLength( Local_Length_Size, None, 1e-6 )
        Deflection = Regular_1D.Deflection1D(4e-5)

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
            viscousLayersSize = 0.02 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
            Viscous_Layers = NETGEN_3D.ViscousLayers(viscousLayersSize,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.FACE_OFFSET)
            #Viscous_Layers = NETGEN_3D.ViscousLayers(viscousLayersSize,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.SURF_OFFSET_SMOOTH)
            #Viscous_Layers = NETGEN_3D.ViscousLayers(viscousLayersSize,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.NODE_OFFSET)
            
            ## Set names of Mesh objects
            #smesh.SetName(Viscous_Layers, 'Viscous Layers')

        #Sub-mesh for leftSide boundary
        Regular_1D_1 = onePoreMesh.Segment(geom=leftSide)
        NETGEN_2D_ONLY_1 = onePoreMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=leftSide)

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
        Projection_1D2D_1 = onePoreMesh.Projection1D2D(geom=rightSide)
        Source_Face_1 = Projection_1D2D_1.SourceFace(leftSide,None,None,None,None,None)

        #Sub-mesh for oppositeSide boundary
        Regular_1D_2 = onePoreMesh.Segment(geom=oppositeSide)
        NETGEN_2D_ONLY_2 = onePoreMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=oppositeSide)

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

        #status = onePoreMesh.AddHypothesis(NETGEN_Parameters_2D_ONLY_1,oppositeSide)

        #Projection oppositeSide's sub-mesh on frontSide boundary
        Projection_1D2D_2 = onePoreMesh.Projection1D2D(geom=frontSide)
        Source_Face_2 = Projection_1D2D_2.SourceFace(oppositeSide,None,None,None,None,None)

        #Sub-mesh for topSide boundary
        Regular_1D_3 = onePoreMesh.Segment(geom=topSide)
        NETGEN_2D_ONLY_3 = onePoreMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=topSide)

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

        #status = onePoreMesh.AddHypothesis(NETGEN_Parameters_2D_ONLY_1,topSide)

        #Projection topSide's sub-mesh on bottomSide boundary
        Projection_1D2D_3 = onePoreMesh.Projection1D2D(geom=bottomSide)
        Source_Face_3 = Projection_1D2D_3.SourceFace(topSide,None,None,None,None,None)

        #Compute the mesh
        isDone = onePoreMesh.Compute()

        leftSideBoundary = onePoreMesh.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
        rightSideBoundary = onePoreMesh.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
        oppositeSideBoundary = onePoreMesh.GroupOnGeom(oppositeSide,'oppositeSide',SMESH.FACE)
        frontSideBoundary = onePoreMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
        topSideBoundary = onePoreMesh.GroupOnGeom(topSide,'topSide',SMESH.FACE)
        bottomSideBoundary = onePoreMesh.GroupOnGeom(bottomSide,'bottomSide',SMESH.FACE)

        Sub_mesh_1 = NETGEN_2D_ONLY_1.GetSubMesh()
        Sub_mesh_2 = Projection_1D2D_1.GetSubMesh()
        Sub_mesh_3 = NETGEN_2D_ONLY_2.GetSubMesh()
        Sub_mesh_4 = Projection_1D2D_2.GetSubMesh()
        Sub_mesh_5 = NETGEN_2D_ONLY_3.GetSubMesh()
        Sub_mesh_6 = Projection_1D2D_3.GetSubMesh()

        ## Set names of Mesh objects
        #smesh.SetName(Local_Length, 'Local Length')
        #smesh.SetName(Deflection, 'Deflection')
        #smesh.SetName(Deflection_1, 'Deflection_1')
        #smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
        #smesh.SetName(Projection_1D2D_1.GetAlgorithm(), 'Projection_1D2D_1')
        #smesh.SetName(Projection_1D2D_2.GetAlgorithm(), 'Projection_1D2D_2')
        #smesh.SetName(Projection_1D2D_3.GetAlgorithm(), 'Projection_1D2D_3')
        #smesh.SetName(NETGEN_2D_ONLY.GetAlgorithm(), 'NETGEN_2D_ONLY')
        #smesh.SetName(NETGEN_2D_ONLY_1.GetAlgorithm(), 'NETGEN_2D_ONLY_1')

        #smesh.SetName(NETGEN_Parameters_3D, 'NETGEN 3D Parameters')
        #smesh.SetName(NETGEN_Parameters_2D_ONLY, 'NETGEN 2D Parameters')
        #smesh.SetName(NETGEN_Parameters_2D_ONLY_1, 'NETGEN 2D Parameters_1')

        #smesh.SetName(Source_Face_1, 'Source Face_1')
        #smesh.SetName(Source_Face_2, 'Source Face_2')
        #smesh.SetName(Source_Face_3, 'Source Face_3')
        
        #smesh.SetName(onePoreMesh.GetMesh(), 'onePoreMesh')

        #smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
        #smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
        #smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
        #smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
        #smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
        #smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser(True)
