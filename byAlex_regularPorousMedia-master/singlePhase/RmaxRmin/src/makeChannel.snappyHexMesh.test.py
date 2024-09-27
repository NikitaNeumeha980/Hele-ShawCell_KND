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
channelLengthMin0 = 10.0
channelLengthMax0 = 100.0
channelLengthStep = 100.0

Rmin0 = 1.0
Rmax0 = 10.0 * Rmin0
Rstep = 1.0

numberOfSegmentsAlongChannelLength0 = 20
#numberOfSegmentsPerUnitChannelLength = 1

tolerance = 1e-07
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

for Rmax in [(Rmin0 + x * Rstep) for x in range(0, int((Rmax0 - Rmin0) / Rstep) + 1)]:

    #For debuging
    print("***********************************")
    print('Rmax=' + str(Rmax))

    for channelLength in [(channelLengthMin0 + x * channelLengthStep) for x in range(0, int((channelLengthMax0 - channelLengthMin0) / channelLengthStep) + 1)]:

        print('channelLength=' + str(channelLength))

        backgroundMeshBoxLength = channelLength
        backgroundMeshBoxWidth = 1.1 * Rmax

        delta = min(backgroundMeshBoxWidth, channelLength / numberOfSegmentsAlongChannelLength0)
        numberOfSegmentsAlongChannelLength = int(channelLength / delta)

        backgroundMeshBoxDepth = channelLength / numberOfSegmentsAlongChannelLength

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

        ##Make middle points on boundaries for searching boundaries
        #backgroundMeshBoxMiddlePointOnLeftFace = geompy.MakeVertexOnSurface(backgroundMeshBoxLeftFace, 0.5, 0.5)
        #backgroundMeshBoxMiddlePointOnRightFace = geompy.MakeVertexOnSurface(backgroundMeshBoxRightFace, 0.5, 0.5)
        #backgroundMeshBoxMiddlePointOnOppositeFace = geompy.MakeVertexOnSurface(backgroundMeshBoxOppositeFace, 0.5, 0.5)
        #backgroundMeshBoxMiddlePointOnFrontFace = geompy.MakeVertexOnSurface(backgroundMeshBoxFrontFace, 0.5, 0.5)
        #backgroundMeshBoxMiddlePointOnTopFace = geompy.MakeVertexOnSurface(backgroundMeshBoxTopFace, 0.5, 0.5)
        #backgroundMeshBoxMiddlePointOnBottomFace = geompy.MakeVertexOnSurface(backgroundMeshBoxBottomFace, 0.5, 0.5)

        ##Make normals for searching boundaries
        #backgroundMeshBoxNormalVectorOnLeftFace = geompy.GetNormal(backgroundMeshBoxLeftFace)
        #backgroundMeshBoxNormalVectorOnRightFace = geompy.GetNormal(backgroundMeshBoxRightFace)
        #backgroundMeshBoxNormalVectorOnOppositeFace = geompy.GetNormal(backgroundMeshBoxOppositeFace)
        #backgroundMeshBoxNormalVectorOnFrontFace = geompy.GetNormal(backgroundMeshBoxFrontFace)
        #backgroundMeshBoxNormalVectorOnTopFace = geompy.GetNormal(backgroundMeshBoxTopFace)
        #backgroundMeshBoxNormalVectorOnBottomFace = geompy.GetNormal(backgroundMeshBoxBottomFace)

        backgroundMeshBox = geompy.MakeHexa(backgroundMeshBoxLeftFace, backgroundMeshBoxRightFace, backgroundMeshBoxOppositeFace, backgroundMeshBoxFrontFace, backgroundMeshBoxTopFace, backgroundMeshBoxBottomFace)

        edgeForSubMeshOnDepth = geompy.GetEdge(backgroundMeshBox, backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_4)
        edgeForSubMeshOnWidth = geompy.GetEdge(backgroundMeshBox, backgroundMeshBoxVertexOnBottom_1, backgroundMeshBoxVertexOnBottom_2)

        #For debuging
        #geompy.addToStudy( backgroundMeshBoxLeftFace, 'backgroundMeshBoxLeftFace' )
        #geompy.addToStudy( backgroundMeshBoxRightFace, 'backgroundMeshBoxRightFace' )
        #geompy.addToStudy( backgroundMeshBoxOppositeFace, 'backgroundMeshBoxOppositeFace' )
        #geompy.addToStudy( backgroundMeshBoxFrontFace, 'backgroundMeshBoxFrontFace' )
        #geompy.addToStudy( backgroundMeshBoxTopFace, 'backgroundMeshBoxTopFace' )
        #geompy.addToStudy( backgroundMeshBoxBottomFace, 'backgroundMeshBoxBottomFace' )

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

        #middleLineOnBottom = geompy.MakeLineTwoPnt(faceVertexOnBottom_1, faceVertexOnBottom_2)
        #middleLineOnTop = geompy.MakeLineTwoPnt(faceVertexOnTop_1, faceVertexOnTop_2)

        #bottomFace = geompy.MakeRevolution2Ways(middleLineOnBottom, OZ, 0.5*sectorAngle*math.pi/180.0)
        #topFace = geompy.MakeRevolution2Ways(middleLineOnTop, OZ, 0.5*sectorAngle*math.pi/180.0)

        #middleFace = geompy.MakeQuad4Vertices(faceVertexOnBottom_1, faceVertexOnBottom_2, faceVertexOnTop_2, faceVertexOnTop_1)
        #leftFace = geompy.MakeRotation(middleFace, OZ, 0.5*sectorAngle*math.pi/180.0)
        #rightFace = geompy.MakeRotation(middleFace, OZ, -0.5*sectorAngle*math.pi/180.0)

        #For debuging
        #geompy.addToStudy( bottomFace, 'bottomFace' )
        #geompy.addToStudy( topFace, 'topFace' )
        #geompy.addToStudy( middleFace, 'middleFace' )
        #geompy.addToStudy( leftFace, 'leftFace' )
        #geompy.addToStudy( rightFace, 'rightFace' )

        #Make channel profile
        channelProfile1D = geompy.MakeCurveParametric("0.5 * " + str(Rmax + Rmin0) + "+ 0.5 * " + str(Rmax - Rmin0) + "* cos(2 * pi * t / " + str(channelLength) + ")", "0", "t", 0, channelLength, 100, GEOM.Interpolation, True)

        channelProfile2D = geompy.MakePrismVecH(channelProfile1D, OY, backgroundMeshBoxDepth)

        #channelProfileFace = geompy.MakeRevolution2Ways(channelProfile1D, OZ, 0.5*sectorAngle*math.pi/180.0)

        sk = geompy.Sketcher3D()
        sk.addPointsAbsolute(Rmax, 0, 0)
        sk.addPointsAbsolute(0, 0, 0)
        sk.addPointsAbsolute(0, 0, channelLength)
        sk.addPointsAbsolute(Rmax, 0, channelLength)
        closingContour = sk.wire()

        channel2D = geompy.MakeFaceWires([channelProfile1D, closingContour], 1)

        planeChannel3D = geompy.MakePrismVecH2Ways(channel2D, OY, 2.0 * backgroundMeshBoxDepth)
        
        geompy.addToStudy( planeChannel3D, 'planeChannel3DRmax' + str(Rmax) + 'channelLength' + str(channelLength) )

        geompy.addToStudyInFather( planeChannel3D, channelProfile2D, 'channelProfile2D' )

        geompy.addToStudyInFather( planeChannel3D, channel2D, 'channel2D' )

        geompy.addToStudyInFather( planeChannel3D, backgroundMeshBox, 'backgroundMeshBox' )

        geompy.addToStudyInFather( backgroundMeshBox, edgeForSubMeshOnDepth, 'edgeForSubMeshOnDepth' )
        geompy.addToStudyInFather( backgroundMeshBox, edgeForSubMeshOnWidth, 'edgeForSubMeshOnWidth' )

        geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxBottomFace, 'backgroundMeshBoxBottomFace' )
        geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxTopFace, 'backgroundMeshBoxTopFace' )
        geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxOppositeFace, 'backgroundMeshBoxOppositeFace' )
        geompy.addToStudyInFather( backgroundMeshBox, backgroundMeshBoxFrontFace, 'backgroundMeshBoxFrontFace' )

        #geompy.addToStudyInFather( channel2D, bottomFace, 'bottomFace' )
        #geompy.addToStudyInFather( channel2D, topFace, 'topFace' )
        #geompy.addToStudyInFather( channel2D, middleFace, 'middleFace' )
        #geompy.addToStudyInFather( channel2D, leftFace, 'leftFace' )
        #geompy.addToStudyInFather( channel2D, rightFace, 'rightFace' )
        #geompy.addToStudyInFather( channel2D, channelProfileFace, 'channelProfileFace' )

        ##Save geometry to file
        #print("Save geometry to " + ABSOLUTE_CASE_PATH + "/channelGeometry.hdf\n")

        #salome.myStudyManager.SaveAs(ABSOLUTE_CASE_PATH + "/channelGeometry.hdf", salome.myStudy, False)

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

        print("done")
        print("***********************************\n")

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser(True)
