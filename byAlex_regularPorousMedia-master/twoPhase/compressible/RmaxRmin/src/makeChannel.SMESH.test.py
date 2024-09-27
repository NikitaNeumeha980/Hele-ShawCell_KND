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
    print('Rmax=' + str(Rmax))

    for channelLength in [(channelLengthMin0 + x * channelLengthStep) for x in range(0, int((channelLengthMax0 - channelLengthMin0) / channelLengthStep) + 1)]:

        print('channelLength=' + str(channelLength))

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

        bottomFace = geompy.MakeRevolution2Ways(middleLineOnBottom, OZ, 0.5*math.pi/180.0)
        topFace = geompy.MakeRevolution2Ways(middleLineOnTop, OZ, 0.5*math.pi/180.0)

        middleFace = geompy.MakeQuad4Vertices(faceVertexOnBottom_1, faceVertexOnBottom_2, faceVertexOnTop_2, faceVertexOnTop_1)
        leftFace = geompy.MakeRotation(middleFace, OZ, 0.5*math.pi/180.0)
        rightFace = geompy.MakeRotation(middleFace, OZ, -0.5*math.pi/180.0)

        #For debuging
        #geompy.addToStudy( bottomFace, 'bottomFace' )
        #geompy.addToStudy( topFace, 'topFace' )
        #geompy.addToStudy( middleFace, 'middleFace' )
        #geompy.addToStudy( leftFace, 'leftFace' )
        #geompy.addToStudy( rightFace, 'rightFace' )

        #Make channel profile
        channelProfile = geompy.MakeCurveParametric("0.5 * " + str(Rmax + Rmin0) + "+ 0.5 * " + str(Rmax - Rmin0) + "* cos(2 * pi * t / " + str(channelLength) + ")", "0", "t", 0, channelLength, 20, GEOM.Interpolation, True)

        channelProfileFace = geompy.MakeRevolution2Ways(channelProfile, OZ, 0.5*math.pi/180.0)

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

        ###
        ### SMESH component
        ###

        smesh = smeshBuilder.New(theStudy)
        channelMesh = smesh.Mesh(channel2D)

        globalRegular1D = channelMesh.Segment()
        globalNumberOfSegments = globalRegular1D.NumberOfSegments(200)

        if (
            channelLength >= 10.0 * Rmin0 and
            Rmax / Rmin0 <= 3.0
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

        channelMesh.RotationSweepObjects( [ channelMesh ], [ channelMesh ], [ channelMesh ], SMESH.AxisStruct( 0, 0, 0, 0, 0, 1 ), 1.0*math.pi/180.0, 1, tolerance, 0 )
        channelMesh.RotateObject( channelMesh, SMESH.AxisStruct( 0, 0, 0, 0, 0, 1 ), -0.5*math.pi/180.0, 0 )

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


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
