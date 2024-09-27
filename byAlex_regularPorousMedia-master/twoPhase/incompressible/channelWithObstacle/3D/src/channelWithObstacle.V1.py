import sys
import math
import numpy
from decimal import *
getcontext().prec=12
import salome

import GEOM
from salome.geom import geomBuilder

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder
from salome.StdMeshers import StdMeshersBuilder

salome.salome_init()
theStudy = salome.myStudy

####################################################
##            Begin of variables section          ##
####################################################
absoluteCasePath_ = sys.argv[1]
drainageLength_ = Decimal('5000.0')
drainageWidth_ = Decimal('200.0')
channel1InOutWidth_ = Decimal('200.0')
channel1MiddleWidth_ = Decimal('300.0')
channel2InOutWidth_ = Decimal('200.0')
channel2MiddleWidth_ = Decimal('200.0')
depth_ = Decimal('200.0')
blobSeparatingDistance_ = Decimal('200.0')
blobRadius_ = Decimal('500.0')

slopeAngle_ = 10*Decimal(str(math.pi))/Decimal('180.0')
delta_ = depth_*Decimal(str(math.tan(slopeAngle_)))
#Not all filleting radius -> rework contour construction!
backContourFilletRadius_ = Decimal('45.0')
frontContourFilletRadius_ = Decimal('25.0')

a1_ =\
    channel1MiddleWidth_ + blobRadius_\
    if float(channel1MiddleWidth_) >= float(channel1InOutWidth_) else\
    (
        ((channel1InOutWidth_ + blobRadius_)**2 - (drainageWidth_/2)**2)
        /(1 - (drainageWidth_/2)**2/(channel1MiddleWidth_ + blobRadius_)**2)
    ).sqrt()

a2_ =\
    channel2MiddleWidth_ + blobRadius_\
    if float(channel2MiddleWidth_) >= float(channel2InOutWidth_) else\
    (
        ((channel2InOutWidth_ + blobRadius_)**2 - (drainageWidth_/2)**2)
        /(1 - (drainageWidth_/2)**2/(channel2MiddleWidth_ + blobRadius_)**2)
    ).sqrt()

b1_ =\
    channel1MiddleWidth_ + blobRadius_\
    if float(channel1MiddleWidth_) < float(channel1InOutWidth_) else\
    (
        ((channel1InOutWidth_ + blobRadius_)**2 - (drainageWidth_/2)**2)
        /(1 - (drainageWidth_/2)**2/(channel1MiddleWidth_ + blobRadius_)**2)
    ).sqrt()

b2_ =\
    channel2MiddleWidth_ + blobRadius_\
    if float(channel2MiddleWidth_) < float(channel2InOutWidth_) else\
    (
        ((channel2InOutWidth_ + blobRadius_)**2 - (drainageWidth_/2)**2)
        /(1 - (drainageWidth_/2)**2/(channel2MiddleWidth_ + blobRadius_)**2)
    ).sqrt()

print("a1_=" + str(a1_))
print("b1_=" + str(b1_))
print("a2_=" + str(a2_))
print("b2_=" + str(b2_))
####################################################
##             End of variables section           ##
####################################################


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

ellipseBackUp =\
    geompy.MakeEllipse(
        O,
        OZ,
        float(a1_),
        float(b1_),
        OY if float(channel1MiddleWidth_) >= float(channel1InOutWidth_) else OX
    )

ellipseBackDown =\
    geompy.MakeEllipse(
        O,
        OZ,
        float(a2_),
        float(b2_),
        OY if float(channel2MiddleWidth_) >= float(channel2InOutWidth_) else OX
    )

ellipseFrontUp =\
    geompy.MakeEllipse(
        O,
        OZ,
        float(a1_ + delta_),
        float(b1_ + delta_),
        OY if float(channel1MiddleWidth_) >= float(channel1InOutWidth_) else OX
    )

ellipseFrontDown =\
    geompy.MakeEllipse(
        O,
        OZ,
        float(a2_ + delta_),
        float(b2_ + delta_),
        OY if float(channel2MiddleWidth_) >= float(channel2InOutWidth_) else OX
    )

#geompy.addToStudy( ellipseBackUp, 'ellipseBackUp' )
#geompy.addToStudy( ellipseBackDown, 'ellipseBackDown' )

#geompy.addToStudy( ellipseFrontUp, 'ellipseFrontUp' )
#geompy.addToStudy( ellipseFrontDown, 'ellipseFrontDown' )

#asd =\
    #geompy.MakeFuseList(
        #[
            #geompy.MakeFace(ellipseBackUp, 1),
            #geompy.TranslateDXDYDZ(
                #geompy.MakeFaceHW(float(max(a1_, a2_) + drainageLength_), float(drainageWidth_), 1),
                #float((max(a1_, a2_) + drainageLength_)/2), 0, 0
            #)
        #],
        #False, False
    #)

#asd =\
    #geompy.MakeFillet2D(
        #asd,
        #float(backContourFilletRadius_),
        #[
            #geompy.GetSubShapeID(
                #asd,
                #geompy.GetVertexNearPoint(asd, geompy.MakeVertex(float(a1_), -float(drainageWidth_/2), 0))
            #),
            #geompy.GetSubShapeID(
                #asd,
                #geompy.GetVertexNearPoint(asd, geompy.MakeVertex(float(a2_), float(drainageWidth_/2), 0))
            #)
        #]
    #)

#geompy.addToStudy( asd, 'asd' )


backContour =\
    geompy.MakeWire(
        geompy.SubShapeAll(
            geompy.MakeFuseList(
                [
                    geompy.MakeCutList(
                        geompy.MakeFace(ellipseBackUp, 1),
                        [
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(float(2*a1_), float(2*a1_), 1),
                                0, -float(a1_), 0
                            ),
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(float(2*a1_), float(2*a1_), 1),
                                -float(a1_), 0, 0
                            )
                        ],
                        True
                    ),
                    geompy.MakeCutList(
                        geompy.MakeFace(ellipseBackDown, 1),
                        [
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(float(2*a2_), float(2*a2_), 1),
                                0, float(a2_), 0
                            ),
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(float(2*a2_), float(2*a2_), 1),
                                -float(a2_), 0, 0
                            )
                        ],
                        True
                    ),
                    geompy.TranslateDXDYDZ(
                        geompy.MakeFaceHW(float(max(a1_, a2_) + drainageLength_), float(drainageWidth_), 1),
                        float((max(a1_, a2_) + drainageLength_)/2), 0, 0
                    )
                ],
                True, True
            ),
            geompy.ShapeType["EDGE"]
        )
    )

filletedBackContour =\
    geompy.MakeFillet1D(
        backContour,
        float(backContourFilletRadius_),
        [
            geompy.GetSubShapeID(
                backContour,
                geompy.GetVertexNearPoint(backContour, geompy.MakeVertex(float(a1_), -float(drainageWidth_/2), 0))
            ),
            geompy.GetSubShapeID(
                backContour,
                geompy.GetVertexNearPoint(backContour, geompy.MakeVertex(float(a2_), float(drainageWidth_/2), 0))
            )
        ]
    )

backFace = geompy.MakeFace(filletedBackContour, 1)

frontContour =\
    geompy.TranslateDXDYDZ(
        geompy.MakeWire(
            geompy.SubShapeAll(
                geompy.MakeFuseList(
                    [
                        geompy.MakeCutList(
                            geompy.MakeFace(ellipseFrontUp, 1),
                            [
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(float(2*(a1_ + delta_)), float(2*(a1_ + delta_)), 1),
                                    0, -float(a1_ + delta_), 0
                                ),
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(float(2*(a1_ + delta_)), float(2*(a1_ + delta_)), 1),
                                    -float(a1_ + delta_), 0, 0
                                )
                            ],
                            True
                        ),
                        geompy.MakeCutList(
                            geompy.MakeFace(ellipseFrontDown, 1),
                            [
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(float(2*(a2_ + delta_)), float(2*(a2_ + delta_)), 1),
                                    0, float(a2_ + delta_), 0
                                ),
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(float(2*(a2_ + delta_)), float(2*(a2_ + delta_)), 1),
                                    -float(a2_ + delta_), 0, 0
                                )
                            ],
                            True
                        ),
                        geompy.TranslateDXDYDZ(
                            geompy.MakeFaceHW(float(max(a1_, a2_) + drainageLength_), float(drainageWidth_ + 2*delta_), 1),
                            float((max(a1_, a2_) + drainageLength_)/2), 0, 0
                        )
                    ],
                    True, True
                ),
                geompy.ShapeType["EDGE"]
            )
        ),
        0, 0, float(depth_)
    )

filletedFrontContour =\
    geompy.MakeFillet1D(
        frontContour,
        float(frontContourFilletRadius_),
        [
            geompy.GetSubShapeID(
                frontContour,
                geompy.GetVertexNearPoint(frontContour, geompy.MakeVertex(float(a1_), -float(drainageWidth_/2), float(depth_)))
            ),
            geompy.GetSubShapeID(
                frontContour,
                geompy.GetVertexNearPoint(frontContour, geompy.MakeVertex(float(a2_), float(drainageWidth_/2), float(depth_)))
            )
        ]
    )

frontFace = geompy.MakeFace(filletedFrontContour, 1)

splittedChannelRightPart =\
    geompy.TranslateDXDYDZ(
        geompy.MakeCut(
            geompy.MakePipeWithDifferentSections(
                [
                    backFace,
                    frontFace
                ],
                [
                    geompy.MakeVertex(0, 0, 0),
                    geompy.MakeVertex(0, 0, float(depth_))
                ],
                geompy.MakeLineTwoPnt(
                    geompy.MakeVertex(0, 0, 0),
                    geompy.MakeVertex(0, 0, float(depth_))
                ),
                0,
                0
            ),
            geompy.Rotate(
                geompy.MakeConeR1R2H(float(blobRadius_), float(blobRadius_ - delta_), float(depth_)),
                OZ,
                math.pi
            ),
            True
        ),
        float(blobSeparatingDistance_/2), 0, 0
    )

splittedChannelLeftPart =\
    geompy.MakeMirrorByPlane(
        splittedChannelRightPart,
        geompy.MakePlane(O, OX, 1000)
    )

splittedChannel =\
    geompy.MakeFuseList(
        [
            splittedChannelRightPart,
            splittedChannelLeftPart,
            geompy.TranslateDXDYDZ(
                geompy.MakeHexa2Faces(
                    geompy.MakeFaceHW(float(blobSeparatingDistance_), float(channel1MiddleWidth_), 1),
                    geompy.TranslateDXDYDZ(
                        geompy.MakeFaceHW(float(blobSeparatingDistance_), float(channel1MiddleWidth_ + 2*delta_), 1),
                        0, 0, float(depth_)
                    )
                ),
                0, float(blobRadius_ + channel1MiddleWidth_/2), 0
            ),
            geompy.TranslateDXDYDZ(
                geompy.MakeHexa2Faces(
                    geompy.MakeFaceHW(float(blobSeparatingDistance_), float(channel2MiddleWidth_), 1),
                    geompy.TranslateDXDYDZ(
                        geompy.MakeFaceHW(float(blobSeparatingDistance_), float(channel2MiddleWidth_ + 2*delta_), 1),
                        0, 0, float(depth_)
                    )
                ),
                0, -float(blobRadius_ + channel2MiddleWidth_/2), 0
            )
        ],
        True, True
    )

#splittedChannel =\
    #geompy.MakePartition(
        #[
            #splittedChannelRightPart,
            #splittedChannelLeftPart,
            #geompy.TranslateDXDYDZ(
                #geompy.MakeHexa2Faces(
                    #geompy.MakeFaceHW(blobSeparatingDistance_, channel1MiddleWidth_, 1),
                    #geompy.TranslateDXDYDZ(
                        #geompy.MakeFaceHW(blobSeparatingDistance_, channel1MiddleWidth_ + 2.0*delta_, 1),
                        #0, 0, depth_
                    #)
                #),
                #0, blobRadius_ + 0.5*channel1MiddleWidth_, 0
            #),
            #geompy.TranslateDXDYDZ(
                #geompy.MakeHexa2Faces(
                    #geompy.MakeFaceHW(blobSeparatingDistance_, channel2MiddleWidth_, 1),
                    #geompy.TranslateDXDYDZ(
                        #geompy.MakeFaceHW(blobSeparatingDistance_, channel2MiddleWidth_ + 2.0*delta_, 1),
                        #0, 0, depth_
                    #)
                #),
                #0, -(blobRadius_ + 0.5*channel2MiddleWidth_), 0
            #)
        #]
    #)

##Create groups for walls extraction
leftSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
rightSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
#frontSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
#backSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
walls = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])

bb = geompy.BoundingBox(splittedChannel)

#bbBox =\
    #geompy.MakeBoxTwoPnt(
        #geompy.MakeVertex(bb[0], bb[2], bb[4]),
        #geompy.MakeVertex(bb[1], bb[3], bb[5])
    #)

#geompy.addToStudy( bbBox, 'bbBox' )

geompy.UnionList(
    leftSide,
    geompy.GetShapesOnShape(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(bb[0], bb[2], bb[4]),
                geompy.MakeVertex(bb[1], bb[3], bb[5])
            ),
            -abs(bb[1] - bb[0]), 0, 0
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
geompy.UnionList(
    rightSide,
    geompy.GetShapesOnShape(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(bb[0], bb[2], bb[4]),
                geompy.MakeVertex(bb[1], bb[3], bb[5])
            ),
            abs(bb[1] - bb[0]), 0, 0
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
#geompy.UnionList(
    #frontSide,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxTwoPnt(
                #geompy.MakeVertex(bb[0], bb[2], bb[4]),
                #geompy.MakeVertex(bb[1], bb[3], bb[5])
            #),
            #0, 0, abs(bb[5] - bb[4])
        #),
        #splittedChannel,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
#geompy.UnionList(
    #backSide,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxTwoPnt(
                #geompy.MakeVertex(bb[0], bb[2], bb[4]),
                #geompy.MakeVertex(bb[1], bb[3], bb[5])
            #),
            #0, 0, -abs(bb[5] - bb[4])
        #),
        #splittedChannel,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
##Walls
geompy.UnionList(
    walls,
    geompy.SubShapeAll(
        splittedChannel,
        geompy.ShapeType["FACE"]
    )
)
geompy.DifferenceList(walls, geompy.SubShapeAll(leftSide, geompy.ShapeType["FACE"]))
geompy.DifferenceList(walls, geompy.SubShapeAll(rightSide, geompy.ShapeType["FACE"]))
#geompy.DifferenceList(walls, geompy.SubShapeAll(frontSide, geompy.ShapeType["FACE"]))
#geompy.DifferenceList(walls, geompy.SubShapeAll(backSide, geompy.ShapeType["FACE"]))

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( ellipseBackUp, 'ellipseBackUp' )
geompy.addToStudy( ellipseBackDown, 'ellipseBackDown' )
geompy.addToStudy( ellipseFrontUp, 'ellipseFrontUp' )
geompy.addToStudy( ellipseFrontDown, 'ellipseFrontDown' )
geompy.addToStudy( backContour, 'backContour' )
geompy.addToStudy( filletedBackContour, 'filletedBackContour' )
geompy.addToStudy( frontContour, 'frontContour' )
geompy.addToStudy( filletedFrontContour, 'filletedFrontContour' )

geompy.addToStudy( splittedChannelRightPart, 'splittedChannelRightPart' )
geompy.addToStudy( splittedChannelLeftPart, 'splittedChannelLeftPart' )
geompy.addToStudy( splittedChannel, 'splittedChannel' )
geompy.addToStudyInFather( splittedChannel, leftSide, 'leftSide' )
geompy.addToStudyInFather( splittedChannel, rightSide, 'rightSide' )
#geompy.addToStudyInFather( splittedChannel, frontSide, 'frontSide' )
#geompy.addToStudyInFather( splittedChannel, backSide, 'backSide' )
geompy.addToStudyInFather( splittedChannel, walls, 'walls' )


###
### SMESH component
###

smesh = smeshBuilder.New()
##smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 ## multiples meshes built in parallel, complex and numerous mesh edition (performance)

splittedChannelMesh = smesh.Mesh(splittedChannel)

splittedChannelMesh.Segment().Adaptive(
    0.02*float(
        min(
            drainageLength_,
            drainageWidth_,
            channel1MiddleWidth_,
            channel2MiddleWidth_,
            channel1InOutWidth_,
            channel2InOutWidth_,
            depth_,
            blobSeparatingDistance_,
            blobRadius_
        )
    ),
    0.1*float(
        min(
            drainageLength_,
            drainageWidth_,
            channel1MiddleWidth_,
            channel2MiddleWidth_,
            channel1InOutWidth_,
            channel2InOutWidth_,
            depth_,
            blobSeparatingDistance_,
            blobRadius_
        )
    ),
    2.0
)

NETGEN_2D_Parameters = splittedChannelMesh.Triangle(algo=smeshBuilder.NETGEN_2D).Parameters()
NETGEN_2D_Parameters.SetMaxSize(
    float(
        min(
            drainageLength_,
            drainageWidth_,
            channel1MiddleWidth_,
            channel2MiddleWidth_,
            channel1InOutWidth_,
            channel2InOutWidth_,
            depth_,
            blobSeparatingDistance_,
            blobRadius_
        )
    )
)
NETGEN_2D_Parameters.SetMinSize(
    0.02*float(
        min(
            drainageLength_,
            drainageWidth_,
            channel1MiddleWidth_,
            channel2MiddleWidth_,
            channel1InOutWidth_,
            channel2InOutWidth_,
            depth_,
            blobSeparatingDistance_,
            blobRadius_
        )
    )
)
NETGEN_2D_Parameters.SetOptimize( 1 )
NETGEN_2D_Parameters.SetFineness( 5 )
NETGEN_2D_Parameters.SetGrowthRate( 0.05 )
NETGEN_2D_Parameters.SetChordalError( -1 )
NETGEN_2D_Parameters.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters.SetWorstElemMeasure( 21851 )
NETGEN_2D_Parameters.SetUseDelauney( 0 )
NETGEN_2D_Parameters.SetCheckChartBoundary( 8 )
NETGEN_2D_Parameters.SetQuadAllowed( 0 )

NETGEN_3D = splittedChannelMesh.Tetrahedron()
NETGEN_3D_Parameters = NETGEN_3D.Parameters()
NETGEN_3D_Parameters.SetMaxSize(
    float(
        min(
            drainageLength_,
            drainageWidth_,
            channel1MiddleWidth_,
            channel2MiddleWidth_,
            channel1InOutWidth_,
            channel2InOutWidth_,
            depth_,
            blobSeparatingDistance_,
            blobRadius_
        )
    )
)
NETGEN_3D_Parameters.SetMinSize(
    0.02*float(
        min(
            drainageLength_,
            drainageWidth_,
            channel1MiddleWidth_,
            channel2MiddleWidth_,
            channel1InOutWidth_,
            channel2InOutWidth_,
            depth_,
            blobSeparatingDistance_,
            blobRadius_
        )
    )
)
NETGEN_3D_Parameters.SetOptimize( 1 )
NETGEN_3D_Parameters.SetFineness( 5 )
NETGEN_3D_Parameters.SetGrowthRate( 0.05 )
NETGEN_3D_Parameters.SetElemSizeWeight( 0 )
NETGEN_3D_Parameters.SetCheckOverlapping( 0 )
NETGEN_3D_Parameters.SetCheckChartBoundary( 8 )
Viscous_Layers = NETGEN_3D.ViscousLayers(
    5, 2, 1,
    [
        geompy.GetSubShapeID(
            splittedChannel,
            leftSide
        ),
        geompy.GetSubShapeID(
            splittedChannel,
            rightSide
        )
    ],
    1,
    smeshBuilder.FACE_OFFSET
    #smeshBuilder.NODE_OFFSET
)

isDone = splittedChannelMesh.Compute()

splittedChannelMesh.GroupOnGeom(leftSide,'left',SMESH.FACE)
splittedChannelMesh.GroupOnGeom(rightSide,'right',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(backSide,'backSide',SMESH.FACE)
splittedChannelWallsMesh = splittedChannelMesh.GroupOnGeom(walls,'walls',SMESH.FACE)

try:
  splittedChannelMesh.ExportUNV( absoluteCasePath_ + "/splittedChannelMesh.unv" )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')

try:
  splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/walls.stl", 1, splittedChannelWallsMesh)
  pass
except:
  print('ExportSTL() failed. Invalid file name?')

salome.myStudy.SaveAs(absoluteCasePath_ + "/splittedChannel.hdf", False, False)





#NbS = splittedChannelMesh.Segment().NumberOfSegments(15, None, [])
#NbS.SetConversionMode( 1 )
#NbS.SetTableFunction(
    #[
        #0, 1,
        #0.1, 0.1,
        #0.9, 0.1,
        #1, 1
    #]
#)

#splittedChannelMesh.Prism()

#NETGEN_1D_2D = splittedChannelMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=frontSide)
#NETGEN_2D_Parameters = NETGEN_1D_2D.Parameters()
#NETGEN_2D_Parameters.SetMaxSize(
    #min(
        #drainageLength_,
        #drainageWidth_,
        #channel1MiddleWidth_,
        #channel2MiddleWidth_,
        #channel1InOutWidth_,
        #channel2InOutWidth_,
        #depth_,
        #blobSeparatingDistance_,
        #blobRadius_
    #)
#)
#NETGEN_2D_Parameters.SetMinSize(
    #0.02*min(
        #drainageLength_,
        #drainageWidth_,
        #channel1MiddleWidth_,
        #channel2MiddleWidth_,
        #channel1InOutWidth_,
        #channel2InOutWidth_,
        #depth_,
        #blobSeparatingDistance_,
        #blobRadius_
    #)
#)
#NETGEN_2D_Parameters.SetSecondOrder( 0 )
#NETGEN_2D_Parameters.SetOptimize( 1 )
#NETGEN_2D_Parameters.SetFineness( 5 )
#NETGEN_2D_Parameters.SetGrowthRate( 0.1 )
#NETGEN_2D_Parameters.SetNbSegPerEdge( 10 )
#NETGEN_2D_Parameters.SetNbSegPerRadius( 10 )
#NETGEN_2D_Parameters.SetChordalError( -1 )
#NETGEN_2D_Parameters.SetChordalErrorEnabled( 0 )
#NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
#NETGEN_2D_Parameters.SetFuseEdges( 1 )
#NETGEN_2D_Parameters.SetWorstElemMeasure( 22083 )
#NETGEN_2D_Parameters.SetUseDelauney( 0 )
#NETGEN_2D_Parameters.SetCheckChartBoundary( 8 )
#NETGEN_2D_Parameters.SetQuadAllowed( 0 )
#Viscous_Layers_2D_1 = NETGEN_1D_2D.ViscousLayers2D(5, 1, 1, [], 1, 'Viscous Layers')

#splittedChannelMesh.Projection1D2D(geom=backSide).SourceFace(frontSide, None, None, None, None, None)

#isDone = splittedChannelMesh.Compute()

#splittedChannelMesh.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(backSide,'backSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(walls,'walls',SMESH.FACE)




#splittedChannelMesh = smesh.Mesh(splittedChannel)
#NETGEN_1D_2D = splittedChannelMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
#NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
#NETGEN_2D_Parameters_1.SetMaxSize(
    #min(
        #drainageLength_,
        #drainageWidth_,
        #channel1MiddleWidth_,
        #channel2MiddleWidth_,
        #channel1InOutWidth_,
        #channel2InOutWidth_,
        #depth_,
        #blobSeparatingDistance_,
        #blobRadius_
    #)
#)
#NETGEN_2D_Parameters_1.SetMinSize(
    #0.02*min(
        #drainageLength_,
        #drainageWidth_,
        #channel1MiddleWidth_,
        #channel2MiddleWidth_,
        #channel1InOutWidth_,
        #channel2InOutWidth_,
        #depth_,
        #blobSeparatingDistance_,
        #blobRadius_
    #)
#)
#NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
#NETGEN_2D_Parameters_1.SetOptimize( 1 )
#NETGEN_2D_Parameters_1.SetFineness( 5 )
#NETGEN_2D_Parameters_1.SetGrowthRate( 0.5 )
#NETGEN_2D_Parameters_1.SetNbSegPerEdge( 3 )
#NETGEN_2D_Parameters_1.SetNbSegPerRadius( 10 )
#NETGEN_2D_Parameters_1.SetChordalError( -1 )
#NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
#NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
#NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
#NETGEN_2D_Parameters_1.SetWorstElemMeasure( 21874 )
#NETGEN_2D_Parameters_1.SetUseDelauney( 0 )
#NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
#NETGEN_2D_Parameters_1.SetCheckChartBoundary( 184 )
#isDone = splittedChannelMesh.Compute()
#splittedChannelWalls = splittedChannelMesh.GroupOnGeom(walls,'walls',SMESH.FACE)

#try:
  #splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/splittedChannel.stl", 1 )
  #splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/splittedChannelWalls.stl", 1, splittedChannelWalls)
  #pass
#except:
  #print('ExportSTL() failed. Invalid file name?')

#salome.myStudy.SaveAs(absoluteCasePath_ + "/splittedChannel.hdf", False, False)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
