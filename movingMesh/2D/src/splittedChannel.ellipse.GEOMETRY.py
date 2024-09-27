import sys
import math
import numpy
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
drainageLength_ = 5000
drainageWidth_ = 200
channel1InOutWidth_ = 100
channel1MiddleWidth_ = 90
channel2InOutWidth_ = 100
channel2MiddleWidth_ = 100
depth_ = 1000
blobSeparatingDistance_ = 100
blobRadius_ = 100

slopeAngle_ = 0.0*math.pi/180.0
delta_ = depth_*math.tan(slopeAngle_)
backContourFilletRadius_ = 50
frontContourFilletRadius_ = 50

a1_ =\
    channel1MiddleWidth_ + blobRadius_\
    if channel1MiddleWidth_ >= channel1InOutWidth_ else\
    round(math.sqrt((channel1InOutWidth_ + blobRadius_)**2 - 0.25*drainageWidth_**2)/math.sqrt(1.0 - 0.25*drainageWidth_**2/(channel1MiddleWidth_ + blobRadius_)**2), 6)

a2_ =\
    channel2MiddleWidth_ + blobRadius_\
    if channel2MiddleWidth_ >= channel2InOutWidth_ else\
    round(math.sqrt((channel2InOutWidth_ + blobRadius_)**2 - 0.25*drainageWidth_**2)/math.sqrt(1.0 - 0.25*drainageWidth_**2/(channel2MiddleWidth_ + blobRadius_)**2), 6)

b1_ =\
    channel1MiddleWidth_ + blobRadius_\
    if channel1MiddleWidth_ < channel1InOutWidth_ else\
    round(math.sqrt((channel1InOutWidth_ + blobRadius_)**2 - 0.25*drainageWidth_**2)/math.sqrt(1.0 - 0.25*drainageWidth_**2/(channel1MiddleWidth_ + blobRadius_)**2), 6)

b2_ =\
    channel2MiddleWidth_ + blobRadius_\
    if channel2MiddleWidth_ < channel2InOutWidth_ else\
    round(math.sqrt((channel2InOutWidth_ + blobRadius_)**2 - 0.25*drainageWidth_**2)/math.sqrt(1.0 - 0.25*drainageWidth_**2/(channel2MiddleWidth_ + blobRadius_)**2), 6)


####################################################
##             For topoSet                        ##

#with open(absoluteCasePath_ + "/forTopoSetBox", 'w') as forTopoSetBox:
    #forTopoSetBox.write(
        #"boxes (("
        #+ str(-max(b1_,b2_)) + " " + 0 + " " + str(-depth_) +
        #")("
        #+ str(max(b1_,b2_)) + " " + str(a1_) + " " + str(depth_)+
        #") ("
        #+ str(-min(b1_,b2_)) + " " + str(-a2_) + " " + str(-depth_) +
        #") ("
        #+ + str(-min(b1_,b2_)) + " " + 0 + " " + str(-depth_) + "));"

with open(absoluteCasePath_ + "/forTopoSetBox", 'w') as forTopoSetBox:
    forTopoSetBox.write(
        "box ("
        + str(-(max(b1_,b2_)+0.5*blobSeparatingDistance_)) + " " + str(-a2_) + " " + str(-depth_) +
        ")("
        + str(max(b1_,b2_)+0.5*blobSeparatingDistance_) + " " + str(a1_) + " " + str(depth_)+
        ");"
    )

####################################################
##             For setFields                      ##

with open(absoluteCasePath_ + "/forSetFieldsBox", 'w') as forSetFieldsBox:
    forSetFieldsBox.write(
        str(max(b1_,b2_)+0.5*blobSeparatingDistance_ + 100)
    )

####################################################
##             For snappyHexMesh                  ##

with open(absoluteCasePath_ + "/system/snappyHexMeshSrc/refinementBox", 'w') as refinementBox:
    refinementBox.write(
            "min (" + str(-(max(b1_,b2_)+0.5*blobSeparatingDistance_)) + " " + str(-a2_) + " " + str(-depth_) + ");" + "\n" +
            "max (" + str(max(b1_,b2_)+0.5*blobSeparatingDistance_) + " " + str(a1_) + " " + str(depth_)+ ");"
    )

####################################################

#print("a1_=" + str(a1_))
#print("b1_=" + str(b1_))
#print("a2_=" + str(a2_))
#print("b2_=" + str(b2_))
####################################################
##             End of variables section           ##
####################################################


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

ellipseBackUp =\
    geompy.MakeEllipse(O, OZ, a1_, b1_, OY)\
    if channel1MiddleWidth_ >= channel1InOutWidth_ else\
    geompy.MakeEllipse(O, OZ, a1_, b1_, OX)

ellipseBackDown =\
    geompy.MakeEllipse(O, OZ, a2_, b2_, OY)\
    if channel2MiddleWidth_ >= channel2InOutWidth_ else\
    geompy.MakeEllipse(O, OZ, a2_, b2_, OX)

ellipseFrontUp =\
    geompy.MakeEllipse(O, OZ, a1_ + delta_, b1_ + delta_, OY)\
    if channel1MiddleWidth_ >= channel1InOutWidth_ else\
    geompy.MakeEllipse(O, OZ, a1_ + delta_, b1_ + delta_, OX)

ellipseFrontDown =\
    geompy.MakeEllipse(O, OZ, a2_ + delta_, b2_ + delta_, OY)\
    if channel2MiddleWidth_ >= channel2InOutWidth_ else\
    geompy.MakeEllipse(O, OZ, a2_ + delta_, b2_ + delta_, OX)

geompy.addToStudy( ellipseBackUp, 'ellipseBackUp' )
geompy.addToStudy( ellipseBackDown, 'ellipseBackDown' )

geompy.addToStudy( ellipseFrontUp, 'ellipseFrontUp' )
geompy.addToStudy( ellipseFrontDown, 'ellipseFrontDown' )

backContour =\
    geompy.MakeWire(
        geompy.SubShapeAll(
            geompy.MakeFuseList(
                [
                    geompy.MakeCutList(
                        geompy.MakeFace(ellipseBackUp, 1),
                        [
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(2.0*a1_, 2.0*a1_, 1),
                                0, -a1_, 0
                            ),
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(2.0*a1_, 2.0*a1_, 1),
                                -a1_, 0, 0
                            )
                        ],
                        True
                    ),
                    geompy.MakeCutList(
                        geompy.MakeFace(ellipseBackDown, 1),
                        [
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(2.0*a2_, 2.0*a2_, 1),
                                0, a2_, 0
                            ),
                            geompy.MakeTranslation(
                                geompy.MakeFaceHW(2.0*a2_, 2.0*a2_, 1),
                                -a2_, 0, 0
                            )
                        ],
                        True
                    ),
                    geompy.TranslateDXDYDZ(
                        geompy.MakeFaceHW(max(a1_, a2_) + drainageLength_, drainageWidth_, 1),
                        0.5*(max(a1_, a2_) + drainageLength_), 0, 0
                    )
                ],
                True,
                True
            ),
            geompy.ShapeType["EDGE"]
        )
    )

filletedBackContour =\
    geompy.MakeFillet1D(
        backContour,
        backContourFilletRadius_,
        [
            geompy.GetSubShapeID(
                backContour,
                geompy.GetVertexNearPoint(backContour, geompy.MakeVertex(a1_, -0.5*drainageWidth_, 0))
            ),
            geompy.GetSubShapeID(
                backContour,
                geompy.GetVertexNearPoint(backContour, geompy.MakeVertex(a2_, 0.5*drainageWidth_, 0))
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
                                    geompy.MakeFaceHW(2.0*(a1_ + delta_), 2.0*(a1_ + delta_), 1),
                                    0, -(a1_ + delta_), 0
                                ),
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(2.0*(a1_ + delta_), 2.0*(a1_ + delta_), 1),
                                    -(a1_ + delta_), 0, 0
                                )
                            ],
                            True
                        ),
                        geompy.MakeCutList(
                            geompy.MakeFace(ellipseFrontDown, 1),
                            [
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(2.0*(a2_ + delta_), 2.0*(a2_ + delta_), 1),
                                    0, a2_ + delta_, 0
                                ),
                                geompy.MakeTranslation(
                                    geompy.MakeFaceHW(2.0*(a2_ + delta_), 2.0*(a2_ + delta_), 1),
                                    -(a2_ + delta_), 0, 0
                                )
                            ],
                            True
                        ),
                        geompy.TranslateDXDYDZ(
                            geompy.MakeFaceHW(max(a1_, a2_) + drainageLength_, drainageWidth_ + 2.0*delta_, 1),
                            0.5*(max(a1_, a2_) + drainageLength_), 0, 0
                        )
                    ],
                    True,
                    True
                ),
                geompy.ShapeType["EDGE"]
            )
        ),
        0, 0, depth_
    )

filletedFrontContour =\
    geompy.MakeFillet1D(
        frontContour,
        frontContourFilletRadius_,
        [
            geompy.GetSubShapeID(
                frontContour,
                geompy.GetVertexNearPoint(frontContour, geompy.MakeVertex(a1_, -0.5*drainageWidth_, depth_))
            ),
            geompy.GetSubShapeID(
                frontContour,
                geompy.GetVertexNearPoint(frontContour, geompy.MakeVertex(a2_, 0.5*drainageWidth_, depth_))
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
                    geompy.MakeVertex(0, 0, depth_)
                ],
                geompy.MakeLineTwoPnt(
                    geompy.MakeVertex(0, 0, 0),
                    geompy.MakeVertex(0, 0, depth_)
                ),
                0,
                0
            ),
            geompy.Rotate(
                geompy.MakeConeR1R2H(blobRadius_, blobRadius_ - delta_, depth_),
                OZ,
                math.pi
            ),
            True
        ),
        0.5*blobSeparatingDistance_,
        0,
        -0.5*depth_
    )

splittedChannelLeftPart =\
    geompy.MakeMirrorByPlane(
        splittedChannelRightPart,
        geompy.MakePlane(O, OX, 1000)
    )

if (blobSeparatingDistance_ <= 0):
    splittedChannel =\
        geompy.MakeFuseList(
            [
                splittedChannelRightPart,
                splittedChannelLeftPart
            ],
            False,
            False
        )
else:
    splittedChannel =\
        geompy.MakeFuseList(
            [
                splittedChannelRightPart,
                splittedChannelLeftPart,
                geompy.TranslateDXDYDZ(
                    geompy.MakeHexa2Faces(
                        geompy.MakeFaceHW(blobSeparatingDistance_, channel1MiddleWidth_, 1),
                        geompy.TranslateDXDYDZ(
                            geompy.MakeFaceHW(blobSeparatingDistance_, channel1MiddleWidth_ + 2.0*delta_, 1),
                            0, 0, depth_
                        )
                    ),
                    0,
                    blobRadius_ + 0.5*channel1MiddleWidth_,
                    -0.5*depth_
                ),
                geompy.TranslateDXDYDZ(
                    geompy.MakeHexa2Faces(
                        geompy.MakeFaceHW(blobSeparatingDistance_, channel2MiddleWidth_, 1),
                        geompy.TranslateDXDYDZ(
                            geompy.MakeFaceHW(blobSeparatingDistance_, channel2MiddleWidth_ + 2.0*delta_, 1),
                            0, 0, depth_
                        )
                    ),
                    0,
                    -(blobRadius_ + 0.5*channel2MiddleWidth_),
                    -0.5*depth_
                )
            ],
            True,
            True
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

#splittedChannelRightPart =\
    #geompy.TranslateDXDYDZ(
        #geompy.MakeCutList(
            #geompy.MakeFuseList(
                #[
                    #geompy.MakePrismVecH(
                        #geompy.MakeCutList(
                            #geompy.MakeFace(ellipseUp, 1),
                            #[
                                #geompy.MakeTranslation(
                                    #geompy.MakeFaceHW(2.0*a1_, 2.0*a1_, 1),
                                    #0, -a1_, 0
                                #),
                                #geompy.MakeTranslation(
                                    #geompy.MakeFaceHW(2.0*a1_, 2.0*a1_, 1),
                                    #-a1_, 0, 0
                                #)
                            #],
                            #True
                        #),
                        #OZ,
                        #depth_
                    #),
                    #geompy.MakePrismVecH(
                        #geompy.MakeCutList(
                            #geompy.MakeFace(ellipseDown, 1),
                            #[
                                #geompy.MakeTranslation(
                                    #geompy.MakeFaceHW(2.0*a2_, 2.0*a2_, 1),
                                    #0, a2_, 0
                                #),
                                #geompy.MakeTranslation(
                                    #geompy.MakeFaceHW(2.0*a2_, 2.0*a2_, 1),
                                    #-a2_, 0, 0
                                #)
                            #],
                            #True
                        #),
                        #OZ,
                        #depth_
                    #),
                    #geompy.TranslateDXDYDZ(
                        #geompy.MakeBoxDXDYDZ(
                            #max(a1_, a2_) + drainageLength_, drainageWidth_, depth_
                        #),
                        #0, -0.5*drainageWidth_, 0
                    #)
                #],
                #True, True
            #),
            #[geompy.MakeCylinderRH(blobRadius_, depth_)],
            #True
        #),
        #0.5*blobSeparatingDistance_, 0, -0.5*depth_
    #)

#splittedChannelLeftPart =\
    #geompy.MakeMirrorByPlane(
        #splittedChannelRightPart,
        #geompy.MakePlane(O, OX, 1000)
    #)

#splittedChannel =\
    #geompy.MakeFuseList(
        #[
            #splittedChannelRightPart,
            #splittedChannelLeftPart,
            #geompy.TranslateDXDYDZ(
                #geompy.MakeBoxDXDYDZ(blobSeparatingDistance_, channel1MiddleWidth_, depth_),
                #-0.5*blobSeparatingDistance_, blobRadius_, -0.5*depth_
            #),
            #geompy.TranslateDXDYDZ(
                #geompy.MakeBoxDXDYDZ(blobSeparatingDistance_, channel2MiddleWidth_, depth_),
                #-0.5*blobSeparatingDistance_, -blobRadius_ - channel2MiddleWidth_, -0.5*depth_
            #)
        #],
        #False, False
    #)

##Create groups for walls extraction
leftSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
rightSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
frontSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
backSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
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
geompy.UnionList(
    frontSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(bb[0], bb[2], bb[4]),
                geompy.MakeVertex(bb[1], bb[3], bb[5])
            ),
            0, 0, abs(bb[5] - bb[4])
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
geompy.UnionList(
    backSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(bb[0], bb[2], bb[4]),
                geompy.MakeVertex(bb[1], bb[3], bb[5])
            ),
            0, 0, -abs(bb[5] - bb[4])
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
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
geompy.DifferenceList(walls, geompy.SubShapeAll(frontSide, geompy.ShapeType["FACE"]))
geompy.DifferenceList(walls, geompy.SubShapeAll(backSide, geompy.ShapeType["FACE"]))

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
geompy.addToStudyInFather( splittedChannel, frontSide, 'frontSide' )
geompy.addToStudyInFather( splittedChannel, backSide, 'backSide' )
geompy.addToStudyInFather( splittedChannel, walls, 'walls' )


###
### SMESH component
###

smesh = smeshBuilder.New()
##smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 ## multiples meshes built in parallel, complex and numerous mesh edition (performance)

splittedChannelMesh = smesh.Mesh(splittedChannel)

splittedChannelMesh.Segment().Adaptive(
    0.02*min(
        drainageLength_,
        drainageWidth_,
        channel1MiddleWidth_,
        channel2MiddleWidth_,
        channel1InOutWidth_,
        channel2InOutWidth_,
        depth_,
        blobRadius_
    ),
    0.1*min(
        drainageLength_,
        drainageWidth_,
        channel1MiddleWidth_,
        channel2MiddleWidth_,
        channel1InOutWidth_,
        channel2InOutWidth_,
        depth_,
        blobRadius_
    ),
    2.0
)

NETGEN_2D_Parameters = splittedChannelMesh.Triangle(algo=smeshBuilder.NETGEN_2D).Parameters()
NETGEN_2D_Parameters.SetMaxSize(
    min(
        drainageLength_,
        drainageWidth_,
        channel1MiddleWidth_,
        channel2MiddleWidth_,
        channel1InOutWidth_,
        channel2InOutWidth_,
        depth_,
        blobRadius_
    )
)
NETGEN_2D_Parameters.SetMinSize(
    0.02*min(
        drainageLength_,
        drainageWidth_,
        channel1MiddleWidth_,
        channel2MiddleWidth_,
        channel1InOutWidth_,
        channel2InOutWidth_,
        depth_,
        blobRadius_
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

#NETGEN_3D = splittedChannelMesh.Tetrahedron()
#NETGEN_3D_Parameters = NETGEN_3D.Parameters()
#NETGEN_3D_Parameters.SetMaxSize(
    #min(
        #drainageLength_,
        #drainageWidth_,
        #channel1MiddleWidth_,
        #channel2MiddleWidth_,
        #channel1InOutWidth_,
        #channel2InOutWidth_,
        #depth_,
        #blobRadius_
    #)
#)
#NETGEN_3D_Parameters.SetMinSize(
    #0.02*min(
        #drainageLength_,
        #drainageWidth_,
        #channel1MiddleWidth_,
        #channel2MiddleWidth_,
        #channel1InOutWidth_,
        #channel2InOutWidth_,
        #depth_,
        #blobRadius_
    #)
#)
#NETGEN_3D_Parameters.SetOptimize( 1 )
#NETGEN_3D_Parameters.SetFineness( 5 )
#NETGEN_3D_Parameters.SetGrowthRate( 0.05 )
#NETGEN_3D_Parameters.SetElemSizeWeight( 0 )
#NETGEN_3D_Parameters.SetCheckOverlapping( 0 )
#NETGEN_3D_Parameters.SetCheckChartBoundary( 8 )
#Viscous_Layers = NETGEN_3D.ViscousLayers(
    #5, 2, 1,
    #[
        #geompy.GetSubShapeID(
            #splittedChannel,
            #leftSide
        #),
        #geompy.GetSubShapeID(
            #splittedChannel,
            #rightSide
        #)
    #],
    #1,
    #smeshBuilder.FACE_OFFSET
    ##smeshBuilder.NODE_OFFSET
#)

isDone = splittedChannelMesh.Compute()
splittedChannelWalls = splittedChannelMesh.GroupOnGeom(walls,'walls',SMESH.FACE)

#splittedChannelMesh.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(backSide,'backSide',SMESH.FACE)
#splittedChannelMesh.GroupOnGeom(walls,'walls',SMESH.FACE)





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


try:
  splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/splittedChannel.stl", 1 )
  splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/splittedChannelWalls.stl", 1, splittedChannelWalls)
  pass
except:
  print('ExportSTL() failed. Invalid file name?')

salome.myStudy.SaveAs(absoluteCasePath_ + "/splittedChannel.hdf", False, False)

#try:
  #splittedChannelMesh.ExportUNV( absoluteCasePath_ + "/splittedChannelMesh.unv" )
  #pass
#except:
  #print('ExportUNV() failed. Invalid file name?')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
