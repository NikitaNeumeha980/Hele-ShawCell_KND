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
drainageWidth_ = 300
drainageOffsetDY_ = 0
channel1Width_ = 310
channel2Width_ = 270
depth_ = 1000
#blobCentreOffsetDX_ = 0
#blobCentreOffsetDY_ = 0.050
blobSeparatingDistance_ = 0
blobRadius_ = 0.5*170
#channel1Radius_ = 200
#channel2Radius_ = 200
####################################################
##             End of variables section           ##
####################################################


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

splittedChannelRightPart =\
    geompy.TranslateDXDYDZ(
        geompy.MakeCutList(
            geompy.MakeFuseList(
                [
                    geompy.MakeCylinderRHA(blobRadius_ + channel1Width_, depth_, 90*math.pi/180.0),
                    geompy.Rotate(geompy.MakeCylinderRHA(blobRadius_ + channel2Width_, depth_, 90*math.pi/180.0), OZ, -90*math.pi/180.0),
                    geompy.TranslateDXDYDZ(
                        geompy.MakeBoxDXDYDZ(
                            max(blobRadius_ + channel1Width_, blobRadius_ + channel2Width_) + drainageLength_,
                            drainageWidth_,
                            depth_
                        ),
                        0,
                        -0.5*drainageWidth_ + drainageOffsetDY_,
                        0
                    )
                ],
                True,
                True
            ),
            [geompy.MakeCylinderRH(blobRadius_, depth_)],
            True
        ),
        0.5*blobSeparatingDistance_,
        0,
        -0.5*depth_
    )

splittedChannelLeftPart =\
    geompy.MakeMirrorByPlane(
        splittedChannelRightPart,
        geompy.MakePlane(O, OX, 2000)
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
                    geompy.MakeBoxDXDYDZ(
                        blobSeparatingDistance_,
                        channel1Width_,
                        depth_
                    ),
                    -0.5*blobSeparatingDistance_,
                    blobRadius_,
                    -0.5*depth_
                ),
                geompy.TranslateDXDYDZ(
                    geompy.MakeBoxDXDYDZ(
                        blobSeparatingDistance_,
                        channel2Width_,
                        depth_
                    ),
                    -0.5*blobSeparatingDistance_,
                    -blobRadius_ - channel2Width_,
                    -0.5*depth_
                )
            ],
            False,
            False
        )

#Create groups for walls extraction
domainBB = geompy.BoundingBox(splittedChannel)

leftSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
rightSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
frontSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
oppositeSide = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])
walls = geompy.CreateGroup(splittedChannel, geompy.ShapeType["FACE"])

geompy.UnionList(
    leftSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(float(domainBB[0]), float(domainBB[2]), float(domainBB[4])),
                geompy.MakeVertex(float(domainBB[1]), float(domainBB[3]), float(domainBB[5]))
            ),
            -float(abs(domainBB[1] - domainBB[0])), 0, 0
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
geompy.UnionList(
    rightSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(float(domainBB[0]), float(domainBB[2]), float(domainBB[4])),
                geompy.MakeVertex(float(domainBB[1]), float(domainBB[3]), float(domainBB[5]))
            ),
            float(abs(domainBB[1] - domainBB[0])), 0, 0
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
                geompy.MakeVertex(float(domainBB[0]), float(domainBB[2]), float(domainBB[4])),
                geompy.MakeVertex(float(domainBB[1]), float(domainBB[3]), float(domainBB[5]))
            ),
            0, 0, float(abs(domainBB[5] - domainBB[4]))
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
geompy.UnionList(
    oppositeSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(float(domainBB[0]), float(domainBB[2]), float(domainBB[4])),
                geompy.MakeVertex(float(domainBB[1]), float(domainBB[3]), float(domainBB[5]))
            ),
            0, 0, -float(abs(domainBB[5] - domainBB[4]))
        ),
        splittedChannel,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
# Walls
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
geompy.DifferenceList(walls, geompy.SubShapeAll(oppositeSide, geompy.ShapeType["FACE"]))

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( splittedChannelRightPart, 'splittedChannelRightPart' )
geompy.addToStudy( splittedChannelLeftPart, 'splittedChannelLeftPart' )
geompy.addToStudy( splittedChannel, 'splittedChannel' )
geompy.addToStudyInFather( splittedChannel, leftSide, 'leftSide' )
geompy.addToStudyInFather( splittedChannel, rightSide, 'rightSide' )
geompy.addToStudyInFather( splittedChannel, frontSide, 'frontSide' )
geompy.addToStudyInFather( splittedChannel, oppositeSide, 'oppositeSide' )
geompy.addToStudyInFather( splittedChannel, walls, 'walls' )


##
## SMESH component
##

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

splittedChannelMesh = smesh.Mesh(splittedChannel)
NETGEN_1D_2D = splittedChannelMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize(
    min(
        drainageLength_,
        drainageWidth_,
        channel1Width_,
        channel2Width_,
        depth_,
        blobRadius_
    )
)
NETGEN_2D_Parameters_1.SetMinSize(
    0.02*min(
        drainageLength_,
        drainageWidth_,
        channel1Width_,
        channel2Width_,
        depth_,
        blobRadius_
    )
)
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 5 )
NETGEN_2D_Parameters_1.SetGrowthRate( 0.5 )
NETGEN_2D_Parameters_1.SetNbSegPerEdge( 3 )
NETGEN_2D_Parameters_1.SetNbSegPerRadius( 10 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 21874 )
NETGEN_2D_Parameters_1.SetUseDelauney( 0 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 184 )
isDone = splittedChannelMesh.Compute()
splittedChannelWalls = splittedChannelMesh.GroupOnGeom(walls,'walls',SMESH.FACE)

try:
  splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/splittedChannel.stl", 1 )
  splittedChannelMesh.ExportSTL( absoluteCasePath_ + "/splittedChannelWalls.stl", 1, splittedChannelWalls)
  pass
except:
  print('ExportSTL() failed. Invalid file name?')

salome.myStudy.SaveAs(absoluteCasePath_ + "/splittedChannel.hdf", False, False)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
