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
collectorLength_ = 2000.0
collectorRadius_ = 1000.0
drainageLength_ = 2000.0
channelLength_ = 4000.0
channelWidth_ = 200.0
channelDepth_ = 1000.0
NbStep_ = 4
step_ = (channelLength_ - channelWidth_)/NbStep_
####################################################
##             End of variables section           ##
####################################################

geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

drainageChannel = []

for index in range(NbStep_ + 1):
    drainageChannel.append(
        geompy.MakeFuseList(
            [
                geompy.MakeRotation(
                    geompy.MakeFuseList(
                        [
                            geompy.MakeTranslation(
                                geompy.MakeBoxDXDYDZ(0.5*channelLength_ - 0.5*channelWidth_, channelWidth_, channelDepth_),
                                0, -0.5*channelWidth_, 0
                            ),
                            geompy.MakeCylinderRH(0.5*channelWidth_, channelDepth_),
                            geompy.MakeTranslation(
                                geompy.MakeCylinderRH(0.5*channelWidth_, channelDepth_),
                                0.5*channelLength_ - 0.5*channelWidth_, 0, 0
                            ),
                        ],
                        True,
                        False
                    ),
                    OZ,
                    math.atan2(
                        (0.5*channelLength_ - 0.5*channelWidth_) - step_*index,
                        math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2)
                    )
                ),
                geompy.MakeTranslation(
                    geompy.MakeFuseList(
                        [
                            geompy.MakeTranslation(
                                geompy.MakeBoxDXDYDZ(
                                    drainageLength_ - math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2),
                                    channelWidth_,
                                    channelDepth_
                                ),
                                0, -0.5*channelWidth_, 0
                            ),
                            geompy.MakeCylinderRH(0.5*channelWidth_, channelDepth_),
                            geompy.MakeTranslation(
                                geompy.MakeCylinderRH(0.5*channelWidth_, channelDepth_),
                                drainageLength_ - math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2), 0, 0
                            ),
                        ],
                        True,
                        False
                    ),
                    math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2),
                    (0.5*channelLength_ - 0.5*channelWidth_) - step_*index,
                    0
                )
            ],
            True,
            False
        )
    )

    #geompy.addToStudy(drainageChannel[index], 'drainageChannel' + str(index))


collectorAndDrainageChannels =\
    geompy.MakeFuseList(
        [
            geompy.MakeFuseList(
                [
                    geompy.MakeTranslation(
                        geompy.MakeBoxDXDYDZ(collectorLength_, 2.0*collectorRadius_, channelDepth_),
                        -collectorLength_, -collectorRadius_, 0
                    ),
                    geompy.MakeCylinderRH(collectorRadius_, channelDepth_),
                    geompy.MakeTranslation(
                        geompy.MakeCylinderRH(collectorRadius_, channelDepth_),
                        -collectorLength_, 0, 0
                    ),
                ],
                True,
                False
            )
        ]
        + drainageChannel,
        True,
        False
    )

lattice =\
    geompy.MakeTranslation(
        geompy.MakeFuseList(
            [
                geompy.MakeMultiTranslation1D(geompy.MakeBoxDXDYDZ(channelLength_, channelWidth_, channelDepth_), OY, step_, NbStep_ + 1),
                geompy.MakeMultiTranslation1D(geompy.MakeBoxDXDYDZ(channelWidth_, channelLength_, channelDepth_), OX, step_, NbStep_ + 1),
                geompy.MakeTranslation(
                    collectorAndDrainageChannels,
                    -(drainageLength_ - math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2)),
                    0.5*channelLength_,
                    0
                ),
                geompy.MakeTranslation(
                    geompy.MakeRotation(
                        collectorAndDrainageChannels,
                        OZ,
                        0.5*math.pi
                    ),
                    0.5*channelLength_,
                    -(drainageLength_ - math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2)),
                    0
                ),
                geompy.MakeTranslation(
                    geompy.MakeRotation(
                        collectorAndDrainageChannels,
                        OZ,
                        math.pi
                    ),
                    channelLength_ + drainageLength_ - math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2),
                    0.5*channelLength_,
                    0
                ),
                geompy.MakeTranslation(
                    geompy.MakeRotation(
                        collectorAndDrainageChannels,
                        OZ,
                        1.5*math.pi
                    ),
                    0.5*channelLength_,
                    channelLength_ + drainageLength_ - math.sqrt((0.5*channelLength_ - 0.5*channelWidth_)**2 - ((0.5*channelLength_ - 0.5*channelWidth_) - step_*index)**2),
                    0
                )
                #geompy.MakeTranslation(
                    #geompy.MakeMultiTranslation1D(geompy.MakeBoxDXDYDZ(drainageLength_, channelWidth_, channelDepth_), OY, step_, NbStep_ + 1),
                    #-drainageLength_, 0, 0
                #),
                #geompy.MakeTranslation(
                    #geompy.MakeMultiTranslation1D(geompy.MakeBoxDXDYDZ(drainageLength_, channelWidth_, channelDepth_), OY, step_, NbStep_ + 1),
                    #channelLength_, 0, 0
                #),
                #geompy.MakeTranslation(
                    #geompy.MakeMultiTranslation1D(geompy.MakeBoxDXDYDZ(channelWidth_, drainageLength_, channelDepth_), OX, step_, NbStep_ + 1),
                    #0, -drainageLength_, 0
                #),
                #geompy.MakeTranslation(
                    #geompy.MakeMultiTranslation1D(geompy.MakeBoxDXDYDZ(channelWidth_, drainageLength_, channelDepth_), OX, step_, NbStep_ + 1),
                    #0, channelLength_, 0
                #)
            ],
            True,
            True
        ),
        -0.5*channelLength_, -0.5*channelLength_, -0.5*channelDepth_
        #0, 0, 0
    )

#Create fillet
#front = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
#opposite = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])

#geompy.UnionList(
    #front,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_ + 2.0*drainageLength_, channelLength_ + 2.0*drainageLength_, channelDepth_),
            #-0.5*channelLength_ - drainageLength_, -0.5*channelLength_ - drainageLength_, -0.5*channelDepth_ + channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
#geompy.UnionList(
    #opposite,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_ + 2.0*drainageLength_, channelLength_ + 2.0*drainageLength_, channelDepth_),
            #-0.5*channelLength_ - drainageLength_, -0.5*channelLength_ - drainageLength_, -0.5*channelDepth_ - channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)

#filletEdges = geompy.CreateGroup(lattice, geompy.ShapeType["EDGE"])

#geompy.UnionList(
    #filletEdges,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_, channelLength_, channelDepth_),
            #-0.5*channelLength_, -0.5*channelLength_, -0.5*channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["EDGE"],
        #GEOM.ST_ONIN
    #)
#)

#geompy.DifferenceList(filletEdges, geompy.SubShapeAll(front, geompy.ShapeType["EDGE"]))
#geompy.DifferenceList(filletEdges, geompy.SubShapeAll(opposite, geompy.ShapeType["EDGE"]))

#lattice =\
    #geompy.MakeFillet(
        #lattice,
        #0.1*channelWidth_,
        #geompy.ShapeType["EDGE"],
        #[geompy.GetSubShapeID(lattice, edgeI) for edgeI in geompy.SubShapeAll(filletEdges, geompy.ShapeType["EDGE"])]
    #)

#Create groups for walls extraction
#leftSide = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
#bottomSide = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
#rightSide = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
#topSide = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
frontSide = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
oppositeSide = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])
walls = geompy.CreateGroup(lattice, geompy.ShapeType["FACE"])

bb = geompy.BoundingBox(lattice)

#bbBox =\
    #geompy.MakeBoxTwoPnt(
        #geompy.MakeVertex(bb[0], bb[2], bb[4]),
        #geompy.MakeVertex(bb[1], bb[3], bb[5])
    #)

#geompy.addToStudy( bbBox, 'bbBox' )

#geompy.UnionList(
    #leftSide,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_ + 2.0*drainageLength_, channelLength_ + 2.0*drainageLength_, channelDepth_),
            #-0.5*channelLength_ - drainageLength_ - channelLength_ - 2.0*drainageLength_, -0.5*channelLength_ - drainageLength_, -0.5*channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
#geompy.UnionList(
    #bottomSide,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_ + 2.0*drainageLength_, channelLength_ + 2.0*drainageLength_, channelDepth_),
            #-0.5*channelLength_ - drainageLength_, -0.5*channelLength_ - drainageLength_ - channelLength_ - 2.0*drainageLength_, -0.5*channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
#geompy.UnionList(
    #rightSide,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_ + 2.0*drainageLength_, channelLength_ + 2.0*drainageLength_, channelDepth_),
            #-0.5*channelLength_ - drainageLength_ + channelLength_ + 2.0*drainageLength_, -0.5*channelLength_ - drainageLength_, -0.5*channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
#geompy.UnionList(
    #topSide,
    #geompy.GetShapesOnBox(
        #geompy.TranslateDXDYDZ(
            #geompy.MakeBoxDXDYDZ(channelLength_ + 2.0*drainageLength_, channelLength_ + 2.0*drainageLength_, channelDepth_),
            #-0.5*channelLength_ - drainageLength_, -0.5*channelLength_ - drainageLength_ + channelLength_ + 2.0*drainageLength_, -0.5*channelDepth_
        #),
        #lattice,
        #geompy.ShapeType["FACE"],
        #GEOM.ST_ON
    #)
#)
geompy.UnionList(
    frontSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(bb[0], bb[2], bb[4]),
                geompy.MakeVertex(bb[1], bb[3], bb[5])
            ),
            0, 0, channelDepth_
        ),
        lattice,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
geompy.UnionList(
    oppositeSide,
    geompy.GetShapesOnBox(
        geompy.TranslateDXDYDZ(
            geompy.MakeBoxTwoPnt(
                geompy.MakeVertex(bb[0], bb[2], bb[4]),
                geompy.MakeVertex(bb[1], bb[3], bb[5])
            ),
            0, 0, -channelDepth_
        ),
        lattice,
        geompy.ShapeType["FACE"],
        GEOM.ST_ON
    )
)
# Walls
geompy.UnionList(
    walls,
    geompy.SubShapeAll(
        lattice,
        geompy.ShapeType["FACE"]
    )
)
#geompy.DifferenceList(walls, geompy.SubShapeAll(leftSide, geompy.ShapeType["FACE"]))
#geompy.DifferenceList(walls, geompy.SubShapeAll(bottomSide, geompy.ShapeType["FACE"]))
#geompy.DifferenceList(walls, geompy.SubShapeAll(rightSide, geompy.ShapeType["FACE"]))
#geompy.DifferenceList(walls, geompy.SubShapeAll(topSide, geompy.ShapeType["FACE"]))
geompy.DifferenceList(walls, geompy.SubShapeAll(frontSide, geompy.ShapeType["FACE"]))
geompy.DifferenceList(walls, geompy.SubShapeAll(oppositeSide, geompy.ShapeType["FACE"]))

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( collectorAndDrainageChannels, 'collectorAndDrainageChannels' )
geompy.addToStudy( lattice, 'lattice' )
#geompy.addToStudyInFather( lattice, filletEdges, 'filletEdges' )
#geompy.addToStudyInFather( lattice, leftSide, 'leftSide' )
#geompy.addToStudyInFather( lattice, bottomSide, 'bottomSide' )
#geompy.addToStudyInFather( lattice, rightSide, 'rightSide' )
#geompy.addToStudyInFather( lattice, topSide, 'topSide' )
geompy.addToStudyInFather( lattice, frontSide, 'frontSide' )
geompy.addToStudyInFather( lattice, oppositeSide, 'oppositeSide' )
geompy.addToStudyInFather( lattice, walls, 'walls' )


##
## SMESH component
##

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

latticeMesh = smesh.Mesh(lattice)
NETGEN_1D_2D = latticeMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize(
    min(
        channelLength_,
        channelWidth_,
        channelDepth_
    )
)
NETGEN_2D_Parameters_1.SetMinSize(
    0.02*min(
        channelLength_,
        channelWidth_,
        channelDepth_
    )
)
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 5 )
NETGEN_2D_Parameters_1.SetGrowthRate( 0.5 )
NETGEN_2D_Parameters_1.SetNbSegPerEdge( 3 )
NETGEN_2D_Parameters_1.SetNbSegPerRadius( 5 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 21874 )
NETGEN_2D_Parameters_1.SetUseDelauney( 0 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 184 )
isDone = latticeMesh.Compute()
latticeWalls = latticeMesh.GroupOnGeom(walls, 'walls', SMESH.FACE)

try:
  latticeMesh.ExportSTL( absoluteCasePath_ + "/lattice.stl", 1 )
  latticeMesh.ExportSTL( absoluteCasePath_ + "/latticeWalls.stl", 1, latticeWalls)
  pass
except:
  print('ExportSTL() failed. Invalid file name?')

salome.myStudy.SaveAs(absoluteCasePath_ + "/lattice.hdf", False, False)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
