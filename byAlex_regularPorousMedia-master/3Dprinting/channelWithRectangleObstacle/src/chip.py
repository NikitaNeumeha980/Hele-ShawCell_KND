## to do list

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
absoluteCasePath_ = ""

printerResolutionXY_ = [2560, 1440]
voxelSizeDXDYDZ_ =\
    [
        Decimal('0.04725'),
        Decimal('0.04725'),
        Decimal('0.01')
    ]

chipDXDYDZ_ =\
    [
        voxelSizeDXDYDZ_[0]*(Decimal('10.0')/voxelSizeDXDYDZ_[0]).quantize(Decimal('1'), ROUND_CEILING),
        voxelSizeDXDYDZ_[1]*(Decimal('40.0')/voxelSizeDXDYDZ_[1]).quantize(Decimal('1'), ROUND_CEILING),
        voxelSizeDXDYDZ_[2]*(Decimal('2.0')/voxelSizeDXDYDZ_[2]).quantize(Decimal('1'), ROUND_CEILING)
    ]






phase12ExpandingChannelDZ_ = voxelSizeDXDYDZ_[2]*math.ceil(0.2/voxelSizeDXDYDZ_[2])

receiverRDZ_ =\
    [
        [
            2.5,
            voxelSizeDXDYDZ_[2]*math.ceil(1.0/voxelSizeDXDYDZ_[2])
        ],
        [
            2.5,
            voxelSizeDXDYDZ_[2]*math.ceil(1.0/voxelSizeDXDYDZ_[2])
        ],
        [
            2.5,
            voxelSizeDXDYDZ_[2]*math.ceil(1.0/voxelSizeDXDYDZ_[2])
        ]
    ]

channelWDZ_ =\
    [
        [
            2.0*voxelSizeDXDYDZ_[0],
            voxelSizeDXDYDZ_[2]*math.ceil(0.1/voxelSizeDXDYDZ_[2])
        ],
        [
            2.0*voxelSizeDXDYDZ_[0],
            voxelSizeDXDYDZ_[2]*math.ceil(0.1/voxelSizeDXDYDZ_[2])
        ],
        [
            2.0*voxelSizeDXDYDZ_[0],
            voxelSizeDXDYDZ_[2]*math.ceil(0.1/voxelSizeDXDYDZ_[2])
        ]
    ]

supportsNbTimes_ = [6, 21]
#supportsRadius_ = 0.25

supportsDXDYDZ_ =\
    [
        voxelSizeDXDYDZ_[0]*math.ceil(0.25/voxelSizeDXDYDZ_[0]),
        voxelSizeDXDYDZ_[1]*math.ceil(0.25/voxelSizeDXDYDZ_[1]),
        voxelSizeDXDYDZ_[2]*math.ceil(5.0/voxelSizeDXDYDZ_[2])
    ]

baseDXDYDZ_ =\
    [
        chipDXDYDZ_[0],
        chipDXDYDZ_[1],
        voxelSizeDXDYDZ_[2]*math.ceil(1.0/voxelSizeDXDYDZ_[2])
    ]

coverDXDYDZ_ =\
    [
        chipDXDYDZ_[0],
        chipDXDYDZ_[1],
        chipDXDYDZ_[2]
    ]

plugRrDZ_ =\
    [
        [
            1.5,
            0.75,
            supportsDXDYDZ_[2]
        ],
        [
            1.5,
            0.75,
            supportsDXDYDZ_[2]
        ],
        [
            1.5,
            0.75,
            supportsDXDYDZ_[2]
        ]
    ]

#tolerance_ = 1e-07
####################################################
##             End of variables section           ##
####################################################


####################################################
##                 GEOM component                 ##
####################################################

geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

receiverCentre =\
    [
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0.5*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0.1*chipDXDYDZ_[1]/voxelSizeDXDYDZ_[1]),
            0
        ),
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0.5*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0.25*chipDXDYDZ_[1]/voxelSizeDXDYDZ_[1]),
            0
        ),
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0.5*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0.9*chipDXDYDZ_[1]/voxelSizeDXDYDZ_[1]),
            0
        )
    ]

# Points positions
# for channel path
#        |
#        |
#        |
#  4---3(IC)---2
#  |     |     |
#  |     |     |
#  |           |
#  5-----0-----1

intersectionCentre =\
    [
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0.3*chipDXDYDZ_[1]/voxelSizeDXDYDZ_[1]),
            voxelSizeDXDYDZ_[2]*math.ceil(0/voxelSizeDXDYDZ_[2])
        )
    ]

phase1ChannelMiddleLineContourPoints =\
    [
        #0
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0/voxelSizeDXDYDZ_[1]),
            voxelSizeDXDYDZ_[2]*math.ceil(0/voxelSizeDXDYDZ_[2])
        ),
        #1
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0.4*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0/voxelSizeDXDYDZ_[1]),
            voxelSizeDXDYDZ_[2]*math.ceil(0/voxelSizeDXDYDZ_[2])
        ),
        #2
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(0.4*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0.3*chipDXDYDZ_[1]/voxelSizeDXDYDZ_[1]),
            voxelSizeDXDYDZ_[2]*math.ceil(0/voxelSizeDXDYDZ_[2])
        ),
        #3
        intersectionCentre[0],
        #4
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(-0.4*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0.3*chipDXDYDZ_[1]/voxelSizeDXDYDZ_[1]),
            voxelSizeDXDYDZ_[2]*math.ceil(0/voxelSizeDXDYDZ_[2])
        ),
        #5
        geompy.MakeVertex(
            voxelSizeDXDYDZ_[0]*math.ceil(-0.4*chipDXDYDZ_[0]/voxelSizeDXDYDZ_[0]),
            voxelSizeDXDYDZ_[1]*math.ceil(0/voxelSizeDXDYDZ_[1]),
            voxelSizeDXDYDZ_[2]*math.ceil(0/voxelSizeDXDYDZ_[2])
        ),
    ]

# Feed channels
sk = geompy.Sketcher2D()
sk.addPoint(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[0])[0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[0])[1] - 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[1])[0] + 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[1])[1] - 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[2])[0] + 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[2])[1] + 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[3])[0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[3])[1] + 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[4])[0] - 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[4])[1] + 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[5])[0] - 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[5])[1] - 0.5*channelWDZ_[0][0]
)
sk.close()

phase1ChannelOutsideContour =\
    sk.wire(
        geompy.MakeMarker(
            geompy.PointCoordinates(receiverCentre[0])[0], geompy.PointCoordinates(receiverCentre[0])[1], chipDXDYDZ_[2] - channelWDZ_[0][1],
            1, 0, 0,
            0, 1, 0
        )
    )

sk = geompy.Sketcher2D()
sk.addPoint(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[0])[0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[0])[1] + 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[1])[0] - 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[1])[1] + 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[2])[0] - 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[2])[1] - 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[3])[0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[3])[1] - 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[4])[0] + 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[4])[1] - 0.5*channelWDZ_[0][0]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[5])[0] + 0.5*channelWDZ_[0][0],
    geompy.PointCoordinates(phase1ChannelMiddleLineContourPoints[5])[1] + 0.5*channelWDZ_[0][0]
)
sk.close()

phase1ChannelInsideContour =\
    sk.wire(
        geompy.MakeMarker(
            geompy.PointCoordinates(receiverCentre[0])[0], geompy.PointCoordinates(receiverCentre[0])[1], chipDXDYDZ_[2] - channelWDZ_[0][1],
            1, 0, 0,
            0, 1, 0
        )
    )

sk = geompy.Sketcher2D()
sk.addPoint(
    0,
    0
)
sk.addSegmentAbsolute(
    0.5*channelWDZ_[1][0],
    0
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(intersectionCentre[0])[0] + 0.5*channelWDZ_[1][0],
    geompy.PointCoordinates(intersectionCentre[0])[1]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(intersectionCentre[0])[0] - 0.5*channelWDZ_[1][0],
    geompy.PointCoordinates(intersectionCentre[0])[1]
)
sk.addSegmentAbsolute(
    0.5*channelWDZ_[1][0],
    0
)
sk.close()

phase2ChannelContour =\
    sk.wire(
        geompy.MakeMarker(
            geompy.PointCoordinates(receiverCentre[1])[0], geompy.PointCoordinates(receiverCentre[1])[1], chipDXDYDZ_[2] - channelWDZ_[1][1],
            1, 0, 0,
            0, 1, 0
        )
    )

# Drainage channel
sk = geompy.Sketcher2D()
sk.addPoint(
    geompy.PointCoordinates(receiverCentre[1])[0] + geompy.PointCoordinates(intersectionCentre[0])[0],
    geompy.PointCoordinates(receiverCentre[1])[1] + geompy.PointCoordinates(intersectionCentre[0])[1]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(receiverCentre[1])[0] + geompy.PointCoordinates(intersectionCentre[0])[0] + 0.5*channelWDZ_[2][0],
    geompy.PointCoordinates(receiverCentre[1])[1] + geompy.PointCoordinates(intersectionCentre[0])[1]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(receiverCentre[2])[0] + 0.5*channelWDZ_[2][0],
    geompy.PointCoordinates(receiverCentre[2])[1]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(receiverCentre[2])[0] - 0.5*channelWDZ_[2][0],
    geompy.PointCoordinates(receiverCentre[2])[1]
)
sk.addSegmentAbsolute(
    geompy.PointCoordinates(receiverCentre[1])[0] + geompy.PointCoordinates(intersectionCentre[0])[0] - 0.5*channelWDZ_[2][0],
    geompy.PointCoordinates(receiverCentre[1])[1] + geompy.PointCoordinates(intersectionCentre[0])[1]
)
sk.close()

phase12ChannelContour =\
    sk.wire(
        geompy.MakeMarker(
            0, 0, chipDXDYDZ_[2] - channelWDZ_[2][1],
            1, 0, 0,
            0, 1, 0
        )
    )

chip =\
    geompy.MakeCutList(
        # Shape
        geompy.MakeBoxDXDYDZ(
            chipDXDYDZ_[0],
            chipDXDYDZ_[1],
            chipDXDYDZ_[2]
        ),
        [
            # Receivers
            geompy.TranslateDXDYDZ(geompy.MakeCylinder(receiverCentre[0], OZ, receiverRDZ_[0][0], receiverRDZ_[0][1]), 0, 0, chipDXDYDZ_[2] - receiverRDZ_[0][1]),
            geompy.TranslateDXDYDZ(geompy.MakeCylinder(receiverCentre[1], OZ, receiverRDZ_[1][0], receiverRDZ_[1][1]), 0, 0, chipDXDYDZ_[2] - receiverRDZ_[1][1]),
            geompy.TranslateDXDYDZ(geompy.MakeCylinder(receiverCentre[2], OZ, receiverRDZ_[2][0], receiverRDZ_[2][1]), 0, 0, chipDXDYDZ_[2] - receiverRDZ_[2][1]),
            # Feed channels
            geompy.MakePrismVecH(
                geompy.MakeFaceWires([phase1ChannelOutsideContour, phase1ChannelInsideContour], 1),
                OZ,
                channelWDZ_[0][1]
            ),
            geompy.MakePrismVecH(
                geompy.MakeFaceWires([phase2ChannelContour], 1),
                OZ,
                channelWDZ_[1][1]
            ),
            # Drainage channel
            geompy.MakePrismVecH(
                geompy.MakeFaceWires([phase12ChannelContour], 1),
                OZ,
                channelWDZ_[2][1]
            ),
            # Drainage channel expansion
            #geompy.MakePrismVecH(
                #geompy.MakeFuseList(
                    #[
                        #geompy.MakeDiskPntVecR(
                            #geompy.MakeVertex(
                                #geompy.PointCoordinates(receiverCentre[0])[0],
                                #geompy.PointCoordinates(receiverCentre[0])[1],
                                #geompy.PointCoordinates(receiverCentre[0])[3]
                            #),
                            #OZ,
                            #0.05*chipDXDYDZ_[0]
                        #),
                        #sphere, box
                    #]
                #),
                #OZ,
                #phase12ExpandingChannelDZ_
            #)
            # Drainage channel expansion (ellipse)
            #geompy.MakePrismVecH(
                #geompy.MakeFaceWires(
                    #[
                        #geompy.MakeEllipse(
                            #geompy.MakeVertex(
                                #geompy.PointCoordinates(receiverCentre[0])[0],
                                #0.5*(
                                    #geompy.PointCoordinates(receiverCentre[2])[1]
                                    #- receiverRDZ_[2][0]
                                    #- (geompy.PointCoordinates(receiverCentre[0])[1] + 0.3*chipDXDYDZ_[1] + 0.5*channelWDZ_[0][0])
                                #)
                                #+ geompy.PointCoordinates(receiverCentre[0])[1] + 0.3*chipDXDYDZ_[1]
                                #+ 0.5*channelWDZ_[0][0],
                                #chipDXDYDZ_[2] - phase12ExpandingChannelDZ_
                            #),
                            #OZ,
                            #0.5*(geompy.PointCoordinates(receiverCentre[2])[1] - receiverRDZ_[2][0]
                            #- (geompy.PointCoordinates(receiverCentre[0])[1] + 0.3*chipDXDYDZ_[1] + 0.5*channelWDZ_[0][0])
                            #- 2.0*channelWDZ_[0][0]),
                            #0.05*chipDXDYDZ_[0],
                            #OY
                        #)
                    #],
                    #1
                #),
                #OZ,
                #phase12ExpandingChannelDZ_
            #)
        ],
        True
    )

geompy.TranslateDXDYDZ(chip, 0, 0, baseDXDYDZ_[2] + supportsDXDYDZ_[2])

supportsForChip =\
    geompy.MakeMultiTranslation2D(
        #geompy.MakeCylinder(
            #O,
            #OZ,
            #supportsRadius_,
            #supportsDXDYDZ_[2]
        #),
        geompy.MakeBoxDXDYDZ(
            supportsDXDYDZ_[0],
            supportsDXDYDZ_[1],
            supportsDXDYDZ_[2]
        ),
        OX, (chipDXDYDZ_[0] - supportsDXDYDZ_[0])/(supportsNbTimes_[0] - 1 + 1e-6), supportsNbTimes_[0],
        OY, (chipDXDYDZ_[1] - supportsDXDYDZ_[1])/(supportsNbTimes_[1] - 1 + 1e-6), supportsNbTimes_[1]
    )

geompy.TranslateDXDYDZ(supportsForChip, 0, 0, baseDXDYDZ_[2])

base =\
    geompy.MakeBoxDXDYDZ(
        baseDXDYDZ_[0],
        baseDXDYDZ_[1],
        baseDXDYDZ_[2]
    )

chipSupportsBase = geompy.MakeCompound([chip, supportsForChip, base])
#chipSupportsBase = geompy.MakeCommonList([chip, supportsForChip, base])

cover =\
    geompy.MakeFuseList(
        [
            geompy.MakeBoxDXDYDZ(
                coverDXDYDZ_[0],
                coverDXDYDZ_[1],
                coverDXDYDZ_[2]
            ),
            #geompy.MakeCylinder(receiverCentre[0], OZ, plugRrDZ_[0][0], plugRrDZ_[0][2]),
            #geompy.MakeCylinder(receiverCentre[1], OZ, plugRrDZ_[1][0], plugRrDZ_[1][2]),
            #geompy.MakeCylinder(receiverCentre[2], OZ, plugRrDZ_[2][0], plugRrDZ_[2][2]),
            geompy.MakeCylinder(receiverCentre[0], OZ, plugRrDZ_[0][0], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
            geompy.MakeCylinder(receiverCentre[1], OZ, plugRrDZ_[1][0], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
            geompy.MakeCylinder(receiverCentre[2], OZ, plugRrDZ_[2][0], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
            geompy.MakeCone(receiverCentre[2], OZ, receiverRDZ_[2][0] + (plugRrDZ_[2][0] - plugRrDZ_[2][1]), plugRrDZ_[2][0], 0.2*plugRrDZ_[2][2])
        ]
    )

cover =\
    geompy.MakeCutList(
        cover,
        [
            #geompy.MakeCylinder(receiverCentre[0], OZ, plugRrDZ_[0][1], plugRrDZ_[0][2]),
            #geompy.MakeCylinder(receiverCentre[1], OZ, plugRrDZ_[1][1], plugRrDZ_[1][2]),
            #geompy.MakeCylinder(receiverCentre[2], OZ, plugRrDZ_[2][1], plugRrDZ_[2][2]),
            geompy.MakeCylinder(receiverCentre[0], OZ, plugRrDZ_[0][1], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
            geompy.MakeCylinder(receiverCentre[1], OZ, plugRrDZ_[1][1], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
            geompy.MakeCylinder(receiverCentre[2], OZ, plugRrDZ_[2][1], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
            geompy.MakeCone(receiverCentre[2], OZ, receiverRDZ_[2][0], plugRrDZ_[2][1], 0.2*plugRrDZ_[2][2])
        ],
        True
    )

geompy.Rotate(
    cover,
    geompy.MakeVector(
        geompy.MakeVertex(0.5*chipDXDYDZ_[0], 0, 0.5*chipDXDYDZ_[2]),
        geompy.MakeVertex(0.5*chipDXDYDZ_[0], chipDXDYDZ_[1], 0.5*chipDXDYDZ_[2])
    ),
    math.pi
)

geompy.TranslateDXDYDZ(cover, 0, 0, baseDXDYDZ_[2] + supportsDXDYDZ_[2])

#supports = []

#for i in range(supportsNbTimes_[0]):
    #for j in range(supportsNbTimes_[1]):
        #if (
            #math.sqrt(
                #(chipDXDYDZ_[0]*i/(supportsNbTimes_[0] - 1 + 1e-6) - geompy.PointCoordinates(receiverCentre[0])[0])**2
                #+ (chipDXDYDZ_[1]*j/(supportsNbTimes_[1] - 1 + 1e-6) - geompy.PointCoordinates(receiverCentre[0])[1])**2
            #) >= plugRrDZ_[0][0] + supportsRadius_
            #and
            #math.sqrt(
                #(chipDXDYDZ_[0]*i/(supportsNbTimes_[0] - 1 + 1e-6) - geompy.PointCoordinates(receiverCentre[1])[0])**2
                #+ (chipDXDYDZ_[1]*j/(supportsNbTimes_[1] - 1 + 1e-6) - geompy.PointCoordinates(receiverCentre[1])[1])**2
            #) >= plugRrDZ_[1][0] + supportsRadius_
            #and
            #math.sqrt(
                #(chipDXDYDZ_[0]*i/(supportsNbTimes_[0] - 1 + 1e-6) - geompy.PointCoordinates(receiverCentre[2])[0])**2
                #+ (chipDXDYDZ_[1]*j/(supportsNbTimes_[1] - 1 + 1e-6) - geompy.PointCoordinates(receiverCentre[2])[1])**2
            #) >= plugRrDZ_[2][0] + supportsRadius_
        #):
            #supports.append(
                #geompy.MakeCylinder(
                    #geompy.MakeVertex(
                        #chipDXDYDZ_[0]*i/(supportsNbTimes_[0] - 1 + 1e-6),
                        #chipDXDYDZ_[1]*j/(supportsNbTimes_[1] - 1 + 1e-6),
                        #0
                    #),
                    #OZ,
                    #supportsRadius_,
                    #supportsDXDYDZ_[2]
                #)
            #)

#for i in range(supportsNbTimes_[0]):
    #for j in range(supportsNbTimes_[1]):
        #if (
            #math.sqrt(
                #((chipDXDYDZ_[0] - supportsDXDYDZ_[0])*i/(supportsNbTimes_[0] - 1 + 1e-6) + 0.5*supportsDXDYDZ_[0] - geompy.PointCoordinates(receiverCentre[0])[0])**2
                #+ ((chipDXDYDZ_[1] - supportsDXDYDZ_[1])*j/(supportsNbTimes_[1] - 1 + 1e-6) + 0.5*supportsDXDYDZ_[1] - geompy.PointCoordinates(receiverCentre[0])[1])**2
            #)
            #- math.sqrt(
                #(supportsDXDYDZ_[0])**2
                #+ (supportsDXDYDZ_[1])**2
            #) >= plugRrDZ_[0][0]
            #and
            #math.sqrt(
                #((chipDXDYDZ_[0] - supportsDXDYDZ_[0])*i/(supportsNbTimes_[0] - 1 + 1e-6) + 0.5*supportsDXDYDZ_[0] - geompy.PointCoordinates(receiverCentre[1])[0])**2
                #+ ((chipDXDYDZ_[1] - supportsDXDYDZ_[1])*j/(supportsNbTimes_[1] - 1 + 1e-6) + 0.5*supportsDXDYDZ_[1] - geompy.PointCoordinates(receiverCentre[1])[1])**2
            #)
            #- math.sqrt(
                #(supportsDXDYDZ_[0])**2
                #+ (supportsDXDYDZ_[1])**2
            #) >= plugRrDZ_[1][0]
            #and
            #math.sqrt(
                #((chipDXDYDZ_[0] - supportsDXDYDZ_[0])*i/(supportsNbTimes_[0] - 1 + 1e-6) + 0.5*supportsDXDYDZ_[0] - geompy.PointCoordinates(receiverCentre[2])[0])**2
                #+ ((chipDXDYDZ_[1] - supportsDXDYDZ_[1])*j/(supportsNbTimes_[1] - 1 + 1e-6) + 0.5*supportsDXDYDZ_[1] - geompy.PointCoordinates(receiverCentre[2])[1])**2
            #)
            #- math.sqrt(
                #(supportsDXDYDZ_[0])**2
                #+ (supportsDXDYDZ_[1])**2
            #) >= plugRrDZ_[2][0]
        #):
            #supports.append(
                #geompy.MakeBox(
                    #(chipDXDYDZ_[0] - supportsDXDYDZ_[0])*i/(supportsNbTimes_[0] - 1 + 1e-6),
                    #(chipDXDYDZ_[1] - supportsDXDYDZ_[1])*j/(supportsNbTimes_[1] - 1 + 1e-6),
                    #0,
                    #(chipDXDYDZ_[0] - supportsDXDYDZ_[0])*i/(supportsNbTimes_[0] - 1 + 1e-6) + supportsDXDYDZ_[0],
                    #(chipDXDYDZ_[1] - supportsDXDYDZ_[1])*j/(supportsNbTimes_[1] - 1 + 1e-6) + supportsDXDYDZ_[1],
                    #supportsDXDYDZ_[2]
                #)
            #)

#supportsForCover = geompy.MakeCompound(supports)

#geompy.TranslateDXDYDZ(supportsForCover, 0, 0, baseDXDYDZ_[2])

supportsForCover =\
    geompy.MakeMultiTranslation2D(
        geompy.MakeBoxDXDYDZ(
            supportsDXDYDZ_[0],
            supportsDXDYDZ_[1],
            supportsDXDYDZ_[2]
        ),
        OX, (chipDXDYDZ_[0] - supportsDXDYDZ_[0])/(supportsNbTimes_[0] - 1 + 1e-6), supportsNbTimes_[0],
        OY, (chipDXDYDZ_[1] - supportsDXDYDZ_[1])/(supportsNbTimes_[1] - 1 + 1e-6), supportsNbTimes_[1]
    )

geompy.TranslateDXDYDZ(supportsForCover, 0, 0, baseDXDYDZ_[2])

supportsForCover =\
    geompy.MakeCommonList(
        [
            supportsForCover,
            geompy.MakeCutList(
                geompy.TranslateDXDYDZ(
                    geompy.MakeBoxDXDYDZ(
                        chipDXDYDZ_[0],
                        chipDXDYDZ_[1],
                        supportsDXDYDZ_[2]
                    ),
                    0,
                    0,
                    baseDXDYDZ_[2]
                ),
                [
                    geompy.MakeCylinder(receiverCentre[0], OZ, 1.1*plugRrDZ_[0][0], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
                    geompy.MakeCylinder(receiverCentre[1], OZ, 1.1*plugRrDZ_[1][0], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
                    geompy.MakeCylinder(receiverCentre[2], OZ, 1.1*plugRrDZ_[2][0], supportsDXDYDZ_[2] + coverDXDYDZ_[2]),
                ]
            )
        ]
    )

coverSupportsBase = geompy.MakeCompound([cover, supportsForCover, base])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

#geompy.addToStudy( chip, 'chip' )
#geompy.addToStudy( supportsForChip, 'supportsForChip' )
#geompy.addToStudy( base, 'base' )
geompy.addToStudy( chipSupportsBase, 'chipSupportsBase' )

#geompy.addToStudy( cover, 'cover' )
#geompy.addToStudy( supportsForCover, 'supportsForCover' )
geompy.addToStudy( coverSupportsBase, 'coverSupportsBase' )

#geompy.ExportSTL(chipSupportsBase, absoluteCasePath_ + "/media/alexshtil/STORAGE/bubbleStuff/3Dprinting/chipSupportsBase.stl", False, 1e-5, True)
#geompy.ExportSTL(coverSupportsBase, absoluteCasePath_ + "/media/alexshtil/STORAGE/bubbleStuff/3Dprinting/coverSupportsBase.stl", False, 1e-5, True)


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()

