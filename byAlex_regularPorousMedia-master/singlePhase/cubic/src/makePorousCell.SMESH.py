## To do
## 1) Remove the scaleFactor
## 2) Remove the increasingScaleFactor
## 3) Get the primitive cell by MakeCommon
## 4) Rework the graines creation part (rotations, transformation vestors)

import sys
import math
import numpy as np
from decimal import *
getcontext().prec=9
import salome

import GEOM
from salome.geom import geomBuilder
import SALOMEDS

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder
from salome.StdMeshers import StdMeshersBuilder

salome.salome_init()

#For SALOME version lower 9.2.0
#theStudy = salome.myStudy

class PorousCell:
    """"""
    def __init__(self, alpha_ = 70, beta_ = 70, gamma_ = 70, zeta_ = 0.01):
        """"""
        self.alpha_ = Decimal(str(math.radians(alpha_)))
        self.beta_ = Decimal(str(math.radians(beta_)))
        self.gamma_ = Decimal(str(math.radians(gamma_)))

        #intersection parameter
        self.zeta_ = Decimal(str(zeta_))

        #theta_ = math.radians(theta_)
        #alpha_ = math.acos(math.cos(theta_)/math.cos(0.5*theta_))

        self.onePoreCellSize_ = Decimal('1.0')
        self.grainSize_0_ = Decimal('0.5')
        self.grainSize_ = self.grainSize_0_/(Decimal('1.0') - self.zeta_)
        self.atomicPlaneSize_ = Decimal('10.0')

        self.filletFactor_ = Decimal('0.15')

        #Basis of primitive lattice
        self.e_1_ =\
            np.array(
                [
                    Decimal(1.0),
                    Decimal(0.0),
                    Decimal(0.0)
                ]
            )
        self.e_2_ =\
            np.array(
                [
                    1*Decimal(math.cos(self.gamma_)),
                    1*Decimal(math.sin(self.gamma_)),
                    1*Decimal(0.0)
                ]
            )
        self.e_3_ =\
            np.array(
                [
                    1*Decimal(math.cos(self.beta_)),
                    1*Decimal((math.cos(self.alpha_) - math.cos(self.beta_)*math.cos(self.gamma_))/math.sin(self.gamma_)),
                    1*Decimal(math.sqrt(1 - math.cos(self.alpha_)**2 - math.cos(self.beta_)**2 - math.cos(self.gamma_)**2 + 2.0*math.cos(self.alpha_)*math.cos(self.beta_)*math.cos(self.gamma_))/math.sin(self.gamma_))
                ]
            )

        #Primitive vectors of Bravais lattice
        self.a_1_ = self.onePoreCellSize_*self.e_1_
        self.a_2_ = self.onePoreCellSize_*self.e_2_
        self.a_3_ = self.onePoreCellSize_*self.e_3_
        self.primitiveUnitCellVolume_ = np.dot(self.a_1_, np.cross(self.a_2_, self.a_3_))

        #Reciprocal lattice basis
        self.b_1_ = np.cross(self.a_2_, self.a_3_)/self.primitiveUnitCellVolume_
        self.b_2_ = np.cross(self.a_3_, self.a_1_)/self.primitiveUnitCellVolume_
        self.b_3_ = np.cross(self.a_1_, self.a_2_)/self.primitiveUnitCellVolume_
        self.reciprocalUnitCellVolume_ = np.dot(self.b_1_, np.cross(self.b_2_, self.b_3_))

        #Miller index
        #(1, 1, 1)
        self.hkl_1_ =\
            np.array(
                [
                    1,
                    1,
                    1
                ]
            )
        self.hkl_2_ =\
            np.array(
                [
                    -1,
                    1,
                    1
                ]
            )
        self.hkl_3_ =\
            np.array(
                [
                    1,
                    -1,
                    1
                ]
            )

        #(1, 1, 0)
        #self.hkl_1_ =\
            #np.array(
                #[
                    #1,
                    #1,
                    #0
                #]
            #)
        #self.hkl_2_ =\
            #np.array(
                #[
                    #-1,
                    #1,
                    #0
                #]
            #)
        #self.hkl_3_ =\
            #np.array(
                #[
                    #0,
                    #0,
                    #1
                #]
            #)
        self.g_hkl_1_ =\
            self.hkl_1_[0]*self.b_1_ + self.hkl_1_[1]*self.b_2_ + self.hkl_1_[2]*self.b_3_
        self.g_hkl_2_ =\
            self.hkl_2_[0]*self.b_1_ + self.hkl_2_[1]*self.b_2_ + self.hkl_2_[2]*self.b_3_
        self.g_hkl_3_ =\
            self.hkl_3_[0]*self.b_1_ + self.hkl_3_[1]*self.b_2_ + self.hkl_3_[2]*self.b_3_

        self.O = geompy.MakeVertex(0, 0, 0)
        self.OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
        self.OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
        self.OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

        #Bravais lattice basis (primitive vectors)
        self.a_1 = geompy.MakeVectorDXDYDZ(float(self.a_1_[0]), float(self.a_1_[1]), float(self.a_1_[2]))
        self.a_2 = geompy.MakeVectorDXDYDZ(float(self.a_2_[0]), float(self.a_2_[1]), float(self.a_2_[2]))
        self.a_3 = geompy.MakeVectorDXDYDZ(float(self.a_3_[0]), float(self.a_3_[1]), float(self.a_3_[2]))
        #Reciprocal lattice basis
        self.b_1 = geompy.MakeVectorDXDYDZ(float(self.b_1_[0]), float(self.b_1_[1]), float(self.b_1_[2]))
        self.b_2 = geompy.MakeVectorDXDYDZ(float(self.b_2_[0]), float(self.b_2_[1]), float(self.b_2_[2]))
        self.b_3 = geompy.MakeVectorDXDYDZ(float(self.b_3_[0]), float(self.b_3_[1]), float(self.b_3_[2]))

        self.g_hkl_1 = geompy.MakeVectorDXDYDZ(float(self.g_hkl_1_[0]), float(self.g_hkl_1_[1]), float(self.g_hkl_1_[2]))
        self.g_hkl_2 = geompy.MakeVectorDXDYDZ(float(self.g_hkl_2_[0]), float(self.g_hkl_2_[1]), float(self.g_hkl_2_[2]))
        self.g_hkl_3 = geompy.MakeVectorDXDYDZ(float(self.g_hkl_3_[0]), float(self.g_hkl_3_[1]), float(self.g_hkl_3_[2]))

        #Make the grain's centeres
        self.grainCentre_000 =\
            geompy.MakeVertex(
                0,
                0,
                0
            )
        self.grainCentre_100 =\
            geompy.MakeVertex(
                float(self.a_1_[0]),
                float(self.a_1_[1]),
                float(self.a_1_[2])
            )
        self.grainCentre_110 =\
            geompy.MakeVertex(
                float(self.a_1_[0] + self.a_2_[0]),
                float(self.a_1_[1] + self.a_2_[1]),
                float(self.a_1_[2] + self.a_2_[2])
            )
        self.grainCentre_010 =\
            geompy.MakeVertex(
                float(self.a_2_[0]),
                float(self.a_2_[1]),
                float(self.a_2_[2])
            )

        self.grainCentre_001 =\
            geompy.MakeVertex(
                float(self.a_3_[0]),
                float(self.a_3_[1]),
                float(self.a_3_[2])
            )
        self.grainCentre_101 =\
            geompy.MakeVertex(
                float(self.a_1_[0] + self.a_3_[0]),
                float(self.a_1_[1] + self.a_3_[1]),
                float(self.a_1_[2] + self.a_3_[2])
            )
        self.grainCentre_111 =\
            geompy.MakeVertex(
                float(self.a_1_[0] + self.a_2_[0] + self.a_3_[0]),
                float(self.a_1_[1] + self.a_2_[1] + self.a_3_[1]),
                float(self.a_1_[2] + self.a_2_[2] + self.a_3_[2])
            )
        self.grainCentre_011 =\
            geompy.MakeVertex(
                float(self.a_2_[0] + self.a_3_[0]),
                float(self.a_2_[1] + self.a_3_[1]),
                float(self.a_2_[2] + self.a_3_[2])
            )

        geompy.addToStudy( self.O, 'O' )
        geompy.addToStudy( self.OX, 'OX' )
        geompy.addToStudy( self.OY, 'OY' )
        geompy.addToStudy( self.OZ, 'OZ' )

        geompy.addToStudy( self.a_1, 'a_1')
        geompy.addToStudy( self.a_2, 'a_2')
        geompy.addToStudy( self.a_3, 'a_3')

        geompy.addToStudy( self.b_1, 'b_1')
        geompy.addToStudy( self.b_2, 'b_2')
        geompy.addToStudy( self.b_3, 'b_3')

        geompy.addToStudy( self.g_hkl_1, 'g_hkl_1')
        geompy.addToStudy( self.g_hkl_2, 'g_hkl_2')
        geompy.addToStudy( self.g_hkl_3, 'g_hkl_3')

    def makeCell(self):
        """"""
        atomicPlane_1 =\
            geompy.MakeDiskPntVecR(
                self.O,
                self.g_hkl_1,
                float(self.atomicPlaneSize_)
            )\
            if (math.sqrt(np.dot(hkl_1_, hkl_1_)) == 1)\
            else\
            geompy.MakeTranslationVectorDistance(
                geompy.MakeDiskPntVecR(
                    self.O,
                    self.g_hkl_1,
                    float(atomicPlaneSize_)
                ),
                self.g_hkl_1,
                float(2*max(hkl_1_)/math.sqrt(np.dot(g_hkl_1_, g_hkl_1_)))
            )
        #atomicPlane_2 =\
            #geompy.MakeTranslationVectorDistance(
                #atomicPlane_1,
                #g_hkl_1,
                #2.0/math.sqrt(np.dot(g_hkl_1_, g_hkl_1_))
            #)

        #atomicPlane_3 =\
            #geompy.MakeDiskPntVecR(
                #O,
                #g_hkl_2,
                #atomicPlaneSize_
            #)
            ##)\
            ##if (math.sqrt(np.dot(hkl_2_, hkl_2_)) == 1)\
            ##else\
            ##geompy.MakeTranslationVectorDistance(
                ##geompy.MakeDiskPntVecR(
                    ##O,
                    ##g_hkl_2,
                    ##atomicPlaneSize_
                ##),
                ##g_hkl_2,
                ##0.0
                ###-max(hkl_2_)*1.0/math.sqrt(np.dot(g_hkl_2_, g_hkl_2_))
            ##)
        #atomicPlane_4 =\
            #geompy.MakeTranslationVectorDistance(
                #atomicPlane_3,
                #g_hkl_2,
                #2.0/math.sqrt(np.dot(g_hkl_2_, g_hkl_2_))
            #)

        #atomicPlane_5 =\
            #geompy.MakeDiskPntVecR(
                #O,
                #g_hkl_3,
                #atomicPlaneSize_
            #)
            ##)\
            ##if (math.sqrt(np.dot(hkl_3_, hkl_3_)) == 1)\
            ##else\
            ##geompy.MakeTranslationVectorDistance(
                ##geompy.MakeDiskPntVecR(
                    ##O,
                    ##g_hkl_3,
                    ##atomicPlaneSize_
                ##),
                ##g_hkl_3,
                ##0.0
                ###-max(hkl_3_)*1.0/math.sqrt(np.dot(g_hkl_3_, g_hkl_3_))
            ##)
        #atomicPlane_6 =\
            #geompy.MakeTranslationVectorDistance(
                #atomicPlane_5,
                #g_hkl_3,
                #2.0/math.sqrt(np.dot(g_hkl_3_, g_hkl_3_))
            #)

        ##block =\
            ##geompy.GetBlockNearPoint(
                ##geompy.MakePartition(
                    ##[geompy.MakeBoxDXDYDZ(5, 5, 5)],
                    ##[atomicPlane_1, atomicPlane_2, atomicPlane_3, atomicPlane_4, atomicPlane_5, atomicPlane_6],
                    ##[], [], geompy.ShapeType["SOLID"], 0, [], 0
                ##),
                ##geompy.MakeVertex(
                    ##a_1_[0] + a_2_[0] + a_3_[0],
                    ##a_1_[1] + a_2_[1] + a_3_[1],
                    ##a_1_[2] + a_2_[2] + a_3_[2]
                ##)
            ##)

        #partition =\
            #geompy.MakePartition(
                #[geompy.MakeBoxDXDYDZ(5, 5, 5)],
                #[atomicPlane_1, atomicPlane_2, atomicPlane_3, atomicPlane_4, atomicPlane_5, atomicPlane_6],
                #[], [], geompy.ShapeType["SOLID"], 0, [], 0
            #)

        #return(partition)

    def makeGrains(self):
        """"""
        #Make the grains with its rotations
        #Grain #1
        grain_000 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_000,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(self.a_1_[0] + self.a_2_[0]),
                    float(self.a_1_[1] + self.a_2_[1]),
                    float(self.a_1_[2] + self.a_2_[2])
                )
            )
            + math.pi
        )
        geompy.Rotate(
            grain_000,
            geompy.MakeVectorDXDYDZ(
                float(self.a_1_[0] + self.a_2_[0]),
                float(self.a_1_[1] + self.a_2_[1]),
                float(self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )

        #Grain #2
        grain_100 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_100,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(-self.a_1_[0] + self.a_2_[0]),
                    float(-self.a_1_[1] + self.a_2_[1]),
                    float(-self.a_1_[2] + self.a_2_[2])
                )
            )
            + math.pi
        )
        geompy.Rotate(
            grain_100,
            geompy.MakeVectorDXDYDZ(
                float(-self.a_1_[0] + self.a_2_[0]),
                float(-self.a_1_[1] + self.a_2_[1]),
                float(-self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_100, self.grainCentre_000, self.grainCentre_100)

        #Grain #3
        grain_110 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_110,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(self.a_1_[0] + self.a_2_[0]),
                    float(self.a_1_[1] + self.a_2_[1]),
                    float(self.a_1_[2] + self.a_2_[2])
                )
            )
        )
        geompy.Rotate(
            grain_110,
            geompy.MakeVectorDXDYDZ(
                float(self.a_1_[0] + self.a_2_[0]),
                float(self.a_1_[1] + self.a_2_[1]),
                float(self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_110, self.grainCentre_000, self.grainCentre_110)

        #Grain #4
        grain_010 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_010,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(-self.a_1_[0] + self.a_2_[0]),
                    float(-self.a_1_[1] + self.a_2_[1]),
                    float(-self.a_1_[2] + self.a_2_[2])
                )
            )
        )
        geompy.Rotate(
            grain_010,
            geompy.MakeVectorDXDYDZ(
                float(-self.a_1_[0] + self.a_2_[0]),
                float(-self.a_1_[1] + self.a_2_[1]),
                float(-self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_010, self.grainCentre_000, self.grainCentre_010)

        #Grain #5
        grain_001 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_001,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(self.a_1_[0] + self.a_2_[0]),
                    float(self.a_1_[1] + self.a_2_[1]),
                    float(self.a_1_[2] + self.a_2_[2])
                )
            )
            + math.pi
        )
        geompy.Rotate(
            grain_001,
            geompy.MakeVectorDXDYDZ(
                float(self.a_1_[0] + self.a_2_[0]),
                float(self.a_1_[1] + self.a_2_[1]),
                float(self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_001, self.grainCentre_000, self.grainCentre_001)

        #Grain #6
        grain_101 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_101,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(-self.a_1_[0] + self.a_2_[0]),
                    float(-self.a_1_[1] + self.a_2_[1]),
                    float(-self.a_1_[2] + self.a_2_[2])
                )
            )
            + math.pi
        )
        geompy.Rotate(
            grain_101,
            geompy.MakeVectorDXDYDZ(
                float(-self.a_1_[0] + self.a_2_[0]),
                float(-self.a_1_[1] + self.a_2_[1]),
                float(-self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_101, self.grainCentre_000, self.grainCentre_101)

        #Grain #7
        grain_111 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_111,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(self.a_1_[0] + self.a_2_[0]),
                    float(self.a_1_[1] + self.a_2_[1]),
                    float(self.a_1_[2] + self.a_2_[2])
                )
            )
        )
        geompy.Rotate(
            grain_111,
            geompy.MakeVectorDXDYDZ(
                float(self.a_1_[0] + self.a_2_[0]),
                float(self.a_1_[1] + self.a_2_[1]),
                float(self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_111, self.grainCentre_000, self.grainCentre_111)

        #Grain #8
        grain_011 = geompy.MakeSpherePntR(self.grainCentre_000, float(self.grainSize_))
        geompy.Rotate(
            grain_011,
            self.b_3,
            geompy.GetAngleRadiansVectors(
                self.a_1,
                geompy.MakeVectorDXDYDZ(
                    float(-self.a_1_[0] + self.a_2_[0]),
                    float(-self.a_1_[1] + self.a_2_[1]),
                    float(-self.a_1_[2] + self.a_2_[2])
                )
            )
        )
        geompy.Rotate(
            grain_011,
            geompy.MakeVectorDXDYDZ(
                float(-self.a_1_[0] + self.a_2_[0]),
                float(-self.a_1_[1] + self.a_2_[1]),
                float(-self.a_1_[2] + self.a_2_[2])
            ),
            0.5*math.pi
        )
        geompy.TranslateTwoPoints(grain_011, self.grainCentre_000, self.grainCentre_011)

        geompy.addToStudy( grain_000, 'grain_000')
        geompy.addToStudy( grain_100, 'grain_100')
        geompy.addToStudy( grain_110, 'grain_110')
        geompy.addToStudy( grain_010, 'grain_010')
        geompy.addToStudy( grain_001, 'grain_001')
        geompy.addToStudy( grain_101, 'grain_101')
        geompy.addToStudy( grain_111, 'grain_111')
        geompy.addToStudy( grain_011, 'grain_011')

        return(
            [
                grain_000,
                grain_100,
                grain_110,
                grain_010,
                grain_001,
                grain_101,
                grain_111,
                grain_011
            ]
        )

    def makePrimitivePore(self, grains):
        """"""
        #Building the skeleton
        #building the porous volume
        #Make faces of porous volume
        onePoreBoxLeftFace = geompy.MakeQuad4Vertices(self.grainCentre_000, self.grainCentre_010, self.grainCentre_011, self.grainCentre_001)
        onePoreBoxRightFace = geompy.MakeQuad4Vertices(self.grainCentre_100, self.grainCentre_110, self.grainCentre_111, self.grainCentre_101)
        onePoreBoxOppositeFace = geompy.MakeQuad4Vertices(self.grainCentre_010, self.grainCentre_110, self.grainCentre_111, self.grainCentre_011)
        onePoreBoxFrontFace = geompy.MakeQuad4Vertices(self.grainCentre_000, self.grainCentre_100, self.grainCentre_101, self.grainCentre_001)
        onePoreBoxTopFace = geompy.MakeQuad4Vertices(self.grainCentre_001, self.grainCentre_101, self.grainCentre_111, self.grainCentre_011)
        onePoreBoxBottomFace = geompy.MakeQuad4Vertices(self.grainCentre_000, self.grainCentre_100, self.grainCentre_110, self.grainCentre_010)

        primitivePoreBox = geompy.MakeHexa(onePoreBoxLeftFace, onePoreBoxRightFace, onePoreBoxOppositeFace, onePoreBoxFrontFace, onePoreBoxTopFace, onePoreBoxBottomFace)

        fusedGrains = geompy.MakeFuseList(grains)

        primitivePore = geompy.MakeCut(primitivePoreBox, fusedGrains)

        #primitivePore = geompy.RemoveExtraEdges(primitivePore, False)

        #Get edges IDs with fillet
        checkFilletSwitcher = 0

        #allEdgesList = []
        allEdgesWithFillet = []
        allEdgesIDsWithFillet = []
        allEdgesWithoutFillet = []

        #allEdgesList = geompy.SubShapeAll(primitivePore, geompy.ShapeType["EDGE"])
        allEdgesWithFillet = geompy.CreateGroup(primitivePore, geompy.ShapeType["EDGE"])

        geompy.UnionList(allEdgesWithFillet, geompy.SubShapeAll(primitivePore, geompy.ShapeType["EDGE"]))

        #Make middle points on boundaries for searching boundaries
        #onePoreBoxMiddlePointOnLeftFace = geompy.MakeVertexOnSurface(onePoreBoxLeftFace, 0.5, 0.5)
        #onePoreBoxMiddlePointOnRightFace = geompy.MakeVertexOnSurface(onePoreBoxRightFace, 0.5, 0.5)
        #onePoreBoxMiddlePointOnOppositeFace = geompy.MakeVertexOnSurface(onePoreBoxOppositeFace, 0.5, 0.5)
        #onePoreBoxMiddlePointOnFrontFace = geompy.MakeVertexOnSurface(onePoreBoxFrontFace, 0.5, 0.5)
        #onePoreBoxMiddlePointOnTopFace = geompy.MakeVertexOnSurface(onePoreBoxTopFace, 0.5, 0.5)
        #onePoreBoxMiddlePointOnBottomFace = geompy.MakeVertexOnSurface(onePoreBoxBottomFace, 0.5, 0.5)

        #Make normals for searching boundaries
        #onePoreBoxNormalVectorOnLeftFace = geompy.GetNormal(onePoreBoxLeftFace)
        #onePoreBoxNormalVectorOnRightFace = geompy.GetNormal(onePoreBoxRightFace)
        #onePoreBoxNormalVectorOnOppositeFace = geompy.GetNormal(onePoreBoxOppositeFace)
        #onePoreBoxNormalVectorOnFrontFace = geompy.GetNormal(onePoreBoxFrontFace)
        #onePoreBoxNormalVectorOnTopFace = geompy.GetNormal(onePoreBoxTopFace)
        #onePoreBoxNormalVectorOnBottomFace = geompy.GetNormal(onePoreBoxBottomFace)

        allEdgesWithoutFillet.extend(
            geompy.GetShapesOnPlaneWithLocation(
                primitivePore,
                geompy.ShapeType["EDGE"],
                geompy.GetNormal(onePoreBoxLeftFace),
                geompy.MakeVertexOnSurface(onePoreBoxLeftFace, 0.5, 0.5),
                GEOM.ST_ON
            )
        )
        allEdgesWithoutFillet.extend(
            geompy.GetShapesOnPlaneWithLocation(
                primitivePore,
                geompy.ShapeType["EDGE"],
                geompy.GetNormal(onePoreBoxRightFace),
                geompy.MakeVertexOnSurface(onePoreBoxRightFace, 0.5, 0.5),
                GEOM.ST_ON
            )
        )
        allEdgesWithoutFillet.extend(
            geompy.GetShapesOnPlaneWithLocation(
                primitivePore,
                geompy.ShapeType["EDGE"],
                geompy.GetNormal(onePoreBoxOppositeFace),
                geompy.MakeVertexOnSurface(onePoreBoxOppositeFace, 0.5, 0.5),
                GEOM.ST_ON
            )
        )
        allEdgesWithoutFillet.extend(
            geompy.GetShapesOnPlaneWithLocation(
                primitivePore,
                geompy.ShapeType["EDGE"],
                geompy.GetNormal(onePoreBoxFrontFace),
                geompy.MakeVertexOnSurface(onePoreBoxFrontFace, 0.5, 0.5),
                GEOM.ST_ON
            )
        )
        allEdgesWithoutFillet.extend(
            geompy.GetShapesOnPlaneWithLocation(
                primitivePore,
                geompy.ShapeType["EDGE"],
                geompy.GetNormal(onePoreBoxTopFace),
                geompy.MakeVertexOnSurface(onePoreBoxTopFace, 0.5, 0.5),
                GEOM.ST_ON
            )
        )
        allEdgesWithoutFillet.extend(
            geompy.GetShapesOnPlaneWithLocation(
                primitivePore,
                geompy.ShapeType["EDGE"],
                geompy.GetNormal(onePoreBoxBottomFace),
                geompy.MakeVertexOnSurface(onePoreBoxBottomFace, 0.5, 0.5),
                GEOM.ST_ON
            )
        )

        geompy.DifferenceList(allEdgesWithFillet, allEdgesWithoutFillet)
        allEdgesIDsWithFillet = geompy.GetObjectIDs(allEdgesWithFillet)
        allEdgesWithFillet = geompy.ExtractShapes(allEdgesWithFillet, geompy.ShapeType["EDGE"], True)

        if (
            float(self.filletFactor_) != 0 and
            intersectionParameter_ != 0
        ):
            #Make fillet
            print("filletRadius=" + str(float(self.filletFactor_) * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])) + "\n")

            #if(
                #float(self.filletFactor_) * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]) > 0.02
            #):
            fusedFilletedGrains = geompy.MakeFilletAll(fusedGrains, float(self.filletFactor_) * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

            #fusedFilletedGrains = geompy.MakeChamferAll(fusedGrains, float(self.filletFactor_) * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

            #For debuging
            #geompy.addToStudy( fusedFilletedGrains, 'fusedFilletedGrains')

            primitivePore = geompy.MakeCut(primitivePoreBox, fusedFilletedGrains)

            #primitivePore = geompy.RemoveExtraEdges(primitivePore, False)
            
            primitivePore = geompy.UnionFaces(primitivePore)

            checkFilletSwitcher = 1
            #else:
                #print("Can not make a fillet!\n")

        geompy.addToStudy( fusedGrains, 'fusedGrains')
        geompy.addToStudy( fusedFilletedGrains, 'fusedFilletedGrains')
        geompy.addToStudy( primitivePoreBox, 'primitivePoreBox')
        geompy.addToStudy( primitivePore, 'primitivePore')

        return(primitivePore)

    def printParameters(self):
        """"""
        print('alpha_ = ' + str(math.degrees(self.alpha_)))
        print('beta_ = ' + str(math.degrees(self.beta_)))
        print('gamma_ = ' + str(math.degrees(self.gamma_)) + '\n')

        print('zeta_ = ' + str(self.zeta_) + '\n')

        print('Generating unit vectors:')
        print('e_1_ = (' + str(self.e_1_[0]) + ', ' + str(self.e_1_[1]) + ', ' + str(self.e_1_[2]) + ')')
        print('e_2_ = (' + str(self.e_2_[0]) + ', ' + str(self.e_2_[1]) + ', ' + str(self.e_2_[2]) + ')')
        print('e_3_ = (' + str(self.e_3_[0]) + ', ' + str(self.e_3_[1]) + ', ' + str(self.e_3_[2]) + ')' + '\n')

        print('Primitive vectors of Bravais lattice:')
        print('a_1_ = (' + str(self.a_1_[0]) + ', ' + str(self.a_1_[1]) + ', ' + str(self.a_1_[2]) + ')')
        print('a_2_ = (' + str(self.a_2_[0]) + ', ' + str(self.a_2_[1]) + ', ' + str(self.a_2_[2]) + ')')
        print('a_3_ = (' + str(self.a_3_[0]) + ', ' + str(self.a_3_[1]) + ', ' + str(self.a_3_[2]) + ')')
        print('primitiveUnitCellVolume_ = ' + str(self.primitiveUnitCellVolume_) + '\n')

        print('Reciprocal cell basis:')
        print('b_1_ = (' + str(self.b_1_[0]) + ', ' + str(self.b_1_[1]) + ', ' + str(self.b_1_[2]) + ')')
        print('b_2_ = (' + str(self.b_2_[0]) + ', ' + str(self.b_2_[1]) + ', ' + str(self.b_2_[2]) + ')')
        print('b_3_ = (' + str(self.b_3_[0]) + ', ' + str(self.b_3_[1]) + ', ' + str(self.b_3_[2]) + ')')
        print('reciprocalUnitCellVolume_ = ' + str(self.reciprocalUnitCellVolume_) + '\n')
        
        print(
            'primitiveUnitCellVolume_*reciprocalUnitCellVolume_ = '
            + str(self.primitiveUnitCellVolume_*self.reciprocalUnitCellVolume_) + '\n'
        )

geompy = geomBuilder.New()

cubic = PorousCell(70, 70, 70, 0.01)

cubic.printParameters()

grains = cubic.makeGrains()

cubic.makePrimitivePore(grains)











####################################################
##            Begin of variables section          ##
####################################################

####################################################
##                  Interface part                ##
####################################################
#absolute_case_path_ = sys.argv[1]
##absolute_case_path_ = sys.argv[0]
##absolute_case_path_ = absolute_case_path_[0:absolute_case_path_.rfind("/")]
##absolute_case_path_ = absolute_case_path_[0:absolute_case_path_.rfind("/")]

#print("\n***********************************")
#print("ABSOLUTE CASE PATH: " + absolute_case_path_)

#theta_ = float(sys.argv[2])
#intersectionParameter_ = float(sys.argv[3])

#theta_ = math.radians(theta_)
#alpha_ = math.acos(math.cos(theta_) / math.cos(0.5 * theta_))

##scaleFactor = float(sys.argv[4])

##increasingScaleFactor = float(sys.argv[5])

##UTranslationDirectionNbTimes_ = int(sys.argv[6])
##VTranslationDirectionNbTimes_ = int(sys.argv[7])
##WTranslationDirectionNbTimes_ = int(sys.argv[8])

#increasingScaleFactor = float(sys.argv[4])

#UTranslationDirectionNbTimes_ = int(sys.argv[5])
#VTranslationDirectionNbTimes_ = int(sys.argv[6])
#WTranslationDirectionNbTimes_ = int(sys.argv[7])

####################################################
##                  Debuging part                 ##
####################################################
#absolute_case_path_ = ""

theta_ = 70.0
intersectionParameter_ = 0.01

theta_ = math.radians(theta_)
alpha_ = math.acos(math.cos(theta_)/math.cos(0.5*theta_))

UTranslationDirectionNbTimes_ = 2
VTranslationDirectionNbTimes_ = 2
WTranslationDirectionNbTimes_ = 2

####################################################
##                   Global part                  ##
####################################################


####################################################
##             End of variables section           ##
####################################################

###
### GEOM component
###

#For SALOME version lower 9.2.0
#geompy = geomBuilder.New(theStudy)

'''

#building the porous volume
#Make faces of porous volume
onePoreBoxLeftFace = geompy.MakeQuad4Vertices(grainCentre_000, grainCentre_010, grainCentre_011, grainCentre_001)
onePoreBoxRightFace = geompy.MakeQuad4Vertices(grainCentre_100, grainCentre_110, grainCentre_111, grainCentre_101)
onePoreBoxOppositeFace = geompy.MakeQuad4Vertices(grainCentre_010, grainCentre_110, grainCentre_111, grainCentre_011)
onePoreBoxFrontFace = geompy.MakeQuad4Vertices(grainCentre_000, grainCentre_100, grainCentre_101, grainCentre_001)
onePoreBoxTopFace = geompy.MakeQuad4Vertices(grainCentre_001, grainCentre_101, grainCentre_111, grainCentre_011)
onePoreBoxBottomFace = geompy.MakeQuad4Vertices(grainCentre_000, grainCentre_100, grainCentre_110, grainCentre_010)

onePoreBox = geompy.MakeHexa(onePoreBoxLeftFace, onePoreBoxRightFace, onePoreBoxOppositeFace, onePoreBoxFrontFace, onePoreBoxTopFace, onePoreBoxBottomFace)

#For debuging
geompy.addToStudy( onePoreBox, 'onePoreBox')

fusedGrains = geompy.MakeFuseList([grain_000, grain_100, grain_110, grain_010, grain_001, grain_101, grain_111, grain_011, mirroredGrain_1, mirroredGrain_2, mirroredGrain_3, mirroredGrain_4, mirroredGrain_5, mirroredGrain_6, mirroredGrain_7, mirroredGrain_8])

#For debuging
geompy.addToStudy( fusedGrains, 'fusedGrains')

onePore = geompy.MakeCut(onePoreBox, fusedGrains)

onePore = geompy.RemoveExtraEdges(onePore, False)

#Get edges IDs with fillet
checkFilletSwitcher = 0

allEdgesList = []
allEdgesWithFillet = []
allEdgesIDsWithFillet = []
allEdgesWithoutFillet = []

allEdgesList = geompy.SubShapeAll(onePore, geompy.ShapeType["EDGE"])
allEdgesWithFillet = geompy.CreateGroup(onePore, geompy.ShapeType["EDGE"])

geompy.UnionList(allEdgesWithFillet, allEdgesList)

#Make middle points on boundaries for searching boundaries
onePoreBoxMiddlePointOnLeftFace = geompy.MakeVertexOnSurface(onePoreBoxLeftFace, 0.5, 0.5)
onePoreBoxMiddlePointOnRightFace = geompy.MakeVertexOnSurface(onePoreBoxRightFace, 0.5, 0.5)
onePoreBoxMiddlePointOnOppositeFace = geompy.MakeVertexOnSurface(onePoreBoxOppositeFace, 0.5, 0.5)
onePoreBoxMiddlePointOnFrontFace = geompy.MakeVertexOnSurface(onePoreBoxFrontFace, 0.5, 0.5)
onePoreBoxMiddlePointOnTopFace = geompy.MakeVertexOnSurface(onePoreBoxTopFace, 0.5, 0.5)
onePoreBoxMiddlePointOnBottomFace = geompy.MakeVertexOnSurface(onePoreBoxBottomFace, 0.5, 0.5)

#Make normals for searching boundaries
onePoreBoxNormalVectorOnLeftFace = geompy.GetNormal(onePoreBoxLeftFace)
onePoreBoxNormalVectorOnRightFace = geompy.GetNormal(onePoreBoxRightFace)
onePoreBoxNormalVectorOnOppositeFace = geompy.GetNormal(onePoreBoxOppositeFace)
onePoreBoxNormalVectorOnFrontFace = geompy.GetNormal(onePoreBoxFrontFace)
onePoreBoxNormalVectorOnTopFace = geompy.GetNormal(onePoreBoxTopFace)
onePoreBoxNormalVectorOnBottomFace = geompy.GetNormal(onePoreBoxBottomFace)

allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnLeftFace, onePoreBoxMiddlePointOnLeftFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnRightFace, onePoreBoxMiddlePointOnRightFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnOppositeFace, onePoreBoxMiddlePointOnOppositeFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnFrontFace, onePoreBoxMiddlePointOnFrontFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnTopFace, onePoreBoxMiddlePointOnTopFace, GEOM.ST_ON))
allEdgesWithoutFillet.extend(geompy.GetShapesOnPlaneWithLocation(onePore, geompy.ShapeType["EDGE"], onePoreBoxNormalVectorOnBottomFace, onePoreBoxMiddlePointOnBottomFace, GEOM.ST_ON))

geompy.DifferenceList(allEdgesWithFillet, allEdgesWithoutFillet)
allEdgesIDsWithFillet = geompy.GetObjectIDs(allEdgesWithFillet)
allEdgesWithFillet = geompy.ExtractShapes(allEdgesWithFillet, geompy.ShapeType["EDGE"], True)

if (
    filletFactor_ != 0 and
    intersectionParameter_ != 0
):
    #Make fillet
    print("filletRadius=" + str(filletFactor_ * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])) + "\n")

    #if(
        #filletFactor_ * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]) > 0.02
    #):
    fusedFilletedGrains = geompy.MakeFilletAll(fusedGrains, filletFactor_ * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

    #fusedFilletedGrains = geompy.MakeChamferAll(fusedGrains, filletFactor_ * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))]))

    #For debuging
    geompy.addToStudy( fusedFilletedGrains, 'fusedFilletedGrains')

    onePore = geompy.MakeCut(onePoreBox, fusedFilletedGrains)

    onePore = geompy.RemoveExtraEdges(onePore, False)

    checkFilletSwitcher = 1
    #else:
        #print("Can not make a fillet!\n")

#For debuging
geompy.addToStudy( onePore, 'onePore')

##One pore centre of mass calculation\
#onePoreCentreOfMassFile = open(absolute_case_path_ + "/onePoreCentreOfMass", "w", 0)
onePoreCentreOfMass = geompy.MakeCDG(onePore)
##onePoreCentreOfMassFile.write(
    ##"centre ("
    ##+ str(scaleFactor / increasingScaleFactor * geompy.PointCoordinates(onePoreCentreOfMass)[0]) + " "
    ##+ str(scaleFactor / increasingScaleFactor * geompy.PointCoordinates(onePoreCentreOfMass)[1]) + " "
    ##+ str(scaleFactor / increasingScaleFactor * geompy.PointCoordinates(onePoreCentreOfMass)[2]) + ");"
##)
#onePoreCentreOfMassFile.write(
    #"centre ("
    #+ str(geompy.PointCoordinates(onePoreCentreOfMass)[0]) + " "
    #+ str(geompy.PointCoordinates(onePoreCentreOfMass)[1]) + " "
    #+ str(geompy.PointCoordinates(onePoreCentreOfMass)[2]) + ");"
#)
#onePoreCentreOfMassFile.close()

#Make translation vectors for building porous cell
UTranslationDirectionVector = geompy.MakeVector(grainCentre_000, grainCentre_100)
VTranslationDirectionVector = geompy.MakeVector(grainCentre_000, grainCentre_010)
WTranslationDirectionVector = geompy.MakeVector(grainCentre_000, grainCentre_001)

#Building porous cell
porousCell = geompy.MakeMultiTranslation2D(onePore, UTranslationDirectionVector, onePoreCellSize_, UTranslationDirectionNbTimes_, VTranslationDirectionVector, onePoreCellSize_, VTranslationDirectionNbTimes_)

porousCell = geompy.MakeMultiTranslation1D(porousCell, WTranslationDirectionVector, onePoreCellSize_, WTranslationDirectionNbTimes_)

listOfSolids = []
listOfSolids = geompy.SubShapeAllSortedCentres(porousCell, geompy.ShapeType["SOLID"])

porousCell = geompy.MakeFuseList(listOfSolids)

porousCell = geompy.MakePartition(
        [porousCell],
        [
            geompy.MakeMultiTranslation1D(
                geompy.MakeMultiTranslation2D(
                    geompy.MakeScaleTransform(
                        onePoreBox,
                        onePoreCentreOfMass,
                        0.85
                    ),
                    UTranslationDirectionVector,
                    onePoreCellSize_,
                    UTranslationDirectionNbTimes_,
                    VTranslationDirectionVector,
                    onePoreCellSize_,
                    VTranslationDirectionNbTimes_
                ),
                WTranslationDirectionVector,
                onePoreCellSize_,
                WTranslationDirectionNbTimes_
            )
        ],
        [],
        [],
        geompy.ShapeType["SOLID"],
        0,
        [],
        0
    )

listOfSolids = geompy.SubShapeAllSortedCentres(porousCell, geompy.ShapeType["SOLID"])

throats = geompy.CreateGroup(porousCell, geompy.ShapeType["SOLID"])
geompy.UnionList(throats, listOfSolids)

geompy.DifferenceList(
    throats,
    [
        geompy.GetShapesNearPoint(porousCell, vertexIndex, geompy.ShapeType["SOLID"])
        for vertexIndex in
        geompy.SubShapeAllSortedCentres(
            geompy.MakeMultiTranslation1D(
                geompy.MakeMultiTranslation2D(
                    onePoreCentreOfMass,
                    UTranslationDirectionVector,
                    onePoreCellSize_,
                    UTranslationDirectionNbTimes_,
                    VTranslationDirectionVector,
                    onePoreCellSize_,
                    VTranslationDirectionNbTimes_
                ),
                WTranslationDirectionVector,
                onePoreCellSize_,
                WTranslationDirectionNbTimes_
            ),
            geompy.ShapeType["VERTEX"]
        )
    ]
)

#geompy.MakeCommon(porousCell, block, True)

#For debuging
geompy.addToStudy( porousCell, 'porousCell')
#geompy.addToStudyInFather( porousCell, throats, 'throats' )

#Make the porous cell box's vertexes
porousCellBoxVertexOnBottom_1 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * 0 * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(0)
)
porousCellBoxVertexOnBottom_2 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(0)
)
porousCellBoxVertexOnBottom_3 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ + UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.cos(theta_),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(theta_),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(0)
)
porousCellBoxVertexOnBottom_4 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.cos(theta_),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(theta_),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(0)
)

porousCellBoxVertexOnTop_1 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.cos(0.5 * theta_) + UTranslationDirectionNbTimes_ * 0 * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.sin(0.5 * theta_) + VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(alpha_)
)
porousCellBoxVertexOnTop_2 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.cos(0.5 * theta_) + UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.cos(0),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.sin(0.5 * theta_) + VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(0),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(alpha_)
)
porousCellBoxVertexOnTop_3 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.cos(0.5 * theta_) + UTranslationDirectionNbTimes_ * onePoreCellSize_ + UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.cos(theta_),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.sin(0.5 * theta_) + VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(theta_),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(alpha_)
)
porousCellBoxVertexOnTop_4 = geompy.MakeVertex(
    UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.cos(0.5 * theta_) + UTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.cos(theta_),
    VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(alpha_) * math.sin(0.5 * theta_) + VTranslationDirectionNbTimes_ * onePoreCellSize_ * math.cos(0) * math.sin(theta_),
    WTranslationDirectionNbTimes_ * onePoreCellSize_ * math.sin(alpha_)
)

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

#For debuging
geompy.addToStudy( porousCellBox, 'porousCellBox')

facesListOnLeftSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnLeftFace, porousCellBoxMiddlePointOnLeftFace, GEOM.ST_ON)
facesListOnRightSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnRightFace, porousCellBoxMiddlePointOnRightFace, GEOM.ST_ON)
facesListOnOppositeSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnOppositeFace, porousCellBoxMiddlePointOnOppositeFace, GEOM.ST_ON)
facesListOnFrontSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnFrontFace, porousCellBoxMiddlePointOnFrontFace, GEOM.ST_ON)
facesListOnTopSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnTopFace, porousCellBoxMiddlePointOnTopFace, GEOM.ST_ON)
facesListOnBottomSide = geompy.GetShapesOnPlaneWithLocation(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnBottomFace, porousCellBoxMiddlePointOnBottomFace, GEOM.ST_ON)

##Get fases IDs without viscous layer
allFacesIDsWithoutViscousLayers = []

allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnLeftFace, porousCellBoxMiddlePointOnLeftFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnRightFace, porousCellBoxMiddlePointOnRightFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnOppositeFace, porousCellBoxMiddlePointOnOppositeFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnFrontFace, porousCellBoxMiddlePointOnFrontFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnTopFace, porousCellBoxMiddlePointOnTopFace, GEOM.ST_ON))
allFacesIDsWithoutViscousLayers.extend(geompy.GetShapesOnPlaneWithLocationIDs(porousCell, geompy.ShapeType["FACE"], porousCellBoxNormalVectorOnBottomFace, porousCellBoxMiddlePointOnBottomFace, GEOM.ST_ON))

#Create groups
#Get faces of non-skeleton wall
nonSkeletonFaces = []

nonSkeletonFaces.extend(facesListOnLeftSide)
nonSkeletonFaces.extend(facesListOnRightSide)
nonSkeletonFaces.extend(facesListOnOppositeSide)
nonSkeletonFaces.extend(facesListOnFrontSide)
nonSkeletonFaces.extend(facesListOnTopSide)
nonSkeletonFaces.extend(facesListOnBottomSide)

listOfSolids = geompy.SubShapeAllSortedCentres(
    geompy.MakeMultiTranslation1D(
        geompy.MakeMultiTranslation2D(
            geompy.MakeScaleTransform(
                onePoreBox,
                onePoreCentreOfMass,
                0.85
            ),
            UTranslationDirectionVector,
            onePoreCellSize_,
            UTranslationDirectionNbTimes_,
            VTranslationDirectionVector,
            onePoreCellSize_,
            VTranslationDirectionNbTimes_
        ),
        WTranslationDirectionVector,
        onePoreCellSize_,
        WTranslationDirectionNbTimes_
    ),
    geompy.ShapeType["SOLID"]
)

nonSkeletonFaces.extend(
    sum(
        [
            geompy.GetShapesOnShape(
                solidIndex,
                porousCell,
                geompy.ShapeType["FACE"],
                GEOM.ST_ON
            )
            for solidIndex in listOfSolids
        ],
        []
    )
)

##Get faces of skeleton wall
allFacesList = []

allFacesList = geompy.SubShapeAll(porousCell, geompy.ShapeType["FACE"])
skeletonWall = geompy.CreateGroup(porousCell, geompy.ShapeType["FACE"])
geompy.UnionList(skeletonWall, allFacesList)
geompy.DifferenceList(skeletonWall, nonSkeletonFaces)

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

geompy.addToStudyInFather( porousCell, skeletonWall, 'skeletonWall' )
geompy.addToStudyInFather( porousCell, leftSide, 'leftSide' )
geompy.addToStudyInFather( porousCell, rightSide, 'rightSide' )
geompy.addToStudyInFather( porousCell, oppositeSide, 'oppositeSide' )
geompy.addToStudyInFather( porousCell, frontSide, 'frontSide' )
geompy.addToStudyInFather( porousCell, topSide, 'topSide' )
geompy.addToStudyInFather( porousCell, bottomSide, 'bottomSide' )

##Save geometry to file
#print("Save geometry to " + absolute_case_path_ + "/porousCellGeometry.hdf\n")

#salome.myStudyManager.SaveAs(absolute_case_path_ + "/porousCellGeometry.hdf", salome.myStudy, False)

####
#### SMESH component
####

#####################################################
###         Begin of mesh variables section        ##
#####################################################
##globalNETGEN_2D_MaxSize = 0.05 * geompy.BasicProperties(onePore)[1] ** (1.0 / 2.0)
##globalNETGEN_2D_MinSize = 0.001 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
##patchesNETGEN_2D_MaxSize = 0.05 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
##patchesNETGEN_2D_MinSize = 0.001 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
##NETGEN_2D_QuadAllowed = 0
##growthRate2D = 0.1

##globalLocalLengthSize = 0.005 * max([geompy.BasicProperties(allEdgesWithFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithFillet))])
##localLengthSize = 0.01 * max([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])

#globalLocalLengthSize = 0.02 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)
#localLengthSize = 0.06 * min([geompy.BasicProperties(allEdgesWithoutFillet[edgeIndex])[0] for edgeIndex in range(len(allEdgesWithoutFillet))])

#globalNETGEN_2D_MaxSize = 2.5 * localLengthSize
#globalNETGEN_2D_MinSize = 1.0 * localLengthSize
#throatsNETGEN_2D_MaxSize = 1.0 * localLengthSize
#throatsNETGEN_2D_MinSize = 1.0 * localLengthSize
#patchesNETGEN_2D_MaxSize = 1.0 * localLengthSize
#patchesNETGEN_2D_MinSize = 1.0 * localLengthSize
#NETGEN_2D_QuadAllowed = 0
#globalGrowthRate2D = 0.02
#throatsGrowthRate2D = 1e-4
#patchesGrowthRate2D = 1e-4

#globalNETGEN_3D_MaxSize = 4.0 * globalLocalLengthSize
#globalNETGEN_3D_MinSize = 1.0 * localLengthSize
#throatsNETGEN_3D_MaxSize = 2.0 * localLengthSize
#throatsNETGEN_3D_MinSize = 1.0 * localLengthSize
#globalGrowthRate3D = 0.05
#throatsGrowthRate3D = 0.02

##viscousLayersSize3D = 10.0 * 1e-3 * onePoreCellSize_
#globalViscousLayersSize3D = min(onePoreCellSize_ / increasingScaleFactor, 0.025 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0))
#globalViscousLayersNumber = 3
#globalViscousLayersGrowth = 1.0

#throatsViscousLayersSize3D = min(onePoreCellSize_ / increasingScaleFactor, 0.025 * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0))
#throatsViscousLayersNumber = 3
#throatsViscousLayersGrowth = 1.0

#####################################################
###          End of mesh variables section         ##
#####################################################

#smesh = smeshBuilder.New(theStudy)
#porousCellMesh = smesh.Mesh(porousCell)

#if(
    #checkFilletSwitcher != 0
#):
    #globalNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)

    #globalNETGEN1D2DParameters = globalNETGEN1D2D.Parameters()
    #globalNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    #globalNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    #globalNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    #globalNETGEN1D2DParameters.SetOptimize( 1 )
    #globalNETGEN1D2DParameters.SetFineness( 5 )
    #globalNETGEN1D2DParameters.SetGrowthRate( globalGrowthRate2D )
    #globalNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    #globalNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    #globalNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    #globalNETGEN1D2DParameters.SetSecondOrder( 0 )
    #globalNETGEN1D2DParameters.SetFuseEdges( 1 )

    ###Sub-mesh for leftSide boundary
    ##leftSideNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=leftSide)

    ##leftSideNETGEN1D2DParameters = leftSideNETGEN1D2D.Parameters()
    ##leftSideNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    ##leftSideNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    ##leftSideNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    ##leftSideNETGEN1D2DParameters.SetOptimize( 1 )
    ##leftSideNETGEN1D2DParameters.SetFineness( 5 )
    ##leftSideNETGEN1D2DParameters.SetGrowthRate( growthRate2D )
    ###leftSideNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    ###leftSideNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    ##leftSideNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    ##leftSideNETGEN1D2DParameters.SetSecondOrder( 0 )
    ##leftSideNETGEN1D2DParameters.SetFuseEdges( 1 )

    ###Sub-mesh for oppositeSide boundary
    ##oppositeSideNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=oppositeSide)

    ##oppositeSideNETGEN1D2DParameters = oppositeSideNETGEN1D2D.Parameters()
    ##oppositeSideNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    ##oppositeSideNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    ##oppositeSideNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    ##oppositeSideNETGEN1D2DParameters.SetOptimize( 1 )
    ##oppositeSideNETGEN1D2DParameters.SetFineness( 5 )
    ##oppositeSideNETGEN1D2DParameters.SetGrowthRate( growthRate2D )
    ###oppositeSideNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    ###oppositeSideNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    ##oppositeSideNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    ##oppositeSideNETGEN1D2DParameters.SetSecondOrder( 0 )
    ##oppositeSideNETGEN1D2DParameters.SetFuseEdges( 1 )

    ###Sub-mesh for topSide boundary
    ##topSideNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=topSide)

    ##topSideNETGEN1D2DParameters = topSideNETGEN1D2D.Parameters()
    ##topSideNETGEN1D2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    ##topSideNETGEN1D2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    ##topSideNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    ##topSideNETGEN1D2DParameters.SetOptimize( 1 )
    ##topSideNETGEN1D2DParameters.SetFineness( 5 )
    ##topSideNETGEN1D2DParameters.SetGrowthRate( growthRate2D )
    ###topSideNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
    ###topSideNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
    ##topSideNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
    ##topSideNETGEN1D2DParameters.SetSecondOrder( 0 )
    ##topSideNETGEN1D2DParameters.SetFuseEdges( 1 )
#else:
    #globalRegular1D = porousCellMesh.Segment()

    #globalLocalLength = globalRegular1D.LocalLength( globalNETGEN_2D_MinSize, None, 1e-6 )

    ##globalFixedPoints = globalRegular1D.FixedPoints1D([ 0.1, 0.9 ],[ 10, 50, 10 ],[])
    ##globalFixedPoints.SetObjectEntry( "porousCell" )

    ##globalDeflection = globalRegular1D.Deflection1D(4e-5)

    #globalNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D)

    #globalNETGEN2DParameters = globalNETGEN2D.Parameters()
    #globalNETGEN2DParameters.SetMaxSize( globalNETGEN_2D_MaxSize )
    #globalNETGEN2DParameters.SetMinSize( globalNETGEN_2D_MinSize )
    #globalNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
    #globalNETGEN2DParameters.SetOptimize( 1 )
    #globalNETGEN2DParameters.SetFineness( 5 )
    #globalNETGEN2DParameters.SetGrowthRate( globalGrowthRate2D )
    ##globalNETGEN2DParameters.SetNbSegPerEdge( 2 )
    ##globalNETGEN2DParameters.SetNbSegPerRadius( 3 )
    #globalNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
    #globalNETGEN2DParameters.SetSecondOrder( 0 )
    #globalNETGEN2DParameters.SetFuseEdges( 1 )

##Sub-mesh for leftSide boundary
#leftSideRegular1D = porousCellMesh.Segment(geom=leftSide)

#leftSideLocalLength = leftSideRegular1D.LocalLength(patchesNETGEN_2D_MinSize, None, 1e-6)
##leftSideDeflection = Regular_1D_1.Deflection1D(1e-5)

#leftSideNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=leftSide)

#leftSideNETGEN2DParameters = leftSideNETGEN2D.Parameters()
#leftSideNETGEN2DParameters.SetMaxSize( patchesNETGEN_2D_MaxSize )
#leftSideNETGEN2DParameters.SetMinSize( patchesNETGEN_2D_MinSize )
#leftSideNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
#leftSideNETGEN2DParameters.SetOptimize( 1 )
#leftSideNETGEN2DParameters.SetFineness( 5 )
#leftSideNETGEN2DParameters.SetGrowthRate( patchesGrowthRate2D )
#leftSideNETGEN2DParameters.SetNbSegPerEdge( 2 )
#leftSideNETGEN2DParameters.SetNbSegPerRadius( 3 )
#leftSideNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
#leftSideNETGEN2DParameters.SetSecondOrder( 0 )
#leftSideNETGEN2DParameters.SetFuseEdges( 1 )

##Sub-mesh for oppositeSide boundary
#oppositeSideRegular1D = porousCellMesh.Segment(geom=oppositeSide)

##oppositeSideDeflection = oppositeSideRegular1D.Deflection1D(1e-5)
#oppositeSideLocalLength = oppositeSideRegular1D.LocalLength(patchesNETGEN_2D_MinSize, None, 1e-6)

#oppositeSideNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=oppositeSide)

#oppositeSideNETGEN2DParameters = oppositeSideNETGEN2D.Parameters()
#oppositeSideNETGEN2DParameters.SetMaxSize( patchesNETGEN_2D_MaxSize )
#oppositeSideNETGEN2DParameters.SetMinSize( patchesNETGEN_2D_MinSize )
#oppositeSideNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
#oppositeSideNETGEN2DParameters.SetOptimize( 1 )
#oppositeSideNETGEN2DParameters.SetFineness( 5 )
#oppositeSideNETGEN2DParameters.SetGrowthRate( patchesGrowthRate2D )
#oppositeSideNETGEN2DParameters.SetNbSegPerEdge( 2 )
#oppositeSideNETGEN2DParameters.SetNbSegPerRadius( 3 )
#oppositeSideNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
#oppositeSideNETGEN2DParameters.SetSecondOrder( 0 )
#oppositeSideNETGEN2DParameters.SetFuseEdges( 1 )

##Sub-mesh for topSide boundary
#topSideRegular1D = porousCellMesh.Segment(geom=topSide)

##topSideDeflection = topSideRegular1D.Deflection1D(1e-5)
#topSideLocalLength = topSideRegular1D.LocalLength(patchesNETGEN_2D_MinSize, None, 1e-6)

#topSideNETGEN2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=topSide)

#topSideNETGEN2DParameters = topSideNETGEN2D.Parameters()
#topSideNETGEN2DParameters.SetMaxSize( patchesNETGEN_2D_MaxSize )
#topSideNETGEN2DParameters.SetMinSize( patchesNETGEN_2D_MinSize )
#topSideNETGEN2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
#topSideNETGEN2DParameters.SetOptimize( 1 )
#topSideNETGEN2DParameters.SetFineness( 5 )
#topSideNETGEN2DParameters.SetGrowthRate( patchesGrowthRate2D )
#topSideNETGEN2DParameters.SetNbSegPerEdge( 2 )
#topSideNETGEN2DParameters.SetNbSegPerRadius( 3 )
#topSideNETGEN2DParameters.SetUseSurfaceCurvature( 1 )
#topSideNETGEN2DParameters.SetSecondOrder( 0 )
#topSideNETGEN2DParameters.SetFuseEdges( 1 )

##Projection leftSide's sub-mesh on rightSide boundary
#rightSideProjection1D2D = porousCellMesh.Projection1D2D(geom=rightSide)
#leftSideSourceFace = rightSideProjection1D2D.SourceFace(leftSide,None,None,None,None,None)

##Projection oppositeSide's sub-mesh on frontSide boundary
#frontSideProjection1D2D = porousCellMesh.Projection1D2D(geom=frontSide)
#oppositeSideSourceFace = frontSideProjection1D2D.SourceFace(oppositeSide,None,None,None,None,None)

##Projection topSide's sub-mesh on bottomSide boundary
#bottomSideProjection1D2D = porousCellMesh.Projection1D2D(geom=bottomSide)
#topSideSourceFace = bottomSideProjection1D2D.SourceFace(topSide,None,None,None,None,None)

#throatsNETGEN1D2D = porousCellMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=throats)

#throatsNETGEN1D2DParameters = throatsNETGEN1D2D.Parameters()
#throatsNETGEN1D2DParameters.SetMaxSize( throatsNETGEN_2D_MaxSize )
#throatsNETGEN1D2DParameters.SetMinSize( throatsNETGEN_2D_MinSize )
#throatsNETGEN1D2DParameters.SetQuadAllowed( NETGEN_2D_QuadAllowed )
#throatsNETGEN1D2DParameters.SetOptimize( 1 )
#throatsNETGEN1D2DParameters.SetFineness( 5 )
#throatsNETGEN1D2DParameters.SetGrowthRate( throatsGrowthRate2D )
#throatsNETGEN1D2DParameters.SetNbSegPerEdge( 2 )
#throatsNETGEN1D2DParameters.SetNbSegPerRadius( 3 )
#throatsNETGEN1D2DParameters.SetUseSurfaceCurvature( 1 )
#throatsNETGEN1D2DParameters.SetSecondOrder( 0 )
#throatsNETGEN1D2DParameters.SetFuseEdges( 1 )

#gloabalNETGEN3D = porousCellMesh.Tetrahedron()

#globalNETGEN3DParameters = gloabalNETGEN3D.Parameters()
#globalNETGEN3DParameters.SetMaxSize( globalNETGEN_3D_MaxSize )
#globalNETGEN3DParameters.SetMinSize( globalNETGEN_3D_MinSize )
#globalNETGEN3DParameters.SetOptimize( 1 )
#globalNETGEN3DParameters.SetFineness( 5 )
#globalNETGEN3DParameters.SetGrowthRate( globalGrowthRate3D )
#globalNETGEN3DParameters.SetUseSurfaceCurvature( 0 )
#globalNETGEN3DParameters.SetSecondOrder( 0 )
#globalNETGEN3DParameters.SetFuseEdges( 1 )

#throatsNETGEN3D = porousCellMesh.Tetrahedron(algo=smeshBuilder.NETGEN_3D,geom=throats)

#throatsNETGEN3DParameters = throatsNETGEN3D.Parameters()
#throatsNETGEN3DParameters.SetMaxSize( throatsNETGEN_3D_MaxSize )
#throatsNETGEN3DParameters.SetMinSize( throatsNETGEN_3D_MinSize )
#throatsNETGEN3DParameters.SetOptimize( 1 )
#throatsNETGEN3DParameters.SetFineness( 5 )
#throatsNETGEN3DParameters.SetGrowthRate( throatsGrowthRate3D )
#throatsNETGEN3DParameters.SetUseSurfaceCurvature( 0 )
#throatsNETGEN3DParameters.SetSecondOrder( 0 )
#throatsNETGEN3DParameters.SetFuseEdges( 1 )

#if (
    #intersectionParameter_ != 0
#):
    ##viscousLayersSize2D = 0.025 / (VTranslationDirectionNbTimes_ * WTranslationDirectionNbTimes_) * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
    ##viscousLayersSize3D = 0.01 * geompy.BasicProperties(onePore)[2] ** (1.0 / 3.0)
    ##viscousLayersSize3D = 0.03 / (VTranslationDirectionNbTimes_ * WTranslationDirectionNbTimes_) * geompy.BasicProperties(leftSide)[1] ** (1.0 / 2.0)
    ##viscousLayersSize3D = 0.5 * localLengthSize
    ##Viscous_Layers_2D_1 = NETGEN_1D_2D_1.ViscousLayers2D(viscousLayersSize2D,viscousLayersNumber,viscousLayersGrowth)
    ##Viscous_Layers_2D_2 = NETGEN_1D_2D_2.ViscousLayers2D(viscousLayersSize2D,viscousLayersNumber,viscousLayersGrowth)
    ##Viscous_Layers_2D_3 = NETGEN_1D_2D_3.ViscousLayers2D(viscousLayersSize2D,viscousLayersNumber,viscousLayersGrowth)
    #globalViscousLayers3D = gloabalNETGEN3D.ViscousLayers(
        #globalViscousLayersSize3D,
        #globalViscousLayersNumber,
        #globalViscousLayersGrowth,
        #allFacesIDsWithoutViscousLayers,
        #1,
        #StdMeshersBuilder.FACE_OFFSET
    #)
    ##throatsViscousLayers3D = throatsNETGEN3D.ViscousLayers(
        ##throatsViscousLayersSize3D,
        ##throatsViscousLayersNumber,
        ##throatsViscousLayersGrowth,
        ##allFacesIDsWithoutViscousLayers,
        ##1,
        ##StdMeshersBuilder.FACE_OFFSET
    ##)
    ##Viscous_Layers_3D = NETGEN_3D.ViscousLayers(viscousLayersSize3D,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.SURF_OFFSET_SMOOTH)
    ##Viscous_Layers_3D = NETGEN_3D.ViscousLayers(viscousLayersSize3D,viscousLayersNumber,viscousLayersGrowth,allFacesIDsWithoutViscousLayers,1,StdMeshersBuilder.NODE_OFFSET)

##Compute the mesh
#isDone = porousCellMesh.Compute()

#skeletonWallBoundary = porousCellMesh.GroupOnGeom(skeletonWall,'skeletonWall',SMESH.FACE)
#leftSideBoundary = porousCellMesh.GroupOnGeom(leftSide,'leftSide',SMESH.FACE)
#rightSideBoundary = porousCellMesh.GroupOnGeom(rightSide,'rightSide',SMESH.FACE)
#oppositeSideBoundary = porousCellMesh.GroupOnGeom(oppositeSide,'oppositeSide',SMESH.FACE)
#frontSideBoundary = porousCellMesh.GroupOnGeom(frontSide,'frontSide',SMESH.FACE)
#topSideBoundary = porousCellMesh.GroupOnGeom(topSide,'topSide',SMESH.FACE)
#bottomSideBoundary = porousCellMesh.GroupOnGeom(bottomSide,'bottomSide',SMESH.FACE)

### Get viscous layers
##aCriteria = []
##aCriterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,skeletonWall,SMESH.FT_Undefined,SMESH.FT_Undefined,0.9 * viscousLayersSize3D)
##aCriteria.append(aCriterion)
##aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
##aFilter_1.SetMesh(porousCellMesh.GetMesh())
##viscousLayers = porousCellMesh.GroupOnFilter( SMESH.VOLUME, 'viscousLayers', aFilter_1 )

### Get film
##aCriteria = []
##aCriterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,skeletonWall,SMESH.FT_Undefined,SMESH.FT_Undefined,1e-7)
##aCriteria.append(aCriterion)
##aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
##aFilter_1.SetMesh(porousCellMesh.GetMesh())
##film = porousCellMesh.GroupOnFilter( SMESH.VOLUME, 'film', aFilter_1 )

#aCriteria = []
#aCriterion = smesh.GetCriterion(SMESH.VOLUME,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,skeletonWall,SMESH.FT_Undefined,SMESH.FT_Undefined,1.1 * globalViscousLayersSize3D)
#aCriteria.append(aCriterion)
#aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
#aFilter_1.SetMesh(porousCellMesh.GetMesh())
#film = porousCellMesh.GroupOnFilter( SMESH.VOLUME, 'film', aFilter_1 )

#porousCellMesh.ExportUNV( absolute_case_path_ + "/porousCell.unv" )

#if salome.sg.hasDesktop():
    #salome.sg.updateObjBrowser(1)

###Save geometry and mesh to file
#print("Save geometry and mesh to " + absolute_case_path_ + "/porousCellGeometryAndMesh.hdf\n")

#salome.myStudyManager.SaveAs(absolute_case_path_ + "/porousCellGeometryAndMesh.hdf", salome.myStudy, False)

#print("done")
#print("***********************************\n")
'''
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
