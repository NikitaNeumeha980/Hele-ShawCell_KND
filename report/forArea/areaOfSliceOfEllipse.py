import matplotlib.pyplot as plt
# from scipy import integrate
import numpy as np
import math

####################################################
##            Begin of variables section          ##
####################################################
# absoluteCasePath_ = sys.argv[1]
drainageLength_ = 5000
drainageWidth_ = 200
channel1InOutWidth_ = 100
channel1MiddleWidth_ = 50 # upside
channel2InOutWidth_ = 100
channel2MiddleWidth_ = 140 #downside
depth_ = 1000
blobSeparatingDistance_ = 0
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
    




####    area    ####
print("a1 , b1= ", a1_," ", b1_)
print("a2 , b2= ", a2_," ", b2_)
# print("second point = ", np.sqrt(4 * blobRadius_**2 + 8 * channel1InOutWidth_* blobRadius_ + 4 * channel1InOutWidth_ - drainageWidth_**2))
# print("arcsin = ", np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2) / (2 * a1_))
# print("side = ", np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2))
# print("r+w1 = ", blobRadius_ + channel1InOutWidth_)
# print("1/2w = ", drainageWidth_ / 2)

S_ellips1_phi = 0.5 * np.pi * a1_ * b1_ - \
    0.5 * np.pi * blobRadius_**2 - \
    (np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2) * drainageWidth_ / 2) / 2 + b1_/(a1_*2) * \
    (   
        (a1_**2 * np.pi / 2) - \
        (np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2) * \
        np.sqrt(a1_**2 - np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2)) + 
        a1_**2 * np.arcsin(np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2) / a1_))
    ) - \
    (blobRadius_**2 * 0.5 * np.arccos((a1_/b1_ * np.sqrt(b1_**2 - drainageWidth_**2 / 4)) / (channel1InOutWidth_ + blobRadius_)))

print("S1 = ",S_ellips1_phi)

S_ellips2_phi = 0.5 * np.pi * a2_ * b2_ - \
    0.5 * np.pi * blobRadius_**2 - \
    (np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2) * drainageWidth_ / 2) / 2 + b2_/(a2_*2) * \
    (   
        (a2_**2 * np.pi / 2) - \
        (np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2) * \
        np.sqrt(a2_**2 - np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2)) + 
        a2_**2 * np.arcsin(np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2) / a2_))
    ) - \
    (blobRadius_**2 * 0.5 * np.arccos(round((a2_/b2_ * np.sqrt(b2_**2 - drainageWidth_**2 / 4)) / (channel2InOutWidth_ + blobRadius_), 1)))

print("S1 = ",S_ellips2_phi)

print("s1/s2= ", S_ellips1_phi/S_ellips2_phi )


print((a2_/b2_ * np.sqrt(b2_**2 - drainageWidth_**2 / 4)) / (channel2InOutWidth_ + blobRadius_))


# b1FromS1 = S_ellips1_phi * 2 * a1_ / (
#     (a1_**2 * np.pi / 2) - 
#     (np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2) * 
#      np.sqrt(a1_**2 - (np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2))**2) + 
#      a1_**2 * np.arcsin(np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ / 2)**2) / a1_))
# )

# print("b1FromS1 = ", b1FromS1)


# print("S1 = ",S_ellips1_phi)

# S_ellips2 = (b2_/(a2_*2)) * (
#     (a2_**2 * np.pi / 2) - 
#     (np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2) * 
#      np.sqrt(a2_**2 - (np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2))) + 
#      a2_**2 * np.arcsin(np.sqrt((blobRadius_ + channel2InOutWidth_)**2 - (drainageWidth_ / 2)**2) / a2_))
# )


# print("S2 = ",S_ellips2)

# print("arc = ", np.arccos((a1_/b1_ * np.sqrt(b1_**2 - drainageWidth_**2 / 4)) / (channel1InOutWidth_ + blobRadius_)) * 180 / np.pi)


# S_ellips1 = b1_/a1_ * integrate(x, b1_, math.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ * 0.5)**2)):

# def S_ellips1(x):
#     return (np.sqrt(b1_**2 - x**2))

# AreaEllips1 = integrate.quad(S_ellips1, b1_, np.sqrt((blobRadius_ + channel1InOutWidth_)**2 - (drainageWidth_ * 0.5)**2))

# print(AreaEllips1)
 
    
    
    
    
    
