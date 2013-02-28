from heliosonic import *
import sys

sfft = SpatialFFT('/home/chrono/data/blork/helio/spatialFFT.npy')
if len(sys.argv)>1:
   sfft.DisplayRing( int(sys.argv[1]) )
else:
   sfft.ExploreWAV()
