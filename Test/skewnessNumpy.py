#!/usr/bin/env python

import SimpleITK as sitk
from scipy import stats
import sys
import numpy as np

if len(sys.argv) < 5 or len(sys.argv) > 6:
  print "Usage: " + sys.argv[0] + "inFile refFile radius outFile [bias]"
  sys.exit(1)

bias=False
if len(sys.argv)==6:
  bias=True
if bias:
  print "Bias correction ON"
else:
  print "Bias correction OFF"
inFile=sys.argv[1]
refFile=sys.argv[2]
radius=int(sys.argv[3])
outFile=sys.argv[4]
image = sitk.ReadImage(inFile)
arr= sitk.GetArrayFromImage(image)
shape=arr.shape
arr_skew=np.zeros(shape)

ref_image = sitk.ReadImage(refFile)
ref_arr = sitk.GetArrayFromImage(ref_image)

print "Shape:"+str(shape)
for x in range(radius,shape[0]-radius-1):
  for y in range(radius,shape[1]-radius-1):
    if len(shape) == 2:
      sub_arr=arr[x-radius:x+radius+1, y-radius:y+radius+1]
      line=sub_arr.reshape([(radius*2+1)**2,1])
      arr_skew[x,y]=stats.skew(line,bias=bias)
    else:
      for z in range(radius,shape[2]-radius-1):
        sub_arr=arr[x-radius:x+radius+1, y-radius:y+radius+1, z-radius:z+radius+1]
        line=sub_arr.reshape([(radius*2+1)**3,1])
        arr_skew[x,y,z]=stats.skew(line,bias=bias)
print "ref_arr shape:"+str(ref_arr.shape)
print "arr_skew shape:"+str(arr_skew.shape)
if len(shape) == 2:
  if (ref_arr[radius:-radius,radius:-radius]
  == arr_skew[radius:-radius,radius:-radius]).all():
    print "Test passed"
    sys.exit(0)
elif (ref_arr[radius:-radius,radius:-radius,radius:-radius]
     == arr_skew[radius:-radius,radius:-radius,radius:-radius]).all():
  print "Test passed"
  sys.exit(0)

outImage = sitk.GetImageFromArray(arr_skew)
outImage.SetOrigin(image.GetOrigin())
outImage.SetSpacing(image.GetSpacing())
outImage.SetDirection(image.GetDirection())
print "Test failed"
print "Saving numpy skewness image as %s"%outFile
sitk.WriteImage(outImage, outFile,True)
