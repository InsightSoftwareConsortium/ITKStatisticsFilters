#!/usr/bin/env python

import SimpleITK as sitk
from scipy import stats
import sys
import numpy as np

if len(sys.argv) < 5 or len(sys.argv) > 7:
  print "Usage: " + sys.argv[0] + "inFile refFile radius outFile [bias] [fisher]"
  sys.exit(1)

bias=False
if len(sys.argv)>=6 and sys.argv[5] == "bias":
  bias=True
fisher=False
if len(sys.argv)==7 and sys.argv[6] == "fisher":
  fisher=True
if bias:
  print "Bias correction ON"
else:
  print "Bias correction OFF"
if fisher:
  print "Fisher correction ON"
else:
  print "Fisher correction OFF"

inFile=sys.argv[1]
refFile=sys.argv[2]
radius=int(sys.argv[3])
outFile=sys.argv[4]
image = sitk.ReadImage(inFile)
arr= sitk.GetArrayFromImage(image)
shape=arr.shape
arr_kurtosis=np.zeros(shape)

ref_image = sitk.ReadImage(refFile)
ref_arr = sitk.GetArrayFromImage(ref_image)

print "Shape:"+str(shape)
for x in range(radius,shape[0]-radius-1):
  for y in range(radius,shape[1]-radius-1):
    if len(shape) == 2:
      sub_arr=arr[x-radius:x+radius+1, y-radius:y+radius+1]
      line=sub_arr.reshape([(radius*2+1)**2,1])
      arr_kurtosis[x,y]=stats.kurtosis(line,bias=bias, fisher=fisher)
    else:
      for z in range(radius,shape[2]-radius-1):
        sub_arr=arr[x-radius:x+radius+1, y-radius:y+radius+1, z-radius:z+radius+1]
        line=sub_arr.reshape([(radius*2+1)**3,1])
        arr_kurtosis[x,y,z]=stats.kurtosis(line,bias=bias, fisher=fisher)
print "ref_arr shape:"+str(ref_arr.shape)
print "arr_kurtosis shape:"+str(arr_kurtosis.shape)
if len(shape) == 2:
  if (ref_arr[radius:-radius,radius:-radius]
  == arr_kurtosis[radius:-radius,radius:-radius]).all():
    print "Test passed"
    sys.exit(0)
elif (ref_arr[radius:-radius,radius:-radius,radius:-radius]
     == arr_kurtosis[radius:-radius,radius:-radius,radius:-radius]).all():
  print "Test passed"
  sys.exit(0)

outImage = sitk.GetImageFromArray(arr_kurtosis)
outImage.SetOrigin(image.GetOrigin())
outImage.SetSpacing(image.GetSpacing())
outImage.SetDirection(image.GetDirection())
print "Test failed"
print "Saving numpy kurtosis image as %s"%outFile
sitk.WriteImage(outImage, outFile,True)
