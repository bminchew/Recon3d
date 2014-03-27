#!/usr/bin/python

"""
plot two num files of the same shape

Usage:  plot_num.py num1 num1 cols

"""
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def main(args):
   fid = open(args[0],'r')
   data1 = np.fromfile(fid,dtype=np.float32).reshape(-1,np.int32(args[2]))
   fid.close()
   fid = open(args[1],'r')
   data2 = np.fromfile(fid,dtype=np.float32).reshape(-1,np.int32(args[2]))
   fid.close()

   fig = plt.figure()
   surf = plt.imshow(data1,cmap=cm.jet)
   plt.clim([0,7])
   cbar = fig.colorbar(surf)

   fig = plt.figure()
   surf = plt.imshow(data2,cmap=cm.jet)
   plt.clim([0,7])
   cbar = fig.colorbar(surf)

   fig = plt.figure()
   surf = plt.imshow(data1-data2,cmap=cm.jet)
   plt.clim([-3,7])
   cbar = fig.colorbar(surf)

   plt.show() 

if __name__=='__main__':
   args = sys.argv[1:]
   if len(args) != 3:
      print(__doc__)
      sys.exit()
   main(args)
