#!/usr/bin/python

import numpy as np


a = np.arange(25,dtype=np.float32).reshape(-1,5)
print a

fid = open('junk_testio.dat','w')
a.flatten().tofile(fid)
fid.close()

b = np.empty_like(a)
fid = open('junk_testio.dat','r')
fid.seek(8*np.shape(a)[0])
for i in np.arange(4):
   b = np.fromfile(fid,dtype=np.float32,count=5)
   print i
   print b


