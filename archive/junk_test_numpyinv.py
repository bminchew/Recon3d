#!/usr/bin/env python

import time
import numpy as np

a = np.random.rand(1.e3,1.e3)
t0 = time.time()
np.linalg.inv(a)
time10_3 = time.time() - t0

a = np.random.rand(1.e4,1.e4)
t0 = time.time()
np.linalg.inv(a)
time10_4 = time.time() - t0


print time10_3,time10_4

tslop = (time10_4-time10_3)/(1.e8-1.e6)

print 'estimated 10_6 time =',tslop*1.e12

