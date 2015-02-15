#!/usr/bin/env python
import numpy as np
from numpy.linalg import norm

def Residuals(Z): return np.sum(np.dot(Z,Z.transpose()))
def OrthoTest(Z): return np.sum(np.dot(Z,Z.transpose())-np.eye(Z.shape[0],Z.shape[1]))
""" Point of this code is a unittest of the gather function which is becoming customized due ot limitations on buffer space """
collec_baseline_f = open('collec_baseline','r')
gather_baseline_f = open('gather_baseline','r')
collect_A_baseline = np.genfromtxt(collec_baseline_f)
gather_A_baseline = np.genfromtxt(gather_baseline_f)

print OrthoTest(gather_A_baseline)
print OrthoTest(collect_A_baseline)
print Residuals(collect_A_baseline - gather_A_baseline)

# gatherings
gather_n2_f = [ open('gather%03d_n2.txt' % (i),'r') for i in range(2) ]
gather_n4_f = [ open('gather%03d_n4.txt' % (i),'r') for i in range(4) ]
gather_A_n2 = [np.genfromtxt(n2) for n2 in gather_n2_f]
gather_A_n4 = [np.genfromtxt(n4) for n4 in gather_n4_f]
gather_n2 = np.sum([n2 for n2 in gather_A_n2],axis=0)
gather_n4 = np.sum([n4 for n4 in gather_A_n4],axis=0)

# collections
collect_n2_f = [ open('collec%03d_n2.txt' % (i),'r') for i in range(2) ]
collect_n4_f = [ open('collec%03d_n4.txt' % (i),'r') for i in range(4) ]
collect_A_n2 = [np.genfromtxt(n2) for n2 in collect_n2_f]
collect_A_n4 = [np.genfromtxt(n4) for n4 in collect_n4_f]
collect_res_n2 = [Residuals(gather_n2-collect_n2) for collect_n2 in collect_A_n2]
collect_res_n4 = [Residuals(gather_n4-collect_n4) for collect_n4 in collect_A_n4]

print 'Error panel:'
print 'gather - ortho test for eigenvectors '
print OrthoTest(gather_n2)
print OrthoTest(gather_n4)

print 'collect'
print collect_res_n2
print (collect_A_n2[0]-collect_A_n2[1])[0]
print (collect_A_n2[0]-collect_A_n2[1])[33]
print collect_res_n4


