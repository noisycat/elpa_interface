#!/usr/bin/env python
import numpy as np
from numpy.linalg import norm
""" Point of this code is a unittest of the gather function which is becoming customized due ot limitations on buffer space """
collec_baseline_f = open('collec_baseline','r')
collect_A_baseline = -np.genfromtxt(collec_baseline_f)
gather_baseline_f = open('gather_baseline','r')
gather_A_baseline = -np.genfromtxt(gather_baseline_f)
# gatherings
gather_n2_f = [ open('gather%03d_n2.txt' % (i),'r') for i in range(2) ]
gather_n4_f = [ open('gather%03d_n4.txt' % (i),'r') for i in range(4) ]
gather_A_n2 = [np.genfromtxt(n2) for n2 in gather_n2_f]
gather_A_n4 = [np.genfromtxt(n4) for n4 in gather_n4_f]
gather_res_n2 = min(norm(gather_A_baseline - sum(gather_A_n2)),norm(gather_A_baseline + sum(gather_A_n2)))
gather_res_n4 = min(norm(gather_A_baseline - sum(gather_A_n4)),norm(gather_A_baseline + sum(gather_A_n4)))

# collections
collect_n2_f = [ open('collec%03d_n2.txt' % (i),'r') for i in range(2) ]
collect_n4_f = [ open('collec%03d_n4.txt' % (i),'r') for i in range(4) ]
collect_A_n2 = [np.genfromtxt(n2) for n2 in collect_n2_f]
collect_A_n4 = [np.genfromtxt(n4) for n4 in collect_n4_f]
collect_res_n2 = [min(norm(collect_A_baseline - A_n2_i),norm(collect_A_baseline+A_n2_i)) for A_n2_i in collect_A_n2]
collect_res_n4 = [min(norm(collect_A_baseline - A_n4_i),norm(collect_A_baseline+A_n4_i)) for A_n4_i in collect_A_n4]

print 'Error panel:'
print 'gather'
print gather_res_n2
print gather_res_n4

print 'collect'
print collect_res_n2
print collect_res_n4

