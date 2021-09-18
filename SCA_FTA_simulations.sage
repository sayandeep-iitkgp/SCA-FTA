#!/usr/bin/env sage

"""
Attack Simulations for ASIACRYPT 2021 paper. 
Author: Sayandeep Saha
"""

import os
import sys
import shutil
import time
import re
import random as pyrandom
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import skew
from scipy.stats import norm

from sage_imports import *

################### Supporting Functions ###########################

def pause():
	programPause = raw_input("Press the <ENTER> key to continue...")

def stop():
	sys.exit("Execution Stopped by the User...")

def bitarrtoint(X, len=4):
	X_int = 0	
	for i in xrange(len):
		if (X[i] == 1):
			X_int = X_int + 2^((len-1)-i)
	return X_int

def inttobitarr(X, len=4):
	l = []
	for i in xrange(len):
		l.append(1 if 2^i & X else 0)
	l = map(GF(2),l)
	X = list(reversed(l))
	return X

def compute_log_likelihood(x, mu, sd):
	ll = 0

	LL = np.sum(norm.logpdf(x, mu, sd))
    

	ll = LL
	return ll


############################# Target SBoxes #######################			
							

# SIFA Protected chi3 S-Box
def Clone (v, v_):
	v_ = v
	return v_

def lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = None):
	#print("here")

	r0 = None
	r1 = None 
	s0 = None
	s1 = None 
	t0 = None 
	t1 = None
	
	
	Rr_ = None
	Rr_ = Clone(Rr, Rr_)
	
	Rt_ = None
	Rt_ = Clone(Rt, Rt_)
	
	Rs = Rr_ + Rt_		#
	
	if (fault_loc == 0):
		c0 = c0 + 1 		# fault
	elif (fault_loc == 1):
		b0 = b0 + 1			# fault
	elif (fault_loc == 2):
		a1 = a1 + 1			# fault

	
	b0_ = None
	b0_ = Clone (b0, b0_)

	c1_ = None
	c1_ = Clone (c1, c1_)	
		
	T0 = (b0_+ 1) * c1_		#
	
	
	a1_ = None
	a1_ = Clone (a1, a1_)	

	b1_ = None
	b1_ = Clone (b1, b1_)		
	 
	T2 = a1_ * b1_		# 
	 
	
	b0_ = None
	b0_ = Clone (b0, b0_)	
	 
	c0_ = None
	c0_ = Clone (c0, c0_)
	
	T1 = (b0_ + 1) * c0_		#
	
	
	a1_ = None
	a1_ = Clone (a1, a1_)	
	 
	b0_ = None
	b0_ = Clone (b0, b0_)
	
	T3 = a1_ * b0_		#	
	
	
	
	Rr_ = None
	Rr_ = Clone(Rr, Rr_)	
	
	r0 = T0 + Rr_		#
	

	Rt_ = None
	Rt_ = Clone(Rt, Rt_)	
	
	t1 = T2 + Rt_		#	


	r0 = r0 + T1		#
	
	t1 = t1 + T3		#
	
	
	
	
	c0_ = None
	c0_ = Clone (c0, c0_)
	
	a1_ = None
	a1_ = Clone (a1, a1_)	
	
	T0 = (c0_ + 1)*a1_	#


	b1_ = None
	b1_ = Clone (b1, b1_)
	
	c1_ = None
	c1_ = Clone (c1, c1_)	
	
	T2 = b1_ * c1_		#
	

	c0_ = None
	c0_ = Clone (c0, c0_)
	
	a0_ = None
	a0_ = Clone (a0, a0_)	
	
	T1 = (c0_ + 1)*a0_	#
	

	b1_ = None
	b1_ = Clone (b1, b1_)
	
	c0_ = None
	c0_ = Clone (c0, c0_)	
	
	T3 = b1_ * c0_		#



	Rs_ = None
	Rs_ = Clone(Rs, Rs_)
	
	s0 = T0 + Rs_		#
	
	r1 = T2 + Rr		#
	
	s0 = s0 + T1		#
	
	r1 = r1 + T3		#
	


	a0_ = None
	a0_ = Clone (a0, a0_)
	
	b1_ = None
	b1_ = Clone (b1, b1_)
	
	T0 = (a0_ + 1)*b1_	#	
	
	c1_ = None
	c1_ = Clone (c1, c1_)
	
	a1_ = None
	a1_ = Clone (a1, a1_)
	
	T2 = c1_ * a1_		#
	

	a0_ = None
	a0_ = Clone (a0, a0_)
	
	b0_ = None
	b0_ = Clone (b0, b0_)
	
	T1 = (a0_ + 1)*b0_	#		
	 

	c1_ = None
	c1_ = Clone (c1, c1_)
	
	a0_ = None
	a0_ = Clone (a0, a0_)
	
	T3 = c1_ * a0_		#
	
	
	

	t0 = T0 + Rt		#
	
	s1 = T2 + Rs		#
	
	t0 = t0 + T1		#
	
	s1 = s1 + T3		#
	
	r0 = r0 + a0		#
	
	t1 = t1 + c1		#
	
	s0 = s0 + b0		#
	
	r1 = r1 + a1		#
	
	t0 = t0 + c0		#
	
	s1 = s1 + b1		#

	r = r0 + r1
	s = s0 + s1
	t = t0 + t1
	
	a = a0 + a1
	b = b0 + b1
	c = c0 + c1
	
	#print((a, b, c), (r, s, t))
	#print(bitarrtoint([a, b, c], len=3), bitarrtoint([r, s, t], len=3))

	return r0, r1, s0, s1, t0, t1
	

# Unmasked PRESENT S-Box
sbox_tab = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2], big_endian=True)
##################################################################


#------------------------------------------------------------------------------
# Testing Unmasked PRESENT S-Box (Noise-free; error-detection on unshared values)
#------------------------------------------------------------------------------

diff_dict = {}
diff_dict_hw = {}
diff_dict_power = {}

f_val_list = [8, 4, 2, 1]
num_measurements = 10


for x in range(16):
	op_diff_list = []
	op_diff_list_hw = []
	for f_val in f_val_list:
		f = f_val		
		
		X = inttobitarr(x, 4)
		F = inttobitarr(f, 4)
		X_f = [0]*4
		for h in range(4):
			X_f[h] = X[h] + F[h]
		x_f = bitarrtoint(X_f)
		
		avg_hw = 0
		target_bit_val = -1
		for n in range(num_measurements):
			corr = sbox_tab(x)
			faulty = sbox_tab(x_f)
	
			# calculate output differential
			C = inttobitarr(corr, 4)
			F = inttobitarr(faulty, 4)
			Op_diff = [0]*4
			for h in range(4):
				Op_diff[h] = C[h] + F[h]
			op_diff = bitarrtoint(Op_diff)	
			#print(Op_diff)		
		
			# calculate the HW of the output differential
			op_diff_hw = 0
			for h in Op_diff:
				if (h == 1):
					op_diff_hw = op_diff_hw + 1	
				
			avg_hw = avg_hw + op_diff_hw
		avg_hw = float(avg_hw)/float(num_measurements)
		op_diff_list.append(op_diff)
		op_diff_list_hw.append(avg_hw)
	diff_dict[x] = tuple(op_diff_list)
	diff_dict_hw[x] = tuple(op_diff_list_hw)

# Constuct the fault template
key_classes = {}
for k, v in diff_dict_hw.iteritems():	
	key_classes[v] = key_classes.get(v, [])
	key_classes[v].append(k)

print("")
print("The template...")
print("")
print(key_classes)
print("")
print("")



# Attack (Template match)
x = 8                 # Change this value in the range (0, 15); this is the unknown in the online phase. 
fault_template = []
cnt = 0
for f in f_val_list:

	X = inttobitarr(x, 4)
	F = inttobitarr(f, 4)
	X_f = [0]*4
	for h in range(4):
		X_f[h] = X[h] + F[h]
	x_f = bitarrtoint(X_f)
	avg_hw = 0
	for n in range(num_measurements):
		corr = sbox_tab(x)
		faulty = sbox_tab(x_f)	
		
		# calculate output differential
		C = inttobitarr(corr, 4)
		F = inttobitarr(faulty, 4)
		Op_diff = [0]*4
		for h in range(4):
			Op_diff[h] = C[h] + F[h]
		op_diff = bitarrtoint(Op_diff)	
		#print(Op_diff)					
		
		# calculate HW of output differential
		op_diff_hw = 0
		for h in Op_diff:
			if (h == 1):
				op_diff_hw = op_diff_hw + 1		
		avg_hw = avg_hw + op_diff_hw
	avg_hw = float(avg_hw)/float(num_measurements)	
	
		
	fault_template.append(avg_hw)	

print ("Pattern constructed during template matching")	
print(fault_template)
print("")
print("Recovered value...")									
print(key_classes[tuple(fault_template)])	
print("")
print("")
#-----------------------------------------------------------------------


#stop()




#--------------------------------------------------------------------------------------------------------
# SIFA countermeasure with noisey HW simulatiom and Noisy faults (detection/correction at each share)
#--------------------------------------------------------------------------------------------------------


input_fault_patt_dict_1 = {}
input_fault_patt_dict_hw_1 = {}
num_measurements = 500 
num_avg_compute = 30  
noise_mean = 33.0
noise_std = 3.0
fault_sig_prob = 0.6



# Build the template		
for x in range(8):
	hw_set_patt = []
	fault_patt = [0]*3	
	for fault_loc in range(3):
		avg_traces_for_template_building = []
		for n in range(num_measurements):
			hw_set = []
			for nav in range(num_avg_compute):
				m1 = pyrandom.randint(0,7)
				rr = pyrandom.randint(0,1)
				rt = pyrandom.randint(0,1)				
					
				M1 = inttobitarr(m1, len=3)	
				X = inttobitarr(x, len=3)
				X1 = [X[i] + M1[i] for i in range(3)]
				X2 = M1
				Rr = inttobitarr(rr, len=1)
				Rt = inttobitarr(rt, len=1)

						 
				a0 = X1[0]
				a1 = X2[0]
				
				b0 = X1[1]
				b1 = X2[1]
				
				c0 = X1[2]
				c1 = X2[2]
				
				Rr = Rr[0]
				Rt = Rt[0]
				
				# Make some noise (fault noise)
				rand_fault_loc = pyrandom.randint(0,3) # choose a fault location randomly
				rn = sage.misc.prandom.random()
				if (rn < fault_sig_prob):		#signal
					r0, r1, s0, s1, t0, t1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = fault_loc)
				else:							#noise
					r0, r1, s0, s1, t0, t1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = rand_fault_loc)	
				
				#r0, r1, s0, s1, t0, t1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = fault_loc)
				r0_c, r1_c, s0_c, s1_c, t0_c, t1_c = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt)
				
				####### Bitwise error detection (uncomment this block if you test error-detection on shares. Also comment the error-correction block) ########
				#dr0 = r0 + r0_c
				#dr1 = r1 + r1_c
				
				#ds0 = s0 + s0_c
				#ds1 = s1 + s1_c
				
				#dt0 = t0 + t0_c
				#dt1 = t1 + t1_c
				
				#DS_inter = [dt1, dt0, ds1, ds0, dr1, dr0]		
				#########################################					

				####### Error Correction (uncomment this block if you test error-correction on shares. Also comment the error-detection block) ############
				
				r0_c1, r1_c1, s0_c1, s1_c1, t0_c1, t1_c1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt)
				
				
				ar0 = r0*r0_c 
				br0 = r0_c*r0_c1
				cr0 = r0*r0_c1
				
				r0_ec = ar0 + br0 + cr0
				
				ar1 = r1*r1_c 
				br1 = r1_c*r1_c1
				cr1 = r1*r1_c1
				
				r1_ec = ar1 + br1 + cr1		
				

				as0 = s0*s0_c 
				bs0 = s0_c*s0_c1
				cs0 = s0*s0_c1
				
				s0_ec = as0 + bs0 + cs0
				
				as1 = s1*s1_c 
				bs1 = s1_c*s1_c1
				cs1 = s1*s1_c1
				
				s1_ec = as1 + bs1 + cs1		
				
				
				at0 = t0*t0_c 
				bt0 = t0_c*t0_c1
				ct0 = t0*t0_c1
				
				t0_ec = at0 + bt0 + ct0
				
				at1 = t1*t1_c 
				bt1 = t1_c*t1_c1
				ct1 = t1*t1_c1
				
				t1_ec = at1 + bt1 + ct1	
				
				DS_inter = [at1, bt1, ct1, at0, bt0, ct0, as1, bs1, cs1, as0, bs0, cs0, ar1, br1, cr1, ar0, br0, cr0]
				
				#print(x)
				#print(fault_loc)
				#print(DS_inter)											
				#pause()
				
				######################################
				hw = 0
				for b in DS_inter:
					if (b == 1):
						hw = hw + 1
						
				# Add noise in HWs (power noise)--------
				hw = hw + np.random.normal(noise_mean, noise_std, 1)
				#pause()	
				hw_set.append(hw)
			#print(hw_set)
			#pause()	
			avg_hw_per_set = np.mean(hw_set)
			avg_traces_for_template_building.append(avg_hw_per_set)
		#print(avg_traces_for_template_building)	
		#pause()
		avg_hw_per_loc = np.mean(avg_traces_for_template_building)
		std_hw_per_loc = np.std(avg_traces_for_template_building)
		#print(avg_hw_per_loc)
		fault_patt[fault_loc] = (avg_hw_per_loc, std_hw_per_loc)	
	print(x)	
	input_fault_patt_dict_hw_1[x] = tuple(fault_patt)

# Constuct the fault template
key_classes_1 = {}
for k, v in input_fault_patt_dict_hw_1.iteritems():	
	key_classes_1[v] = key_classes_1.get(v, [])
	key_classes_1[v].append(k)

print("")
print("The template...")
print("")
print(key_classes_1)
print("")
print("")
print("")
#pause()




# Template matching

x = 4						# this is the unknown value to be exposed; change this in the range (0,7) to test template matching
num_measurements = 31
num_avg_compute = 35


# Euclidean Distance Based
fault_patt = [0]*3
for fault_loc in range(3):
	avg_traces_for_template_matching = []
	for n in range(num_measurements):
		hw_set = []
		for nav in range(num_avg_compute):
			m1 = pyrandom.randint(0,7)
			rr = pyrandom.randint(0,1)
			rt = pyrandom.randint(0,1)				
				
			M1 = inttobitarr(m1, len=3)	
			X = inttobitarr(x, len=3)
			X1 = [X[i] + M1[i] for i in range(3)]
			X2 = M1
			Rr = inttobitarr(rr, len=1)
			Rt = inttobitarr(rt, len=1)

					 
			a0 = X1[0]
			a1 = X2[0]
			
			b0 = X1[1]
			b1 = X2[1]
			
			c0 = X1[2]
			c1 = X2[2]
			
			Rr = Rr[0]
			Rt = Rt[0]
			
			
			# Make some noise (fault noise)
			rand_fault_loc = pyrandom.randint(0,3) # choose a fault location randomly
			rn = sage.misc.prandom.random()
			if (rn < fault_sig_prob):		#signal
				r0, r1, s0, s1, t0, t1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = fault_loc)
			else:							#noise
				r0, r1, s0, s1, t0, t1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = rand_fault_loc)	
			#r0, r1, s0, s1, t0, t1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt, fault_loc = fault_loc)
			r0_c, r1_c, s0_c, s1_c, t0_c, t1_c = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt)
			
			########## Error Detection (uncomment this block if you test error-detection on shares. Also comment the error-correction block) ##############
			#dr0 = r0 + r0_c
			#dr1 = r1 + r1_c
			
			#ds0 = s0 + s0_c
			#ds1 = s1 + s1_c
			
			#dt0 = t0 + t0_c
			#dt1 = t1 + t1_c
			
			#DS_inter = [dt1, dt0, ds1, ds0, dr1, dr0]
			#########################################
			
			####### Error Correction (uncomment this block if you test error-correction on shares. Also comment the error-detection block) ############
			
			r0_c1, r1_c1, s0_c1, s1_c1, t0_c1, t1_c1 = lam_S_box (a0, a1, b0, b1, c0, c1, Rr, Rt)
			
			ar0 = r0*r0_c 
			br0 = r0_c*r0_c1
			cr0 = r0*r0_c1
			
			r0_ec = ar0 + br0 + cr0
			
			ar1 = r1*r1_c 
			br1 = r1_c*r1_c1
			cr1 = r1*r1_c1
			
			r1_ec = ar1 + br1 + cr1		
			

			as0 = s0*s0_c 
			bs0 = s0_c*s0_c1
			cs0 = s0*s0_c1
			
			s0_ec = as0 + bs0 + cs0
			
			as1 = s1*s1_c 
			bs1 = s1_c*s1_c1
			cs1 = s1*s1_c1
			
			s1_ec = as1 + bs1 + cs1		
			
			
			at0 = t0*t0_c 
			bt0 = t0_c*t0_c1
			ct0 = t0*t0_c1
			
			t0_ec = at0 + bt0 + ct0
			
			at1 = t1*t1_c 
			bt1 = t1_c*t1_c1
			ct1 = t1*t1_c1
			
			t1_ec = at1 + bt1 + ct1	
			
			DS_inter = [at1, bt1, ct1, at0, bt0, ct0, as1, bs1, cs1, as0, bs0, cs0, ar1, br1, cr1, ar0, br0, cr0]
			
			#print(x)
			#print(fault_loc)
			#print(DS_inter)											
			#pause()
			
			######################################			
			
			hw = 0
			for b in DS_inter:
				if (b == 1):
					hw = hw + 1
			
			# Add noise in HWs (power noise)--------
			hw = hw + np.random.normal(noise_mean, noise_std, 1)
			hw_set.append(hw)
		#print(hw_set)
		#pause()	
		avg_hw_per_set = np.mean(hw_set)
		avg_traces_for_template_matching.append(avg_hw_per_set)
	#print(avg_traces_for_template_matching)	
	#pause()
	#-----MLE estimation for template matching---------
	std_match = noise_std
	max_ll = -9999
	max_ll_idx = -1
	for k, v in key_classes_1.iteritems():
		tem_mean = k[fault_loc][0]
		tem_std = k[fault_loc][1]
		ll = compute_log_likelihood(avg_traces_for_template_matching, mu=tem_mean, sd=tem_std)
		#print(ll)
		if (ll >= max_ll):
			max_ll = ll
			#max_ll_idx = tem_cls
			max_ll_idx = k[fault_loc]
	fault_patt[fault_loc] = max_ll_idx
	#print(fault_patt)
	#pause()
	#print(pat)
	#--------------------------------------------------

# Distance-based template matching----------------

#print(fault_patt)   

min_corr = 4
fault_patt_mean = [fault_patt[0][0], fault_patt[1][0], fault_patt[2][0]]
hw_set_patt_template = None
for k , v in key_classes_1.iteritems():
	sqr_sum = 0
	K_mean = [k[0][0], k[1][0], k[2][0]]
	#print(K_mean)
	#pause()
	for t in range(len(K_mean)):
		sqr_sum = sqr_sum + (K_mean[t] - fault_patt_mean[t])^2		
	Euclidean_dist = np.sqrt(sqr_sum)
	#print(Euclidean_dist)
	if (Euclidean_dist < min_corr):
		min_corr = Euclidean_dist
		hw_set_patt_template = k	

print("Pattern constructed in the template matching...")
print(hw_set_patt_template)   

print("")
print("Recovered value...")
print(key_classes_1[tuple(hw_set_patt_template)])
