from __future__ import division
import numpy as np
import scipy as sp
from numpy import linalg
#from numpy import *



def ellip(A):
	e, v = linalg.eig(A)
	e = sort(sqrt(e)) # Smallest to largest
	a, b, c = e[2], e[1], e[0]
	S = c/a
	T = (a**2 - b**2)/(a**2 - c**2)
	return S, T

def s_and_q(MATRIX):
	vals, vecs 	= sp.linalg.eigh(MATRIX)
	abc			= np.sqrt(vals)  ## Some negative values here, why?
	order		= np.argsort(abc)
	a, b, c		= abc[order][0], abc[order][1], abc[order][2]
	s 			= a/c
	q			= b/c
	return s, q	

def conv(r, r1):
	return abs(1. - (r1/r))

#############################################################################################################################################################################################################################################################


def reduced(weight, locs, r_sph):
	M = np.ones((3,3))
	M[:] = np.nan
	R_matrix 	= np.array(locs[:,:,np.newaxis] * locs[:,np.newaxis,:])
	locs_sqr 	= np.square(locs)
	r_tilde_sqr = np.sum(np.divide(locs_sqr, np.square(r_sph)), axis = 1) # Working
	WEIGHT		= np.divide(weight, r_tilde_sqr)
	corr  		= r_tilde_sqr <= 1.
	if sum(corr) > 2:
		mass_weighted_mat = np.squeeze(np.multiply(R_matrix, WEIGHT[:,np.newaxis, np.newaxis]))
		M 			=  np.sum(mass_weighted_mat[corr], axis = 0)
		sum_weight	=  np.sum(WEIGHT[corr], axis = 0) ### ERROR HERE?
		M 			=  np.divide(M, sum_weight) 
	return M

def reduced_iterative(weight, locs, r_sph):
	R_matrix 		= np.array(locs[:,:,np.newaxis] * locs[:,np.newaxis,:])
	M			 	= reduced(weight, locs, r_sph)  # TIME INTENSIVE
	final_M 		= np.empty((3,3,))
	final_M[:]		= np.nan
	max_it  		= np.nan
	sub_length		= np.nan
	if (np.isnan(M).any() == False) * (np.isinf(M).any() == False) == True:
		s, q		= s_and_q(M)
		for iteration in range(1, 100):
			eigenvalues, eigenvectors = sp.linalg.eigh(M)
			order			= np.argsort(eigenvalues).astype('int64')
			a, b, c			= np.sqrt(eigenvalues[order])   # Negative value in sqrt returns error 
			scale_term		= r_sph * (a*b*c)**(-1./3.)
			A, B, C			= a*scale_term, b*scale_term, c*scale_term
			P_ab			= eigenvectors[:,order].T
			new_locs 		= np.squeeze(np.matmul(P_ab, locs[:,:, np.newaxis]))
			nlocs_sqr		= np.square(new_locs)
			a_sqrs			= np.repeat(A**2, len(new_locs))
			b_sqrs			= np.repeat(B**2, len(new_locs))
			c_sqrs			= np.repeat(C**2, len(new_locs))
			if len(new_locs) < 2:
				break
			else:
				divisors		= np.zeros(np.shape(nlocs_sqr))
				divisors[:,0], divisors[:,1], divisors[:,2] = a_sqrs, b_sqrs, c_sqrs
				r_tilde_sqr 	= np.sum(nlocs_sqr/divisors, axis = 1)	
				WEIGHT			= weight/r_tilde_sqr
				corr 			= r_tilde_sqr <= 1.
				length			= len(WEIGHT[(corr)*(WEIGHT > 0)])
			if length < 10:
				break
			else:
				mass_weighted_mat = np.squeeze(np.multiply(R_matrix, WEIGHT[:,np.newaxis, np.newaxis]))
				M1 				=  np.sum(mass_weighted_mat[corr], axis = 0)
				sum_weight		=  np.sum(WEIGHT[corr], axis = 0) ### ERROR HERE?
				M1 				=  M1/sum_weight
				if np.isnan(M1).any() == True:
					break
				if np.isinf(M1).any() == True:
					break
				s1, q1 =  s_and_q(M1)		
				convergence		= (conv(s, s1) < 0.01) * (conv(q, q1) < 0.01)
				if convergence == True:
					final_M 	= M1
					max_it 		= iteration
					sub_length	= length
					break					
				else:
					M, s, q = M1, s1, q1
	return final_M, sub_length

def simple(weight_pr, locs_pr, r_sph):
	M = np.ones((3,3))
	M[:] = np.nan
	dists = np.linalg.norm(locs_pr, axis = 1)
	aperture = dists < r_sph
	weight, locs = weight_pr[aperture], locs_pr[aperture]
	if len(weight)>20:
		R_matrix 	= np.array(locs[:,:,np.newaxis] * locs[:,np.newaxis,:])
		mass_weighted_mat = np.squeeze(np.multiply(R_matrix, weight[:,np.newaxis, np.newaxis]))
		M 			=  np.sum(mass_weighted_mat, axis = 0)
		sum_weight	=  np.sum(weight, axis = 0)
		M 			=  np.divide(M, sum_weight) 
	return M

def simple_iterative(weight, locs, r_sph):
	R_matrix 		= np.array(locs[:,:,np.newaxis] * locs[:,np.newaxis,:])
	M			 	= simple(weight, locs, r_sph)  # TIME INTENSIVE
	final_M 		= np.empty((3,3,))
	final_M[:]		= np.nan
	max_it  		= np.nan
	sub_length		= np.nan
	if (np.isnan(M).any() == False) * (np.isinf(M).any() == False) == True:
		s, q		= s_and_q(M)
		for iteration in range(1, 100):
			eigenvalues, eigenvectors = sp.linalg.eigh(M)
			order			= np.argsort(eigenvalues).astype('int64')
			a, b, c			= np.sqrt(eigenvalues[order])   # Negative value in sqrt returns error 
			scale_term		= r_sph * (a*b*c)**(-1./3.)
			A, B, C			= a*scale_term, b*scale_term, c*scale_term
			P_ab			= eigenvectors[:,order].T
			new_locs 		= np.squeeze(np.matmul(P_ab, locs[:,:, np.newaxis]))
			nlocs_sqr		= np.square(new_locs)
			a_sqrs			= np.repeat(A**2, len(new_locs))
			b_sqrs			= np.repeat(B**2, len(new_locs))
			c_sqrs			= np.repeat(C**2, len(new_locs))
			if len(new_locs) < 2:
				break
			else:
				divisors		= np.zeros(np.shape(nlocs_sqr))
				divisors[:,0], divisors[:,1], divisors[:,2] = a_sqrs, b_sqrs, c_sqrs
				r_tilde_sqr 	= np.sum(nlocs_sqr/divisors, axis = 1)	
				WEIGHT			= weight/r_tilde_sqr
				corr 			= r_tilde_sqr <= 1.
				length			= len(WEIGHT[(corr)*(WEIGHT > 0)])
			if length < 10:
				break
			else:
				mass_weighted_mat = np.squeeze(np.multiply(R_matrix, weight[:,np.newaxis, np.newaxis]))
				M1 				=  np.sum(mass_weighted_mat[corr], axis = 0)
				sum_weight		=  np.sum(weight[corr], axis = 0) ### ERROR HERE?
				M1 				=  M1/sum_weight
				if np.isnan(M1).any() == True:
					break
				if np.isinf(M1).any() == True:
					break
				s1, q1 =  s_and_q(M1)		
				convergence		= (conv(s, s1) < 0.01) * (conv(q, q1) < 0.01)
				if convergence == True:
					final_M 	= M1
					max_it 		= iteration
					sub_length	= length
					break					
				else:
					M, s, q = M1, s1, q1
	return final_M, sub_length
