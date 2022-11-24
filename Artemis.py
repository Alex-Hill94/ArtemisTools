from __future__ import division
import numpy as np
from scipy.linalg import norm as norm
from scipy.stats import binned_statistic as binned
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from multiprocessing import Pool
import matplotlib.gridspec as gridspec
from LoadData import LoadData_ARTEMIS as LoadArtemis
from tens_3d import reduced_iterative
from astropy import units as u
from astropy import constants as c
from scipy.optimize import curve_fit
import pandas as pd

'''
Reading script, note EAGLE and BAHAMAS can be read in h-less and physical, Illustris can't
EAGLE and BAHAMAS are Mpc
Illustris kpc
'''

class Central_Galaxy_and_Halo():

	def __init__(self, box, tag):
		self.box = box
		self.tag = tag

	def load_data(self):

		'''
		Uses the LoadArtemis packaage to read in data
		'''
		e = LoadArtemis(box = self.box, tag = self.tag)

		# Grab header information
		e.query_header(qois = ('HubbleParam', 'ExpansionFactor', 'BoxSize'))
		h, a_0, box_size = e.header_data[0], e.header_data[1], e.header_data[2]

		# Boxsize in h-less physical units
		BS = box_size*a_0/h

		# Mass of a DM particle
		dm_partmass = 1.2*1e5#*h**(-1)

		# Grab the star particle data: coordinates, metallicity, group number, subgroup number, mass, formation time
		e.query_partdata(qois = ('PartType4/Coordinates', 'PartType4/Metallicity','PartType4/GroupNumber', 'PartType4/SubGroupNumber', 'PartType4/Mass', 'PartType4/StellarFormationTime'))
		Star_Coordinates, Star_Metallicity, Star_GroupNumber, Star_SubGroupNumber, Star_Mass, Star_FormationTime = e.partdata[0], e.partdata[1], e.partdata[2], e.partdata[3], e.partdata[4], e.partdata[5]

		# Grab the DM particle data: coordinates, velocity, group number, subgroup number
		e.query_partdata(qois = ('PartType1/Coordinates', 'PartType1/Velocity', 'PartType1/GroupNumber', 'PartType1/SubGroupNumber'))
		DM_Coordinates, DM_Velocity, DM_GroupNumber, DM_SubGroupNumber = e.partdata[0], e.partdata[1], e.partdata[2], e.partdata[3]
		DM_mass = np.ones(len(DM_SubGroupNumber)) * dm_partmass

		# Grab the star particle data: coordinates, sfr, group number, subgroup number
		e.query_partdata(qois = ('PartType0/Coordinates', 'PartType0/StarFormationRate', 'PartType0/GroupNumber', 'PartType0/SubGroupNumber'))
		Gas_Coordinates, Gas_SFR, Gas_GroupNumber, Gas_SubGroupNumber = e.partdata[0], e.partdata[1], e.partdata[2], e.partdata[3]

		# Loads the subfind data, which relates to the subhalo: coordinates, mass, group number, subgroup number, effective radius
		e.query_subfind(qois = ('CentreOfPotential', 'Mass', 'GroupNumber', 'SubGroupNumber', 'HalfMassRad'))
		Sub_COP, Sub_Mass, Sub_GN, Sub_SN, Sub_Reff, Sub_DMRad = e.subfind_data[0], e.subfind_data[1], e.subfind_data[2], e.subfind_data[3], e.subfind_data[4][:,4], e.subfind_data[4][:,1]

		# Loads friends of friends data: r200, m200
		e.query_fof(qois = ('Group_R_Crit200', 'Group_M_Crit200'))
		FOF_r200, FOF_m200 = e.fof_data[0], e.fof_data[1]


		Sub_R200, Sub_M200 = FOF_r200[Sub_GN - 1], FOF_m200[Sub_GN - 1]

		# I'm only interested in the central subhalo of the zoom simulation, so I select the biggest FOF halo and the associated central
		central_sub = (Sub_GN == 1) * (Sub_SN == 0)

		# The centre of potential of the galaxy, as well as its radius and mass
		subhalo_centre = Sub_COP[central_sub]
		subhalo_radius = Sub_R200[central_sub]
		subhalo_mass   = Sub_M200[central_sub]*1e10

		# Reduces the partdata to just those associated with the galaxy of interest
		scs 		= Star_Coordinates[(abs(Star_GroupNumber) == 1)* (Star_SubGroupNumber == 0)]
		star_mass 	= Star_Mass[(abs(Star_GroupNumber) == 1)* (Star_SubGroupNumber == 0)]
		star_ft 	= Star_FormationTime[(abs(Star_GroupNumber) == 1)* (Star_SubGroupNumber == 0)]
		dms 		= DM_Coordinates[(abs(DM_GroupNumber) == 1)* (DM_SubGroupNumber == 0)]
		dm_mass 	= DM_mass[(abs(DM_GroupNumber) == 1)* (DM_SubGroupNumber == 0)]
		sfgs		= Gas_Coordinates[(abs(Gas_GroupNumber) == 1)* (Gas_SubGroupNumber == 0) * (Gas_SFR > 0)]
		sfg_sfr		= Gas_SFR[(abs(Gas_GroupNumber) == 1)* (Gas_SubGroupNumber == 0) * (Gas_SFR > 0)]

		# Centre all coordiantes with respect to the centre of potential of the main galaxy
		# Inclusion of boxsize stuff relates to the possibility of wrapping: because
		# we're dealing with a periodic volume, you can have a bound galaxy with particles
		# on either side of the volume. This calculation ensures that this is accounted for
		SCS = (scs - subhalo_centre+BS/2.)%BS - BS/2.
		DMS = (dms - subhalo_centre+BS/2.)%BS - BS/2.
		SFG_S = (sfgs - subhalo_centre+BS/2.)%BS - BS/2.

		# Data I'm interested in using going forward
		self.h		   	= h
		self.star_locs 	= SCS
		self.star_mass 	= star_mass
		self.star_ft 	= star_ft
		self.sfg_locs   = SFG_S
		self.sfg_mass   = sfg_sfr		
		self.dm_locs   	= DMS
		self.dm_mass   	= dm_mass

		self.halo_r200 	= subhalo_radius
		self.halo_m200 	= subhalo_mass

	def get_dm_profile(self, n_bins = 1000, show_fig = True, save_fig = False):
		'''
		Computes the NFW profile of the halo
		'''
		radii = norm(self.dm_locs, axis = 1)
		anulus_bins = np.linspace(0., self.halo_r200, n_bins)
		anulus_vol = np.zeros(len(anulus_bins) - 1)

		anulus_bins_r200 = anulus_bins/self.halo_r200
		r_r200_plot = (anulus_bins_r200[1:] + anulus_bins_r200[:-1]) / 2

		def vol_sphere(r):
			return 4./3. * np.pi * r **3

		def vol_diff(r1, r2):
			v1 = vol_sphere(r1)
			v2 = vol_sphere(r2)
			return abs(v2 - v1)

		def NFW(r_r200, delta, c):
			rho_rhocrit = (delta)/((r_r200 * c) * (1 + r_r200 * c)**2)
			return rho_rhocrit


		for j in range(0, len(anulus_bins)-1):
			anulus_vol[j] = vol_diff(anulus_bins[j], anulus_bins[j+1])

		vals, be1, be2 = binned(radii, values = self.dm_mass, bins = anulus_bins.T, statistic = np.sum)

		densities = vals/anulus_vol
	
		G_prime = c.G.to(u.Mpc**3/(u.Msun*u.s**2))

		hubble = self.h * (u.km/u.s)/u.Mpc
		h_prime = hubble.to((u.Mpc/u.s)/u.Mpc)

		# Probably wrong units
		rho_crit = 3.*(100.*h_prime)**2./(8. * np.pi * G_prime)
		self.r_r200_plot = r_r200_plot
		self.ys = densities/rho_crit.value


		def NFW(r, c, rho_0):
			# r is in terms of Mpc
			R_s = self.halo_r200[0]/c
			num = rho_0
			denom = (r/R_s) * (1. + r/R_s)**2.
			rho_r = num/denom
			return rho_r

		x = self.r_r200_plot * self.halo_r200

		blanks = np.zeros(len(x))
		blanksy = np.zeros(len(x))

		blanks[:] = x.T[0][:]

		blanksy[:] = (densities)[:]

		popt, pcov = curve_fit(NFW, blanks, blanksy, p0 = [15., 1e15])

		fit_vals = NFW(x.T[0], c = popt[0], rho_0 = popt[1])

		rho_0 = np.format_float_scientific(popt[1], precision = 1)
		rho_0 = np.format_float_scientific(popt[1], precision = 1)

		if show_fig:
			FS = 15
			matplotlib.rcParams.update({'font.size': FS})
			fig, ax1 = plt.subplots(1,1, figsize = [6., 6.])
			ax1.plot(r_r200_plot, densities/rho_crit.value, label = '$\mathrm{Data}$')
			ax1.plot(r_r200_plot, fit_vals//rho_crit.value, label = '$\mathrm{NFW:} c = %s, $$\\rho_{0} = %s$' % (round(popt[0], 3), rho_0))
			ax1.set_yscale('log')
			ax1.set_xscale('log')
			ax1.set_xlabel('$r/R_{\mathrm{200,\ crit}}$')
			ax1.set_ylabel('$\\rho(r)/\\rho_{\mathrm{crit}}$')
			ax1.tick_params(which = 'both', direction= 'in', right =True, top = True)
			plt.legend()
			plt.subplots_adjust(wspace=0.0)
			plt.show()
		if save_fig:
			FS = 15
			matplotlib.rcParams.update({'font.size': FS})
			fig, ax1 = plt.subplots(1,1, figsize = [6., 6.])
			ax1.plot(r_r200_plot, densities/rho_crit.value, label = '$\mathrm{Data}$')
			ax1.plot(r_r200_plot, fit_vals//rho_crit.value, label = '$\mathrm{NFW:} c = %s, $$\\rho_{0} = %s$' % (round(popt[0], 3), rho_0))
			ax1.set_yscale('log')
			ax1.set_xscale('log')
			ax1.set_xlabel('$r/R_{\mathrm{200,\ crit}}$')
			ax1.set_ylabel('$\\rho(r)/\\rho_{\mathrm{crit}}$')
			ax1.tick_params(which = 'both', direction= 'in', right =True, top = True)
			plt.legend()
			plt.subplots_adjust(wspace=0.0)
			plt.savefig('%s_%s_dm_profile.png' % (self.box, self.tag) )
		
		self.density 	  = densities
		self.r_r200_plot  = r_r200_plot
		self.anulus_bins  = anulus_bins
		self.rho_crit 	  = rho_crit.value
		self.nfw_fit	  = popt

	def rotate_frame(self, aperture = 0.03, plot = False, n_parts = None):
		'''
		Rotates the particles to be in the frame of reference defined by the iterative reduced inertia
		tensor. A default aperture of 30pkpc is applied.
		'''

		M, sl = reduced_iterative(self.star_mass, self.star_locs, aperture)
		eigval, eigvec = np.linalg.eig(M)
		order = np.argsort(eigval)
		
		# The transformation matrix is made of the eigenvalues of the tensor
		trans_mat = eigvec[:,order].T
		self.star_locs_rot = np.squeeze(np.matmul(trans_mat, self.star_locs[:,:, np.newaxis]))
		self.dm_locs_rot = np.squeeze(np.matmul(trans_mat, self.dm_locs[:,:, np.newaxis]))
		self.sfg_locs_rot = np.squeeze(np.matmul(trans_mat, self.sfg_locs[:,:, np.newaxis]))

		if plot:
			if n_parts is None:
				n_parts = 10000
			
			star_rands = np.random.rand(len(self.star_locs_rot))
			star_frac		= n_parts/len(self.star_locs_rot)
			star_plot_locs 	= self.star_locs_rot[star_rands < star_frac]
			fig, axs = plt.subplots(1, 3, figsize = [12, 4])

			axs[0].scatter(star_plot_locs[:,1], star_plot_locs[:,0], s = 0.01)
			axs[1].scatter(star_plot_locs[:,2], star_plot_locs[:,0], s = 0.01)
			axs[2].scatter(star_plot_locs[:,2], star_plot_locs[:,1], s = 0.01)

			axs[0].set_title('Edge On (Intermediate/Minor)')
			axs[1].set_title('Edge On (Major/Minor)')
			axs[2].set_title('Face On (Major/Intermediate)')
			axs[0].set_xlim([-0.03, 0.03])
			axs[1].set_xlim([-0.03, 0.03])			
			axs[2].set_xlim([-0.03, 0.03])
			axs[0].set_ylim([-0.03, 0.03])
			axs[1].set_ylim([-0.03, 0.03])			
			axs[2].set_ylim([-0.03, 0.03])
			plt.savefig('test.png')
			
			plt.close()

	def get_shapes(self, star_aperture = 0.03, dm_inner_ap = False):
		'''
		Returns the dimensions of the best-fitting iterative reduced inertia tensor. 
		'''
		
		if dm_inner_ap:
			dm_aperture = 0.03
		else:
			dm_aperture = self.halo_r200

		tensor_comp = reduced_iterative

		def get_ellipsoid(locs, weights, aperture, scale = True):
			## scale = True scales the axis lengths such that the major axis 
			# length is the same as the initial aperture
			M, _ = tensor_comp(locs, weights, aperture)
			eigval, eigvec = np.linalg.eig(M)
			order = np.argsort(eigval)
			trans_mat = eigvec[:,order].T
			a,b,c = eigval[order]
			if scale == True:
				ratio = c/aperture
				a, b, c = a/ratio, b/ratio, c/ratio
			return (a, b, c), trans_mat

		self.dm_aperture = dm_aperture
		DM_axis_lens, DM_vectors = get_ellipsoid(self.dm_mass, self.dm_locs_rot, dm_aperture[0])
		Star_axis_lens, Star_vectors = get_ellipsoid(self.star_mass, self.star_locs_rot, star_aperture)
		#SFG_axis_lens, SFG_vectors = get_ellipsoid(self.sfg_mass, self.sfg_locs_rot, star_aperture)


		self.star_aperture = star_aperture
		self.DM_ellipsoid = (DM_axis_lens, DM_vectors)
		self.Star_ellipsoid = (Star_axis_lens, Star_vectors)
		#self.SSF_ellipsoid = (SFG_axis_lens, SFG_vectors)

	def data_save(self, reduce_particles = 10000, grab_young = True, frac_young_stars = 0.75):

		exp_f = 0.9318
		
		dm_rands = np.random.rand(len(self.dm_locs_rot))
		dm_frac		= reduce_particles/len(self.dm_locs_rot)

		sfg_rands 		= np.random.rand(len(self.sfg_locs_rot))
		sfg_frac		= reduce_particles/len(self.sfg_locs_rot)


		if grab_young:
			young_locs = self.star_locs_rot[self.star_ft > exp_f]
			old_locs 	= self.star_locs_rot[self.star_ft < exp_f]
			reduce_young = reduce_particles * frac_young_stars
			reduce_old = reduce_particles * (1-frac_young_stars)

			star_young_rands = np.random.rand(len(young_locs))
			star_young_frac   = reduce_young/len(young_locs)
			young_locs_red 				= young_locs[star_young_rands < star_young_frac]

			star_old_rands = np.random.rand(len(old_locs))
			star_old_frac   = reduce_old/len(old_locs)
			old_locs_red 	= old_locs[star_old_rands < star_old_frac]

			#self.STAR_LOCATIONS 	= self.star_locs_rot[star_rands < star_frac]
			self.STAR_AGED_LOCS 	= np.vstack((young_locs_red, old_locs_red))
		else:
			star_rands = np.random.rand(len(self.star_locs_rot))
			star_frac		= reduce_particles/len(self.star_locs_rot)
			self.STAR_AGED_LOCS 		= self.star_locs_rot[star_rands < star_frac]


		self.SFG_LOCATIONS 		= self.sfg_locs_rot[sfg_rands < sfg_frac]


		self.DM_LOCATIONS 		= self.dm_locs_rot[dm_rands < dm_frac]
		self.DM_axis_lengths 	= self.DM_ellipsoid[0]
		self.DM_vector_dirs 	= self.DM_ellipsoid[1]
		self.STAR_axis_lengths 	= self.Star_ellipsoid[0]
		self.STAR_vector_dirs 	= self.Star_ellipsoid[1]
		self.DM_density_params  = self.nfw_fit
		self.DM_radius			= self.halo_r200
		self.STAR_radius		= self.star_aperture

		dm_locs = pd.DataFrame(data = self.DM_LOCATIONS, columns = ('DM x', 'DM y', 'DM z'))
		star_locs = pd.DataFrame(data = self.STAR_AGED_LOCS, columns = ('Star x', 'Star y', 'Star z'))
		sfg_locs = pd.DataFrame(data = self.SFG_LOCATIONS, columns = ('SFG x', 'SFG y', 'SFG z'))

		#dm_locs.to_csv('dm.txt', index=False)  
		#star_locs.to_csv('star.txt', index=False)  

		self.dm_halo_information = 'DM Halo Axis Lengths = %s. DM Halo Vectors = %s. DM Radius = %s Mpc. c = %s, rho_0 = %s. ################################################' % (self.DM_axis_lengths, self.DM_vector_dirs, self.DM_radius, self.DM_density_params[0], self.DM_density_params[1])
		self.galaxy_information = 'Star Axis Lengths = %s. Star Vectors = %s. Star Radius = %s Mpc. ################################################' % (self.STAR_axis_lengths, self.STAR_vector_dirs, self.STAR_radius)

		dm_filename = 'dm_halo_info_%s_%s.txt' % (self.box, self.tag)
		star_filename = 'star_info_%s_%s.txt' % (self.box, self.tag)
		sfg_filename = 'sfg_info_%s_%s.txt' % (self.box, self.tag)

		np.savetxt(fname = dm_filename, X = dm_locs, header = self.dm_halo_information)
		np.savetxt(fname = star_filename, X = star_locs, header = self.galaxy_information)
		np.savetxt(fname = sfg_filename, X = sfg_locs)#, header = self.galaxy_information)
		
		plt.figure()
		plt.scatter(self.DM_LOCATIONS[:,0], self.DM_LOCATIONS[:,1], s = 0.01, label = 'stars')
		plt.scatter(self.STAR_AGED_LOCS[:,0], self.STAR_AGED_LOCS[:,1], s = 0.01, label = 'dm')
		plt.legend()
		plt.savefig('%s_%s_scatterplot.png' % (self.box, self.tag))
		plt.close()




if __name__ == '__main__':
	tags = ("009_z005p037", "016_z002p012", "020_z001p004", "024_z000p503", "029_z000p000")
	tag = "029_z000p000"
	G = Central_Galaxy_and_Halo("G42", tag)
	G.load_data()
	#G.get_dm_profile(n_bins = 50, show_fig = False, save_fig = True)
	G.rotate_frame(plot = True)
	#G.get_shapes()
	#G.data_save(reduce_particles=10000, grab_young=False)


