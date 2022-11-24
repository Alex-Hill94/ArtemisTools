from __future__ import division
import numpy as np
import eagle_IO.eagle_IO as E
import scipy as sp
import h5py as h5
import illustris_python as il

class LoadData_EAGLE():
    
    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation = "/beegfs2/hpcdata0_backup/simulations/EAGLE/%s/data" % box
        self.alexs_simulation_tensor_path = '/hpcdata7/ariahill/data/%s/%s/tensors5.hdf5' % (box, tag)
        self.alexs_simulation_density_path = '/hpcdata7/ariahill/data/%s/%s/rho_measures_temp.hdf5' % (box, tag)
        self.dmo = dmo

    def query_eagle_table_array(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.data_array = data
    
    def query_eagle_subfind_group(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.data_subfind_group = data

    def query_eagle_table_header(self, qois):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_header("PARTDATA", self.simulation, self.tag, "/Header/%s" % qoi)
            data.append(data_temp)
        self.data_attribute = data

    def query_alexs_tensor_data(self, qois, matter, scale, aperture):
        data = []
        u = h5.File(self.alexs_simulation_tensor_path, 'r')
        for i in range(0, len(qois)):
            qoi = qois[i]            
            data_temp = u['%s/%s/%s/%s' % (matter, scale, aperture, qoi)][:].astype('float64')
            data.append(data_temp)
        u.close()
        self.data_alex_tensor = data

    def query_alexs_density_data(self, qois, sampling, value):
        data = []
        u = h5.File(self.alexs_simulation_density_path, 'r')
        for i in range(0, len(qois)):
            qoi = qois[i]            
            data_temp = u['%s/cut_%s/%s' % (sampling, value, qoi)][:].astype('float64')
            data.append(data_temp)
        u.close()
        self.data_alex_density = data

    def query_rho_data_density_sample(self,qois, aperture_30kpc = True, focus_llim = 8.5, ref_llim = 9.5, centrals_only = False):
        ''' 
        d_3_to_6                 Dataset {4, 22722}
        m_star_focus             Dataset {22722}
        m_star_reference         Dataset {7314}
        n_halo_3_to_6            Dataset {4, 22722}
        read_me                  Dataset {SCALAR}
        rho_3_to_6               Dataset {4, 22722}
        subhalo_IDs_focus        Dataset {22722}
        subhalo_IDs_reference    Dataset {7314}
        '''
        self.centrals_only = centrals_only

        if self.centrals_only:
            cent_path  = '/centrals_only'
            print('Centrals Only')
        elif self.centrals_only == False:
            cent_path = '/all_objects'
            print('All Objects')

        u = h5.File(self.alexs_simulation_density_path, 'r')
        if self.dmo == False:
            print('Loading rho_N data for density sample, selected above %s for the focus sample and above %s for the reference sample. Aperture_30kpc = %s' % (focus_llim, ref_llim, aperture_30kpc))
            data = []
            if aperture_30kpc:
                mass_arg = 'm_star_within_30kpc'
            else:
                mass_arg = 'm_star_within_30kpc'            
            for i in range(0, len(qois)):
                qoi = qois[i]            
                data_temp = u['galaxies/density_selected%s/%s/masses_from_%s_%s/%s' % (cent_path, mass_arg, focus_llim, ref_llim, qoi)][:].astype('float64') 
                data.append(data_temp)
        elif self.dmo == True:
            print('Loading rho_N data for density sample (DMO), selected above %s for the focus sample and above %s for the reference sample.' % (focus_llim, ref_llim))
            data = []
            mass_arg = 'm_subhalo'
            for i in range(0, len(qois)):
                qoi = qois[i]            
                data_temp = u['galaxies/density_selected%s/%s/masses_from_%s_%s/%s' % (cent_path, mass_arg, focus_llim, ref_llim, qoi)][:].astype('float64') 
                data.append(data_temp)
        u.close()
        self.data_alex_density_samp = data

    def query_rho_data_mass_sample(self,qois, aperture_30kpc = True, focus_llim = 8.5):
        '''
        d_3_to_6                 Dataset {4, 22722}
        m_star                   Dataset {22722}
        n_halo_3_to_6            Dataset {4, 22722}
        read_me                  Dataset {SCALAR}
        rho_3_to_6               Dataset {4, 22722}
        subhalo_IDs              Dataset {22722}
        '''
        print('Loading rho_N data for mass sample, selected above %s. Aperture_30kpc = %s' % (focus_llim, aperture_30kpc))
        data = []
        if aperture_30kpc:
            mass_arg = 'm_star_within_30kpc'
        else:
            mass_arg = 'm_star_within_30kpc'
        u = h5.File(self.alexs_simulation_density_path, 'r')
        for i in range(0, len(qois)):
            qoi = qois[i]            
            data_temp = u['galaxies/mass_selected/%s/%s/%s' % (mass_arg, focus_llim, qoi)][:].astype('float64') 
            data.append(data_temp)
        u.close()
        self.data_alex_mass_samp = data

    def query_partdata(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("PARTDATA", self.simulation, self.tag, "/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.partdata = data

    def query_subfind(self, qois, noH = True, physicalUnits = True, verbose = True, call_subnums = False):
        if call_subnums:
            sub_len = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/NumOfSubhalos", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            self.subgrpnum = np.concatenate([np.arange(ii) for ii in sub_len])
            grpnum = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/GroupNumber", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            if min(grpnum) == 1:
                self.grpnum = grpnum
            elif min(grpnum) == 0:
                self.grpnum = grpnum + 1
            else:
                'Weird Group Numbers'
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.subfind_data = data

class LoadData_BAHAMAS():
    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation = "/beegfs2/hpcdata0_backup/simulations/BAHAMAS/%s" % box

    def query_eagle_table_header(self, qois):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_header("PARTDATA", self.simulation, self.tag, "/Header/%s" % qoi)
            data.append(data_temp)
        self.data_attribute = data

    def query_eagle_table_array(self, qois, noH = True, physicalUnits = True, verbose = True, call_subnums = False):
        if call_subnums:
            sub_len = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/NumOfSubhalos", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            self.subgrpnum = np.concatenate([np.arange(ii) for ii in sub_len])
            grpnum = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/GroupNumber", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            if min(grpnum) == 1:
                self.grpnum = grpnum
            elif min(grpnum) == 0:
                self.grpnum = grpnum + 1
            else:
                'Weird Group Numbers'
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.data_array = data
    
    def query_eagle_subfind_group(self, qois, noH = True, physicalUnits = True, verbose = True):
        
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.data_subfind_group = data

class LoadData_BAHAMAS_XL():
    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation = "/beegfs2/hpcdata0_backup/simulations/BAHAMAS_XL/%s/data" % box

    def query_eagle_table_array(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.data_array = data
    
    def query_eagle_subfind_group(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.data_subfind_group = data

    def query_eagle_table_header(self, qois):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_header("PARTDATA", self.simulation, self.tag, "/Header/%s" % qoi)
            data.append(data_temp)
        self.data_attribute = data

class LoadData_TNG():
    '''
    Loads IllustrisTNG data. Note that the units are conventionally in ckpc/h. 
    There is no flag to change this, so the user will need to do this themselves
    '''
    def __init__(self, box, tag):
        self.box = box
        self.tag = tag
        self.simulation = '/beegfs2/hpcdata0_backup/simulations/IllustrisTNG/'+self.box+'/output/'
        header = il.groupcat.loadHeader(self.simulation,99)
        self.h = header['HubbleParam']
        self.z = header['Redshift']
        self.a = 1./(1. + self.z)
        self.boxsize = header['BoxSize']*self.a/(1000. * self.h)

    def query_subhalo_data(self, qois, calc_subnum = True):
        if calc_subnum:
            sub_len = il.groupcat.loadHalos(self.simulation,self.tag,fields=['GroupNsubs'])
            self.subgrpnum = np.concatenate([np.arange(ii) for ii in sub_len])
        fields_sub = qois + ['SubhaloGrNr']
        subhaloes = il.groupcat.loadSubhalos(self.simulation,self.tag,fields=fields_sub)
        self.grpnum = subhaloes['SubhaloGrNr'] + 1
        self.data_sub = subhaloes

    def query_group_data(self, qois):
        fields_sub = qois
        haloes = il.groupcat.loadHalos(self.simulation,self.tag,fields=fields_sub)
        self.data_grp = haloes

    def query_jon_data(self, qois):
        path = "/beegfs2/hpcdata0_backup/arijdav/TNG_catalogues/%s/snap_0%s/catalogue.hdf5" % (self.box, self.tag)
        u = h5.File(path, 'r')
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]            
            data_temp = u['%s' % qoi][:].astype('float64')
            data.append(data_temp)
        self.jon_data = data

class LoadDataBAHAMAS_Special():
    # What is this for?
    # Includes partdata and snapshot
    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation_subfind = "/beegfs2/hpcdata0_backup/simulations/BAHAMAS/%s/Data/EagleSubGroups_5r200" % box
        self.simulation_snapshot = "/beegfs2/hpcdata0_backup/simulations/BAHAMAS/%s/Data/Snapshots" % box

    def query_header(self, qois):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_header("PARTDATA", self.simulation_subfind, self.tag, "/Header/%s" % qoi)
            data.append(data_temp)
        self.header_data = data

    def query_subfind(self, qois, noH = True, physicalUnits = True, verbose = True, call_subnums = False):
        if call_subnums:
            sub_len = E.read_array("SUBFIND_GROUP", self.simulation_subfind, self.tag, "/FOF/NumOfSubhalos", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            self.subgrpnum = np.concatenate([np.arange(ii) for ii in sub_len])
            grpnum = E.read_array("SUBFIND", self.simulation_subfind, self.tag, "/Subhalo/GroupNumber", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            if min(grpnum) == 1:
                self.grpnum = grpnum
            elif min(grpnum) == 0:
                self.grpnum = grpnum + 1
            else:
                'Weird Group Numbers'
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation_subfind, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.subfind_data = data
    
    def query_fof(self, qois, noH = True, physicalUnits = True, verbose = True):
        
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation_subfind, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SUBFIND_GROUP", self.simulation_subfind, self.tag, "/FOF/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.fof_data = data

    def query_partdata(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation_subfind, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("PARTDATA", self.simulation_subfind, self.tag, "/PartType1/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.partdata = data

    def query_snapshot(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation_snapshot, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SNAPSHOT", self.simulation_snapshot, self.tag, "/PartType1/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.snapshot_data = data

class LoadData_ARTEMIS():

    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation = '/hpcdata4/simulations/ARTEMIS/Gal_Names/%s/data' % box

    def query_header(self, qois):
        data = []
        file = self.simulation + '/snapshot_029_z000p000/snap_029_z000p000.0.hdf5'
        print(file)
        hf = h5.File(file, 'r')
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = hf['Header'].attrs[qoi]
            #data_temp = E.read_header("PARTDATA", self.simulation, self.tag, "/Header/%s" % qoi)
            data.append(data_temp)
        self.header_data = data

    def query_subfind(self, qois, noH = True, physicalUnits = True, verbose = True, call_subnums = False):
        if call_subnums:
            sub_len = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/NumOfSubhalos", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            self.subgrpnum = np.concatenate([np.arange(ii) for ii in sub_len])
            grpnum = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/GroupNumber", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            if min(grpnum) == 1:
                self.grpnum = grpnum
            elif min(grpnum) == 0:
                self.grpnum = grpnum + 1
            else:
                'Weird Group Numbers'
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.subfind_data = data
    
    def query_fof(self, qois, noH = True, physicalUnits = True, verbose = True):
        
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.fof_data = data

    def query_partdata(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("PARTDATA", self.simulation, self.tag, "/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.partdata = data

    def query_snapshot(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SNAPSHOT", self.simulation, self.tag, "/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.snapshot_data = data

class LoadData_C_EAGLE():

    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation = '/beegfs2/hpcdata0_backup/simulations/C-EAGLE/%s/data' % box

    def query_header(self, qois):
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_header("PARTDATA", self.simulation, self.tag, "/Header/%s" % qoi)
            data.append(data_temp)
        self.header_data = data

    def query_subfind(self, qois, noH = True, physicalUnits = True, verbose = True, call_subnums = False):
        if call_subnums:
            sub_len = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/NumOfSubhalos", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            self.subgrpnum = np.concatenate([np.arange(ii) for ii in sub_len])
            grpnum = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/GroupNumber", noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            if min(grpnum) == 1:
                self.grpnum = grpnum
            elif min(grpnum) == 0:
                self.grpnum = grpnum + 1
            else:
                'Weird Group Numbers'
        data = []
        for i in range(0, len(qois)):
            qoi = qois[i]
            data_temp = E.read_array("SUBFIND", self.simulation, self.tag, "/Subhalo/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.subfind_data = data
    
    def query_fof(self, qois, noH = True, physicalUnits = True, verbose = True):
        
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SUBFIND_GROUP", self.simulation, self.tag, "/FOF/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.fof_data = data

    def query_partdata(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("PARTDATA", self.simulation, self.tag, "/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.partdata = data

    def query_snapshot(self, qois, noH = True, physicalUnits = True, verbose = True):
        data = []
        qois = np.array((qois))
        for i in range(0, len(qois)):
            qoi = qois[i]
            print(self.simulation, self.tag, "/FOF/%s" % qoi)
            data_temp = E.read_array("SNAPSHOT", self.simulation, self.tag, "/%s" % qoi, noH = noH, physicalUnits = physicalUnits, verbose = verbose)
            data.append(data_temp)
        self.snapshot_data = data
