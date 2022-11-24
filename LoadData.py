from __future__ import division
import numpy as np
import eagle_IO.eagle_IO as E
import scipy as sp
import h5py as h5
#import illustris_python as il

class LoadData_EAGLE():
    
    def __init__(self, box, tag, dmo = False):
        self.box = box
        self.tag = tag
        self.simulation = "/beegfs2/hpcdata0_backup/simulations/EAGLE/%s/data" % box
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
