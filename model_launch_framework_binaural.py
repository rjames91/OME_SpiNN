import spinnaker_graph_front_end as g

from OME_vertex import OMEVertex
from DRNL_vertex import DRNLVertex
from IHCAN_vertex import IHCANVertex
from MCack_vertex import MCackVertex
import model_binaries

from pacman.model.constraints.placer_constraints\
    .chip_and_core_constraint import ChipAndCoreConstraint
from pacman.model.graphs.machine import MachineEdge

from spinnman.model.enums.cpu_state import CPUState

from spinn_front_end_common.utilities import globals_variables

import numpy
import logging
import time
import math

logger = logging.getLogger(__name__)

def calculate_additional_chips(num_rows,num_macks):
    acc=0.
    for i in range(num_rows-1):
        acc += num_macks**i
    return int(numpy.ceil(acc/16.))

def run_model(
        data, n_chips=None,n_drnl=0,pole_freqs=[4000],n_ihcan=0,
        fs=44100,resample_factor=1,num_macks=4,seg_size=96,bitfield=True,
        rt=True,profile=True):
    if data.shape[0]>1:
        n_ears = 2
    else:
        n_ears = 1

    # calculate number of mack tree rows above OME (final row is all DRNL instances)
    n_mack_tree_rows = int(numpy.ceil(math.log(len(pole_freqs),num_macks)))
    #calculate remainder drnls, if non-zero once this many are generated all following mack-->drnl connections will be 1 to 1
    remainder_macks =(num_macks**(n_mack_tree_rows)) - len(pole_freqs)
    # remainder_macks =(math.factorial(n_mack_tree_rows)*num_macks) - len(pole_freqs)

    additional_chips = calculate_additional_chips(num_macks=num_macks, num_rows=n_mack_tree_rows)
    additional_chips *= n_ears

    # Set up the simulation
    if n_ears>1:
        #ensure at least 2 boards are allocated (1 ethernet chip per ear)
        requested_n_chips = max(n_chips+additional_chips,2*48)
    else:
        requested_n_chips = n_chips+additional_chips
    print "requested_n_boards:{}".format(numpy.ceil(requested_n_chips/48.))
    g.setup(n_chips_required=requested_n_chips, model_binary_module=model_binaries)
    # Get the number of cores available for use
    machine = g.machine()
    boards = dict()

#repeat single ear partitioning etc for both ears - the only thing that should differ is the data input to OME
    binaural_drnls = []
    left_ear_chips = []
    binaural_omes = []
    for ear in range(n_ears):
        #changed to lists to ensure data is read back in the same order that vertices are instantiated
        # omes=list()
        drnls=list()
        ihcans=list()
        macks=list()

        delays =  numpy.round(len(pole_freqs)*numpy.random.rand(len(pole_freqs)))
        if n_ihcan>0:
            #create unique seeds for IHCAN instances
            n_ihcans = n_chips*n_drnl*n_ihcan
            random_range = numpy.arange(n_ihcans*4,dtype=numpy.uint32)
            seeds = numpy.random.choice(random_range,int(n_ihcans*4),replace=False)
        seed_index = 0
        cf_index = 0
        count = 0

        #OME is on ethernet chip for live streaming
        for chip in machine.ethernet_connected_chips:
            if (chip.x,chip.y) not in left_ear_chips:
                # create OME
                ome = OMEVertex(data[ear], fs, len(pole_freqs),rt=rt,profile=profile)
                g.add_machine_vertex_instance(ome)
                # constrain placement to local chip
                ome.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
                binaural_omes.append(ome)
                boards[chip.x, chip.y] = chip.ip_address
                ome_chip_xy = [(chip.x, chip.y)]
                if ear == 0:
                    left_ear_chips.append((chip.x, chip.y))
                break

        for chip in machine.chips:
            if count >= n_chips/n_ears:
                break
            #else:
            elif (chip.x, chip.y) not in ome_chip_xy and (chip.x, chip.y) not in left_ear_chips:
                #create DRNLs
                for i in range(n_drnl):

                    CF=pole_freqs[cf_index]
                    drnl=DRNLVertex(ome,CF,delays[cf_index],profile=profile)
                    g.add_machine_vertex_instance(drnl)
                    # constrain placement to local chip
                    drnl.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
                    drnls.append(drnl)
                    cf_index += 1

                    for j in range(n_ihcan):
                        ihcan = IHCANVertex(drnl,resample_factor,
                                            seeds[seed_index:seed_index+4],
                                            bitfield=bitfield,profile=profile)
                        seed_index += 4
                        g.add_machine_vertex_instance(ihcan)
                        # constrain placement to local chip
                        ihcan.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
                        ihcans.append(ihcan)
                        # Add an edge from the DRNL to the IHCAN, to send the data
                        g.add_machine_edge_instance(
                            MachineEdge(drnl, ihcan),
                            drnl.data_partition_name)

                        # Add an edge from the IHCAN to the DRNL,
                        # to send acknowledgement
                        g.add_machine_edge_instance(
                            MachineEdge(ihcan, drnl),
                            drnl.acknowledge_partition_name)

                    # Add an edge from the OME to the DRNL, to send the data
                    g.add_machine_edge_instance(
                        MachineEdge(ome, drnl),
                        ome.data_partition_name)
                if ear == 0:
                    left_ear_chips.append((chip.x, chip.y))
                count+=1

        # here we define the MC ack tree
        row_macks = drnls
        for k in range(n_mack_tree_rows):
            parents_look_up = []
            parents = []
            for mack_index,mack in enumerate(row_macks):
                parent_index = int(mack_index/num_macks)
                if parent_index not in parents_look_up:
                    if len(row_macks) <= num_macks:
                        parent = ome
                    else:
                        parent = MCackVertex()
                        g.add_machine_vertex_instance(parent)
                    parents.append(parent)
                    parents_look_up.append(parent_index)
                else:
                    parent = parents[parent_index]

                parent.register_mack_processor(mack)
                mack.register_parent_processor(parent)
                # Add an edge to send r2s data
                g.add_machine_edge_instance(
                    MachineEdge(parent, mack),
                    parent.command_partition_name)
                # Add an edge to send ack
                g.add_machine_edge_instance(
                    MachineEdge(mack, parent),
                    parent.acknowledge_partition_name)

            row_macks[:] = parents

        binaural_drnls.append(drnls)
##end for

    # Run the simulation
    if rt:
        latency = 3 * (96. / fs)
        duration = 1000.*((data.size / fs) + latency)
    else:
        latency = 3 * (10000e-6)
        duration = 1000.*((10000e-6 * (data.size / 96.)) + latency)

    g.transceiver().set_watch_dog(False)
    g.run(duration)
    # g.run_until_complete()

##begin for
    samples = list()
    for ear,drnls in enumerate(binaural_drnls):
        samples.append([])
        # Get the data back
        for drnl in drnls:
            samples[ear].append(drnl.read_samples(g.buffer_manager()))
        samples[ear] = numpy.hstack(samples[ear])

        print "ear{} channels running:{}".format(ear,len(drnls))
        print "output data: {} fibres with length {}".format(len(ihcans)*2,len(samples[ear]))
        if(len(samples[ear]) != len(ihcans)*2*numpy.floor(len(data[ear])/seg_size)*seg_size):
            print "samples length {} isn't expected size {}".format(len(samples[ear]),len(ihcans)*2*numpy.floor(len(data[ear])/seg_size)*seg_size)
##end for
    profiles = [[] for __ in range(n_ears)]
    if profile:
        for ear,ome in enumerate(binaural_omes):
            drnl_profiles=[]
            ihc_profiles=[]
            profiles[ear].append(ome._process_profile_times)
            for drnl in ome._drnl_vertices:
                drnl_profiles.append(drnl._process_profile_times)
                for ihc in drnl._ihcan_vertices:
                    ihc_profiles.append(ihc._process_profile_times)
            profiles[ear].append(numpy.asarray(drnl_profiles))
            profiles[ear].append(numpy.asarray(ihc_profiles))
    # Close the machine
    g.stop()

    return samples,numpy.asarray(profiles)
