from pacman.model.graphs.application.application_vertex import ApplicationVertex
from spynnaker.pyNN.models.abstract_models.abstract_accepts_incoming_synapses import AbstractAcceptsIncomingSynapses
from pacman.model.graphs.common.constrained_object import ConstrainedObject
from pacman.model.decorators.overrides import overrides
from pacman.model.resources import ResourceContainer,SDRAMResource,CPUCyclesPerTickResource,ReverseIPtagResource,DTCMResource
from pacman.model.constraints.partitioner_constraints.fixed_vertex_atoms_constraint import FixedVertexAtomsConstraint
from pacman.model.graphs.machine.machine_edge import MachineEdge
from pacman.model.constraints.placer_constraints import SameChipAsConstraint,EarConstraint
from pacman.model.placements import Placement, Placements
from spinn_utilities.progress_bar import ProgressBar


from spinn_front_end_common.utilities import globals_variables,helpful_functions,constants
from spinn_front_end_common.utilities import constants as front_end_common_constants
from spinn_front_end_common.utilities.constants import (
    DEFAULT_BUFFER_SIZE_BEFORE_RECEIVE, MAX_SIZE_OF_BUFFERED_REGION_ON_CHIP,
    SDP_PORTS)

from pacman.executor.injection_decorator import inject_items
from data_specification.enums.data_type import DataType
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities

from spinn_front_end_common.abstract_models \
    .abstract_generates_data_specification \
    import AbstractGeneratesDataSpecification
from spinn_front_end_common.abstract_models.abstract_has_associated_binary \
    import AbstractHasAssociatedBinary
from spinn_front_end_common.abstract_models. \
    abstract_provides_outgoing_partition_constraints import \
    AbstractProvidesOutgoingPartitionConstraints
from spynnaker.pyNN.models.common.simple_population_settable \
    import SimplePopulationSettable

from spinnak_ear.OME_vertex import OMEVertex
from spinnak_ear.DRNL_vertex import DRNLVertex
from spinnak_ear.IHCAN_vertex import IHCANVertex
from spinnak_ear.MCack_vertex import MCackVertex
from spinnak_ear.AN_group_vertex import ANGroupVertex

import numpy as np
import math
import logging
logger = logging.getLogger(__name__)

# **HACK** for Projection to connect a synapse type is required
class SpiNNakEarSynapseType(object):
    def get_synapse_id_by_target(self, target):
        return 0

acknowledge_partition_dict = { 'drnl':'DRNLDataAck',
                             'ome':'OMEAck',
                             'mack':'MCackDataAck'}

data_partition_dict = {'drnl': 'DRNLData',
                     'ome':'OMEData'}

command_partition_dict = {'ome':'OMECommand',
                          'mack':'MCackData'}

#TODO: find out how we can constrain OMEs to ethernet chips

class SpiNNakEarVertex(ApplicationVertex,
                       AbstractAcceptsIncomingSynapses,
                       # ConstrainedObject,
                       # AbstractGeneratesDataSpecification,
                       # AbstractHasAssociatedBinary,
                       # AbstractProvidesOutgoingPartitionConstraints,
                       SimplePopulationSettable,
                       ):
    # The data type of each data element
    _DATA_ELEMENT_TYPE = DataType.FLOAT_64#DataType.FLOAT_32#
    # The data type of the data count
    _DATA_COUNT_TYPE = DataType.UINT32
    # The data type of the keys
    _KEY_ELEMENT_TYPE = DataType.UINT32

    SPIKE_RECORDING_REGION_ID = 0
    _N_POPULATION_RECORDING_REGIONS = 1
    _MAX_N_ATOMS_PER_CORE = 2#256
    _N_FIBRES_PER_IHCAN = 2

    def __init__(
            self, n_neurons, audio_input,fs,n_channels,pole_freqs,param_file,ear_index,
            constraints, label,max_atoms_per_core, model):
        self._model_name = "SpikeSourceSpiNNakEar"
        self._model = model
        self.param_file = param_file
        self.audio_input = audio_input
        self._data_size_bytes = (
            (self.audio_input.size * self._DATA_ELEMENT_TYPE.size) +
            self._DATA_COUNT_TYPE.size
        )
        self.fs = fs
        self.n_channels = int(n_channels)
        self._n_mack = 4 #number of mack children per parent
        self._n_ihc = 5 #number of ihcs per parent drnl
        self._sdram_resource_bytes = audio_input.dtype.itemsize * audio_input.size
        self._todo_edges = []
        self._todo_mack_reg = []
        self._ear_index=ear_index
        self._n_group_tree_rows = int(np.ceil(math.log((n_channels*10)/self._N_FIBRES_PER_IHCAN, self._MAX_N_ATOMS_PER_CORE)))
        self._max_n_atoms_per_group_tree_row = (self._MAX_N_ATOMS_PER_CORE ** np.arange(1,self._n_group_tree_rows+1)) * self._N_FIBRES_PER_IHCAN
        # self._max_n_atoms_per_group_tree_row = self._max_n_atoms_per_group_tree_row[self._max_n_atoms_per_group_tree_row < 256]
        self._max_n_atoms_per_group_tree_row = self._max_n_atoms_per_group_tree_row[self._max_n_atoms_per_group_tree_row <= 256]
        self._n_group_tree_rows = self._max_n_atoms_per_group_tree_row.size
        self._is_recording_spikes = False
        self._is_recording_moc = False

        self._mv_list = []#append to each time create_machine_vertex is called
        if pole_freqs is None:
            max_power = min([np.log10(self.fs/2.),4.25])
            self.pole_freqs = np.logspace(np.log10(30),max_power,self.n_channels)
        else:
            self.pole_freqs = pole_freqs
        self._seed_index = 0
        self._pole_index = 0

        config = globals_variables.get_simulator().config
        self._time_scale_factor = helpful_functions.read_config_int(config,"Machine","time_scale_factor")
        if self.fs / self._time_scale_factor > 22050:
            raise Exception("The input sampling frequency is too high for the chosen simulation time scale."
                            "Please reduce Fs or increase the time scale factor in the config file")
        try:
            pre_gen_vars = np.load(self.param_file)
            self._n_atoms=pre_gen_vars['n_atoms']
            self._mv_index_list=pre_gen_vars['mv_index_list']
            self._parent_index_list=pre_gen_vars['parent_index_list']
            self._edge_index_list=pre_gen_vars['edge_index_list']
            self._ihc_seeds=pre_gen_vars['ihc_seeds']
            self._ome_indices=pre_gen_vars['ome_indices']

        except:
            self._n_atoms,self._mv_index_list,self._parent_index_list,\
            self._edge_index_list,self._ihc_seeds,self._ome_indices = self.calculate_n_atoms(self.n_channels,self._n_group_tree_rows,
                                                                  n_macks=self._n_mack,n_ihcs=self._n_ihc)
            if self.param_file is not None:
                self.save_pre_gen_vars(self.param_file)

        self._size = n_neurons
        self._new_chip_indices = []
        drnl_count = 0
        for i,vertex_name in enumerate(self._mv_index_list):
            if vertex_name == "drnl":
                if drnl_count%2 == 0:
                    self._new_chip_indices.append(i)
                drnl_count+=1
        # Superclasses
        ApplicationVertex.__init__(
            self, label, constraints, max_atoms_per_core)

    def get_synapse_id_by_target(self, target):
        return 0

    def set_synapse_dynamics(self, synapse_dynamics):
        pass

    def get_maximum_delay_supported_in_ms(self, machine_time_step):
        return 1 * machine_time_step

    def add_pre_run_connection_holder(self, connection_holder, projection_edge, synapse_information):
        super(SpiNNakEarVertex, self).add_pre_run_connection_holder(connection_holder, projection_edge, synapse_information)

    def save_pre_gen_vars(self,filepath):
        np.savez_compressed(filepath,
                            n_atoms=self._n_atoms,mv_index_list=self._mv_index_list, parent_index_list= self._parent_index_list,
                            edge_index_list=self._edge_index_list,ihc_seeds=self._ihc_seeds,ome_indices=self._ome_indices)
    def record(self,variables):
        if not isinstance(variables,list):
            variables = [variables]
        if len(variables)==0:
            variables.append("all")
        for variable in variables:
            if variable == "spikes":
                self._is_recording_spikes = True
            elif variable == "moc":
                self._is_recording_moc = True
            elif variable == "all":
                self._is_recording_spikes = True
                self._is_recording_moc = True
            else:
                raise Exception("recording of " + variable + " not supported by SpiNNak-Ear!")

    def get_data(self,variables):
        b_manager = globals_variables.get_simulator().buffer_manager
        output_data = {}
        if not isinstance(variables,list):
            variables = [variables]
        if len(variables)==0:
            variables.append("all")
        if "all" in variables:
            variables = ['spikes','moc']
        for variable in variables:
            recorded_output = []
            drnl_indices = [i for i,label in enumerate(self._mv_index_list) if label=="drnl"]
            progress = ProgressBar(len(drnl_indices), "reading ear {} ".format(self._ear_index) + variable)
            for drnl in drnl_indices:
                drnl_vertex = self._mv_list[drnl]
                channel_fibres = (drnl_vertex.read_samples(b_manager,variable))
                progress.update()
                for fibre in channel_fibres:
                    recorded_output.append(fibre)
            output_data[variable]=np.asarray(recorded_output)
            progress.end()
        return output_data

    # def get_binary_start_type(self):
    #     super(SpiNNakEarVertex, self).get_binary_start_type()
    #
    # def requires_mapping(self):
    #     pass
    def get_connections_from_machine(self, transceiver, placement, edge, graph_mapper,
                                     routing_infos, synapse_information, machine_time_step):

        super(SpiNNakEarVertex, self).get_connections_from_machine(transceiver, placement, edge,
                                                           graph_mapper, routing_infos,
                                                           synapse_information,
                                                           machine_time_step)
    def clear_connection_cache(self):
        pass

    def _max_spikes_per_ts(self, n_machine_time_steps, machine_time_step):
        return 1

    @inject_items({
        "n_machine_time_steps": "TotalMachineTimeSteps",
        "machine_time_step": "MachineTimeStep"
    })
    @overrides(
        ApplicationVertex.get_resources_used_by_atoms,
        additional_arguments={"n_machine_time_steps", "machine_time_step"}
    )
    def get_resources_used_by_atoms(
            self, vertex_slice, n_machine_time_steps, machine_time_step):
        vertex_label = self._mv_index_list[vertex_slice.lo_atom]
        if vertex_label == "ome":
            sdram_resource_bytes = (9*4) + (6*8) + self._data_size_bytes
            drnl_vertices = [i for i in self._mv_index_list if i == "drnl"]
            sdram_resource_bytes += len(drnl_vertices) * self._KEY_ELEMENT_TYPE.size
            # # Live input parameters
            # reverse_iptags = [ReverseIPtagResource(
            #     port=None, sdp_port=SDP_PORTS.INPUT_BUFFERING_SDP_PORT.value,
            #     tag=None)]

        elif vertex_label == "mack":
            sdram_resource_bytes = 2*4 + 4 * self._KEY_ELEMENT_TYPE.size
            reverse_iptags = None

        elif vertex_label == "drnl":
            sdram_resource_bytes = 14*4
            sdram_resource_bytes += constants.SYSTEM_BYTES_REQUIREMENT + 8
            sdram_resource_bytes += 512 * 12#key mask tab
            sdram_resource_bytes += 8 * 8
            sdram_resource_bytes += 256 #max n bytes for conn_lut
            if self._is_recording_moc:
                sdram_resource_bytes += self._data_size_bytes
            reverse_iptags = None

        elif vertex_label == "ihc":
            if self._is_recording_spikes:
                sdram_resource_bytes = 12*4 + 1 * self._KEY_ELEMENT_TYPE.size + self._N_FIBRES_PER_IHCAN * np.ceil(self.audio_input.size/8.) * 4
            else:
                sdram_resource_bytes = 12*4 + 1 * self._KEY_ELEMENT_TYPE.size
            reverse_iptags = None
        else:#angroup
            child_vertices = [self._mv_list[vertex_index] for vertex_index in self._parent_index_list[vertex_slice.lo_atom]]
            n_child_keys = len(child_vertices)
            sdram_resource_bytes = 4*4 + 12 * n_child_keys

        container = ResourceContainer(
            sdram=SDRAMResource(
                sdram_resource_bytes +
                front_end_common_constants.SYSTEM_BYTES_REQUIREMENT),
            dtcm=DTCMResource(0),
            cpu_cycles=CPUCyclesPerTickResource(0))
            #,reverse_iptags=reverse_iptags)

        return container


    @overrides(SimplePopulationSettable.set_value)
    def set_value(self, key, value):
        SimplePopulationSettable.set_value(self, key, value)
        self._change_requires_neuron_parameters_reload = True
    def describe(self):
        """ Returns a human-readable description of the cell or synapse type.

        The output may be customised by specifying a different template\
        together with an associated template engine\
        (see ``pyNN.descriptions``).

        If template is None, then a dictionary containing the template\
        context will be returned.
        """

        parameters = dict()
        for parameter_name in self._model.default_parameters:
            parameters[parameter_name] = self.get_value(parameter_name)

        context = {
            "name": self._model_name,
            "default_parameters": self._model.default_parameters,
            "default_initial_values": self._model.default_parameters,
            "parameters": parameters,
        }
        return context
    #
    # @overrides(SimplePopulationSettable.get_value)
    # def get_value(self, key):


    @inject_items({
        "n_machine_time_steps": "TotalMachineTimeSteps",
        "machine_time_step": "MachineTimeStep"
    })
    @overrides(
        ApplicationVertex.create_machine_vertex,
        additional_arguments={"n_machine_time_steps", "machine_time_step"}
    )
    def create_machine_vertex(
            self, vertex_slice,
            resources_required,n_machine_time_steps,
            machine_time_step,  # @UnusedVariable
            label=None, constraints=None):
        #lookup relevant mv type, parent mv and edges associated with this atom
        mv_type = self._mv_index_list[vertex_slice.lo_atom]
        parent_mvs = self._parent_index_list[vertex_slice.lo_atom]
        parent_mvs.sort()#ensure lowest parent index (ome) will be first in list
        mv_edges = self._edge_index_list[vertex_slice.lo_atom]

        if mv_type == 'ome':
            vertex = OMEVertex(self.audio_input, self.fs,self.n_channels,
                               time_scale=self._time_scale_factor, profile=False)
            vertex.add_constraint(EarConstraint())

        elif mv_type == 'mack':
            vertex = MCackVertex()
            for parent in parent_mvs:
                self._mv_list[parent].register_mack_processor(vertex)
                vertex.register_parent_processor(self._mv_list[parent])
            vertex.add_constraint(EarConstraint())


        elif mv_type == 'drnl':
            for parent_index,parent in enumerate(parent_mvs):
                #first parent will be ome
                if parent_index in self._ome_indices:
                    ome = self._mv_list[parent]
                    vertex = DRNLVertex(ome,self.pole_freqs[self._pole_index],0.,is_recording=self._is_recording_moc,
                                        profile=False,drnl_index=self._pole_index)
                    self._pole_index +=1
                else:#will be a mack vertex
                    self._mv_list[parent].register_mack_processor(vertex)
                    vertex.register_parent_processor(self._mv_list[parent])
            vertex.add_constraint(EarConstraint())

        elif mv_type == 'ihc':
            for parent in parent_mvs:
                vertex = IHCANVertex(self._mv_list[parent], 1,
                                    self._ihc_seeds[self._seed_index:self._seed_index + 4],self._is_recording_spikes,
                                     ear_index=self._ear_index,bitfield=True, profile=False)
                self._seed_index += 4
                # ensure placement is on the same chip as the parent DRNL
                vertex.add_constraint(SameChipAsConstraint(self._mv_list[parent]))

        else:#an_group
            child_vertices = [self._mv_list[vertex_index] for vertex_index in parent_mvs]
            row_index = int(mv_type[-1])
            max_n_atoms = self._max_n_atoms_per_group_tree_row[row_index]

            if len(mv_edges)==0:
                #final row AN group
                is_final_row = True
            else:
                is_final_row = False

            vertex = ANGroupVertex(child_vertices,max_n_atoms = max_n_atoms,is_final_row=is_final_row)
            vertex.add_constraint(EarConstraint())

        globals_variables.get_simulator().add_machine_vertex(vertex)
        if len(mv_edges)>0:
            for (j,partition_name) in mv_edges:
                if j>vertex_slice.lo_atom: # vertex not already built
                    #add to the "to build" edge list
                    self._todo_edges.append((vertex_slice.lo_atom,j,partition_name))
                else:
                    #add edge instance
                    try:
                        globals_variables.get_simulator().add_machine_edge(MachineEdge(vertex,self._mv_list[j],label="spinnakear"),partition_name)
                    except IndexError:
                        print

        self._mv_list.append(vertex)

        if vertex_slice.lo_atom == self._n_atoms -1:
            #all vertices have been generated so add all incomplete edges
            for (source,target,partition_name) in self._todo_edges:
                globals_variables.get_simulator().add_machine_edge(MachineEdge(self._mv_list[source],self._mv_list[target],label="spinnakear"),partition_name)
        return vertex

    @property
    @overrides(ApplicationVertex.n_atoms)
    def n_atoms(self):

        return self._n_atoms

    def calculate_n_atoms(self,n_channels,n_group_tree_rows, n_macks=4, n_ihcs=5):
        # list indices correspond to atom index
        mv_index_list = []
        edge_index_list = []  # each entry is a list of tuples containnig mv indices the mv connects to and the data partion name for the edge
        parent_index_list = []  # each entry is the ID of the parent vertex - used to obtain parent spike IDs
        ome_indices = []
        ome_index = len(mv_index_list)
        ome_indices.append(ome_index)
        mv_index_list.append('ome')
        parent_index_list.append([])
        edge_index_list.append([])

        # generate MCACK tree
        mack_edge_list = list()
        mack_mv_list = list()
        mack_parent_list = list()
        n_mack_tree_rows = int(np.ceil(math.log(n_channels, n_macks)))
        row_macks = range(n_ihcs, n_channels * (n_ihcs + 1) + n_ihcs, n_ihcs + 1)
        for _ in row_macks:
            # Generate the corresponding IHC mv indices
            for j in range(n_ihcs):
                mack_mv_list.append('ihc')
                mack_edge_list.append([])
                mack_parent_list.append([])
            mack_mv_list.append('drnl')
            mack_edge_list.append([])
            mack_parent_list.append([])

        for k in range(n_mack_tree_rows):
            parents_look_up = []
            parents = []
            for mack_index, mack in enumerate(row_macks):
                parent_index = int(mack_index / n_macks)
                if len(row_macks) <= n_macks:
                    parent = ome_index
                elif parent_index not in parents_look_up:
                    parents_look_up.append(parent_index)
                    parent = len(mack_mv_list)
                    # create parent
                    mack_mv_list.append('mack')
                    mack_edge_list.append([])
                    mack_parent_list.append([])
                    parents.append(parent)
                else:
                    parent = parents[parent_index]

                try:
                    mack_parent_list[mack].append(parent)
                except IndexError:
                    print "index error"
                if parent == ome_index:
                    # add parent index to mack edge entry
                    mack_edge_list[mack].append((parent, acknowledge_partition_dict['ome']))
                else:
                    # add parent index to mack edge entry
                    mack_edge_list[mack].append((parent, acknowledge_partition_dict['mack']))

                # add mack index to ome edge entry
                edge_index_list[ome_index].append((mack, command_partition_dict['ome']))

            row_macks[:] = parents
        # reverse the order of the mack lists so it's bottom up (OME->DRNLs)
        mack_edge_list.reverse()

        mc_list_offset = len(mv_index_list)
        # reverse previously calculated indices in ome entry of the edge_index_list
        ome_edges = edge_index_list[ome_index][:]
        edge_index_list[ome_index] = []
        for (j, partition_name) in ome_edges:
            rev_index = len(mack_mv_list) - 1 - j
            edge_index_list[ome_index].append((mc_list_offset + rev_index, partition_name))
        # add all previously calculated entries to the edge_index_list
        for entry in mack_edge_list:
            edge_index_list.append([])
            for (j, partition_name) in entry:
                if partition_name == 'OMEAck':
                    edge_index_list[-1].append((ome_index, partition_name))
                else:  # reverse previously calculated indices in this list
                    rev_index = len(mack_mv_list) - 1 - j
                    rev_entry = (mc_list_offset + rev_index, partition_name)
                    edge_index_list[-1].append(rev_entry)

        mack_mv_list.reverse()
        for entry in mack_mv_list:
            mv_index_list.append(entry)

        mack_parent_list.reverse()
        for i, entry in enumerate(mack_parent_list):
            parent_index_list.append([])
            if len(entry) > 0:
                for ps in entry:
                    if ps == ome_index:
                        parent_index_list[-1].append(ps)
                    else:
                        parent_index_list[-1].append(mc_list_offset + len(mack_mv_list) - 1 - ps)

        # DRNLs should already be in the mv list so now we need to add the relevant ihc mvs
        drnl_indices = [i for i, j in enumerate(mv_index_list) if j == 'drnl']
        for i in drnl_indices:
            # Add the data edges (OME->DRNLs) to the ome entry in the edge list
            edge_index_list[ome_index].append((i, data_partition_dict['ome']))
            parent_index_list[i].append(ome_index)
            #Add the r2s edges from OME to DRNL
            edge_index_list[ome_index].append((i, command_partition_dict['ome']))
            # Generate the corresponding IHC mv indices
            for j in range(n_ihcs):
                # add the IHC mv index to the DRNL edge list entries
                edge_index_list[i].append((i + j + 1, data_partition_dict['drnl']))
                # add the DRNL index to the IHC edge list
                edge_index_list[i + j + 1].append((i, acknowledge_partition_dict['drnl']))
                # add the drnl parent index to the ihc
                parent_index_list[i + j + 1].append(i)

        # generate ihc seeds
        n_ihcans = n_channels * n_ihcs
        random_range = np.arange(n_ihcans * 4, dtype=np.uint32)
        ihc_seeds = np.random.choice(random_range, int(n_ihcans * 4), replace=False)

        # now add on the AN Group vertices
        n_child_per_group = self._MAX_N_ATOMS_PER_CORE#4#int(128./)
        n_angs = n_ihcans

        for row_index in range(n_group_tree_rows):
            n_row_angs = int(np.ceil(float(n_angs) / n_child_per_group))
            if row_index>0:
                ang_indices = [i for i, label in enumerate(mv_index_list) if label == "inter_{}".format(row_index-1)]
            else:
                ang_indices = [i for i, label in enumerate(mv_index_list) if label == "ihc"]
            for an in range(n_row_angs):
                if row_index==n_group_tree_rows-1:
                    mv_index_list.append("group_{}".format(row_index))
                else:
                    mv_index_list.append("inter_{}".format(row_index))
                edge_index_list.append([])
                ang_index = len(mv_index_list) - 1
                # find child ihcans
                child_indices = ang_indices[an * n_child_per_group:an * n_child_per_group + n_child_per_group]
                parent_index_list.append(child_indices)
                for i in child_indices:
                    edge_index_list[i].append((ang_index, 'AN'))
                #Add r2s connection from OME node
                edge_index_list[ome_index].append((ang_index, command_partition_dict['ome']))
            n_angs = n_row_angs

        return len(mv_index_list), mv_index_list, parent_index_list, edge_index_list, ihc_seeds, ome_indices