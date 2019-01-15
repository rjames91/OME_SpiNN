from pacman.model.graphs.application.application_vertex import ApplicationVertex
from spynnaker.pyNN.models.abstract_models.abstract_accepts_incoming_synapses import AbstractAcceptsIncomingSynapses
from pacman.model.graphs.common.constrained_object import ConstrainedObject
from pacman.model.decorators.overrides import overrides
from pacman.model.resources.resource_container import ResourceContainer
from pacman.model.resources.sdram_resource import SDRAMResource
from pacman.model.resources.cpu_cycles_per_tick_resource import CPUCyclesPerTickResource
from pacman.model.resources.dtcm_resource import DTCMResource
from pacman.model.constraints.partitioner_constraints.fixed_vertex_atoms_constraint import FixedVertexAtomsConstraint
from pacman.model.graphs.machine.machine_edge import MachineEdge
from pacman.model.constraints.placer_constraints import SameChipAsConstraint

from spinn_front_end_common.utilities import globals_variables,helpful_functions
from spinn_front_end_common.utilities import constants as front_end_common_constants

from data_specification.enums.data_type import DataType

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

from OME_vertex import OMEVertex
from DRNL_vertex import DRNLVertex
from IHCAN_vertex import IHCANVertex
from MCack_vertex import MCackVertex

import numpy as np
import math
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

def calculate_n_atoms(n_channels,n_ears,n_macks=4,n_ihcs=5):
    #list indices correspond to atom index
    mv_index_list = []
    edge_index_list = []#each entry is a list of tuples containnig mv indices the mv connects to and the data partion name for the edge
    parent_index_list = [] #each entry is the ID of the parent vertex - used to obtain parent spike IDs
    for ear in range(n_ears):
        ome_index = len(mv_index_list)
        mv_index_list.append('ome')
        parent_index_list.append([])
        edge_index_list.append([])

        #generate MCACK tree
        mack_edge_list=[]
        mack_mv_list = []
        mack_parent_list = []
        n_mack_tree_rows = int(np.ceil(math.log(n_channels, n_macks)))
        row_macks = range(ome_index+n_ihcs,n_channels*(n_ihcs+1)+ome_index+n_ihcs,n_ihcs+1)
        for _ in row_macks:
            #Generate the corresponding IHC mv indices
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
            for mack_index,mack in enumerate(row_macks):
                parent_index = int(mack_index / n_macks)
                if len(row_macks) <= n_macks:
                    parent = ome_index
                elif parent_index not in parents_look_up:
                    parents_look_up.append(parent_index)
                    parent = ome_index+len(mack_mv_list)
                    #create parent
                    mack_mv_list.append('mack')
                    mack_edge_list.append([])
                    mack_parent_list.append([])
                    parents.append(parent)
                else:
                    parent = parents[parent_index]

                # mack_parent_list[mack]=parent #will overwrite empty list initial value
                mack_parent_list[mack].append(parent)
                if parent == ome_index:
                    # add parent index to mack edge entry
                    mack_edge_list[mack].append((parent, acknowledge_partition_dict['ome']))
                    # add mack index to parent edge entry
                    edge_index_list[parent].append((mack, command_partition_dict['ome']))
                else:
                    # add parent index to mack edge entry
                    mack_edge_list[mack].append((parent, acknowledge_partition_dict['mack']))
                    # add mack index to parent edge entry
                    mack_edge_list[parent-ome_index].append((mack, command_partition_dict['mack']))


            row_macks[:] = parents
        #reverse the order of the mack lists so it's bottom up (OME->DRNLs)
        mack_edge_list.reverse()

        mc_list_offset = len(mv_index_list)
        #reverse previously calculated indices in ome entry of the edge_index_list
        ome_edges = edge_index_list[ome_index][:]
        edge_index_list[ome_index]=[]
        for(j,partition_name) in ome_edges:
            rev_index = len(mack_mv_list)-1-j
            edge_index_list[ome_index].append((mc_list_offset+rev_index,partition_name))
        #add all previously calculated entries to the edge_index_list
        for entry in mack_edge_list:
            edge_index_list.append([])
            for (j,partition_name) in entry:
                if partition_name == 'OMEAck':
                    edge_index_list[-1].append((ome_index,partition_name))
                else:#reverse previously calculated indices in this list
                    rev_index = len(mack_mv_list)-1-j
                    rev_entry = (mc_list_offset+rev_index,partition_name)
                    edge_index_list[-1].append(rev_entry)

        mack_mv_list.reverse()
        for entry in mack_mv_list:
            mv_index_list.append(entry)

        mack_parent_list.reverse()
        for i,entry in enumerate(mack_parent_list):
            parent_index_list.append([])
            if len(entry)>0:
                for ps in entry:
                    if ps == ome_index:
                        parent_index_list[-1].append(ps)
                    else:
                        parent_index_list[-1].append(mc_list_offset+len(mack_mv_list)-1-ps)

        #DRNLs should already be in the mv list so now we need to add the relevant ihc mvs
        drnl_indices = [i for i,j in enumerate(mv_index_list) if j=='drnl']

        ihc_index = len(mv_index_list)
        for i in drnl_indices:
            # Add the data edges (OME->DRNLs) to the ome entry in the edge list
            edge_index_list[ome_index].append((i,data_partition_dict['ome']))
            parent_index_list[i].append(ome_index)
            #Generate the corresponding IHC mv indices
            for j in range(n_ihcs):
                #add the IHC mv index to the DRNL edge list entries
                edge_index_list[i].append((i+j+1,data_partition_dict['drnl']))
                #add the DRNL index to the IHC edge list
                edge_index_list[i+j+1].append((i,acknowledge_partition_dict['drnl']))
                #add the drnl parent index to the ihc
                parent_index_list[i+j+1].append(i)

    #generate ihc seeds
    n_ihcans = n_channels * n_ihcs
    random_range = np.arange(n_ihcans * 4, dtype=np.uint32)
    ihc_seeds = np.random.choice(random_range, int(n_ihcans * 4), replace=False)
    return len(mv_index_list),mv_index_list,parent_index_list,edge_index_list,ihc_seeds

#TODO: find out how we can constrain OMEs to ethernet chips

class SpiNNakEarVertex(ApplicationVertex,
                       AbstractAcceptsIncomingSynapses,
                       # ConstrainedObject,
                       # AbstractGeneratesDataSpecification,
                       # AbstractHasAssociatedBinary,
                       # AbstractProvidesOutgoingPartitionConstraints,
                       # SimplePopulationSettable,
                       ):
    # The data type of each data element
    _DATA_ELEMENT_TYPE = DataType.FLOAT_64#DataType.FLOAT_32#
    # The data type of the data count
    _DATA_COUNT_TYPE = DataType.UINT32
    # The data type of the keys
    _KEY_ELEMENT_TYPE = DataType.UINT32
    def __init__(
            self, n_neurons, audio_input,fs,n_channels,n_ears,
            port, tag,   ip_address,board_address,
            max_on_chip_memory_usage_for_spikes_in_bytes,
            space_before_notification, constraints, label,
            spike_recorder_buffer_size, buffer_size_before_receive,
            max_atoms_per_core, model):
        self._model_name = "SpikeSourceSpiNNakEar"
        self._model = model
        self._audio_input = audio_input
        self._data_size = (
            (self._audio_input.size * self._DATA_ELEMENT_TYPE.size) +
            self._DATA_COUNT_TYPE.size
        )
        self._fs = fs
        self._n_channels = n_channels
        self._n_ears = n_ears
        self._n_mack = 4 #number of mack children per parent
        self._n_ihc = 5 #number of ihcs per parent drnl
        self._ear_index = 0
        self._sdram_resource_bytes = audio_input.dtype.itemsize * audio_input.size
        self._todo_edges = []
        self._todo_mack_reg = []

        self._mv_list = []#append to each time create_machine_vertex is called
        max_power = min([np.log10(self._fs/2.),4.25])
        self._pole_freqs = np.logspace(np.log10(30),max_power,self._n_channels)
        self._seed_index = 0
        self._pole_index = 0

        config = globals_variables.get_simulator().config
        self._ip_address = ip_address
        if ip_address is None:
            self._ip_address = config.get("Buffers", "receive_buffer_host")
        self._port = port
        if port is None:
            self._port = helpful_functions.read_config_int(
                config, "Buffers", "receive_buffer_port")
        self._time_scale_factor = helpful_functions.read_config_int(config,"Machine","time_scale_factor")
        if self._fs / self._time_scale_factor > 22050:
            raise Exception("The input sampling frequency is too high for the chosen simulation time scale."
                            "Please reduce Fs or increase the time scale factor in the config file")
        self._n_atoms,self._mv_index_list,self._parent_index_list,\
        self._edge_index_list,self._ihc_seeds = calculate_n_atoms(n_channels,n_ears,
                                                                  n_macks=self._n_mack,n_ihcs=self._n_ihc)
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

    # **HACK** for Projection to connect a synapse type is required
    # synapse_type = SpiNNakEarSynapseType()
    def get_synapse_id_by_target(self, target):
        return 0

    def set_synapse_dynamics(self, synapse_dynamics):
        pass

    def get_maximum_delay_supported_in_ms(self, machine_time_step):
        return 1 * machine_time_step

    def add_pre_run_connection_holder(self, connection_holder, projection_edge, synapse_information):
        super(SpiNNakEarVertex, self).add_pre_run_connection_holder(connection_holder, projection_edge, synapse_information)

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

    @overrides(ApplicationVertex.get_resources_used_by_atoms)
    def get_resources_used_by_atoms(self, vertex_slice):
        # **HACK** only way to force no partitioning is to zero dtcm and cpu
        vertex_label = self._mv_index_list[vertex_slice.lo_atom]
        if vertex_label == "ome":
            sdram_resource_bytes = (9*4) + (6*8) + self._data_size
            drnl_vertices = [i for i in self._mv_index_list if i == "drnl"]
            sdram_resource_bytes += len(drnl_vertices) * self._KEY_ELEMENT_TYPE.size
        elif vertex_label == "mack":
            sdram_resource_bytes = 3*4
        elif vertex_label == "drnl":
            sdram_resource_bytes = 10*4
            sdram_resource_bytes += self._n_ihc * self._KEY_ELEMENT_TYPE.size

        else:#ihc
            sdram_resource_bytes = 10*4 + 1 * self._KEY_ELEMENT_TYPE.size + self._data_size

        container = ResourceContainer(
            sdram=SDRAMResource(
                sdram_resource_bytes +
                front_end_common_constants.SYSTEM_BYTES_REQUIREMENT),
            dtcm=DTCMResource(0),
            cpu_cycles=CPUCyclesPerTickResource(0))

        return container
    # ------------------------------------------------------------------------
    # AbstractHasAssociatedBinary overrides TODO: allow correct binary and start type to be chosen depending on slice atom index
    # ------------------------------------------------------------------------
    # @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    # def get_binary_file_name(self):
    #     return "breakout.aplx"
    #
    # @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    # def get_binary_start_type(self):
    #     # return ExecutableStartType.USES_SIMULATION_INTERFACE
    #     return ExecutableType.USES_SIMULATION_INTERFACE

    # ------------------------------------------------------------------------
    # AbstractProvidesOutgoingPartitionConstraints overrides
    # ------------------------------------------------------------------------
    # @overrides(AbstractProvidesOutgoingPartitionConstraints.
    #            get_outgoing_partition_constraints)
    # def get_outgoing_partition_constraints(self, partition):
    #     return [ContiguousKeyRangeContraint()]
    #
    # @property
    # @overrides(AbstractChangableAfterRun.requires_mapping)
    # def requires_mapping(self):
    #     return self._change_requires_mapping
    #
    # @overrides(AbstractChangableAfterRun.mark_no_changes)
    # def mark_no_changes(self):
    #     self._change_requires_mapping = False

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

    @overrides(ApplicationVertex.create_machine_vertex)
    def create_machine_vertex(
            self, vertex_slice,
            resources_required,  # @UnusedVariable
            label=None, constraints=None):
        # send_buffer_times = self._send_buffer_times
        # if send_buffer_times is not None and len(send_buffer_times):
        #     if hasattr(send_buffer_times[0], "__len__"):
        #         send_buffer_times = send_buffer_times[
        #                             vertex_slice.lo_atom:vertex_slice.hi_atom + 1]
        #lookup relevant mv type, parent mv and edges associated with this atom
        mv_type = self._mv_index_list[vertex_slice.lo_atom]
        parent_mvs = self._parent_index_list[vertex_slice.lo_atom]
        parent_mvs.sort()#ensure lowest parent index (ome) will be first in list
        mv_edges = self._edge_index_list[vertex_slice.lo_atom]

        if mv_type == 'ome':
            vertex = OMEVertex(self._audio_input[self._ear_index], self._fs,self._n_channels,
                               time_scale=self._time_scale_factor, profile=False)
            self._ear_index +=1

        elif mv_type == 'mack':
            vertex = MCackVertex()
            for parent in parent_mvs:
                self._mv_list[parent].register_mack_processor(vertex)
                vertex.register_parent_processor(self._mv_list[parent])
        elif mv_type == 'drnl':
            for parent_index,parent in enumerate(parent_mvs):
                #first parent will be ome
                if parent_index == 0:
                    ome = self._mv_list[parent]
                    vertex = DRNLVertex(ome,self._pole_freqs[self._pole_index],0.,profile=False)
                    self._pole_index +=1
                else:#will be a mack vertex
                    self._mv_list[parent].register_mack_processor(vertex)
                    vertex.register_parent_processor(self._mv_list[parent])
        else:#ihcans
            for parent in parent_mvs:
                vertex = IHCANVertex(self._mv_list[parent], 1,
                                    self._ihc_seeds[self._seed_index:self._seed_index + 4],
                                    bitfield=True, profile=False)
                self._seed_index += 4
                # ensure placement is on the same chip as the parent DRNL
                vertex.add_constraint(SameChipAsConstraint(self._mv_list[parent]))

        globals_variables.get_simulator().add_machine_vertex(vertex)
        for (j,partition_name) in mv_edges:
            if j>vertex_slice.lo_atom: # vertex not already built
                #add to the "to build" edge list
                self._todo_edges.append((vertex_slice.lo_atom,j,partition_name))
            else:
                #add edge instance
                globals_variables.get_simulator().add_machine_edge(MachineEdge(vertex,self._mv_list[j],label="spinnakear"),partition_name)

        # vertex = ReverseIPTagMulticastSourceMachineVertex(
        #     n_keys=vertex_slice.n_atoms,
        #     label=label, constraints=constraints,
        #     board_address=self._board_address,
        #     receive_port=self._receive_port,
        #     receive_sdp_port=self._receive_sdp_port,
        #     receive_tag=self._receive_tag,
        #     receive_rate=self._receive_rate,
        #     virtual_key=self._virtual_key, prefix=self._prefix,
        #     prefix_type=self._prefix_type, check_keys=self._check_keys,
        #     send_buffer_times=send_buffer_times,
        #     send_buffer_partition_id=self._send_buffer_partition_id,
        #     send_buffer_max_space=self._send_buffer_max_space,
        #     send_buffer_space_before_notify=(
        #         self._send_buffer_space_before_notify),
        #     buffer_notification_ip_address=(
        #         self._buffer_notification_ip_address),
        #     buffer_notification_port=self._buffer_notification_port,
        #     buffer_notification_tag=self._buffer_notification_tag,
        #     reserve_reverse_ip_tag=self._reserve_reverse_ip_tag)
        # if self._record_buffer_size > 0:
        #     vertex.enable_recording(
        #         self._record_buffer_size,
        #         self._record_buffer_size_before_receive,
        #         self._record_time_between_requests)
        # self._machine_vertices.append((vertex_slice, vertex))
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

