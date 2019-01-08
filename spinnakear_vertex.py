from pacman.model.graphs.application.application_vertex import ApplicationVertex
from spynnaker.pyNN.models.abstract_models.abstract_accepts_incoming_synapses import AbstractAcceptsIncomingSynapses
from pacman.model.graphs.common.constrained_object import ConstrainedObject
from pacman.model.decorators.overrides import overrides

from spinn_front_end_common.utilities import globals_variables
from spinn_front_end_common.utilities import helpful_functions

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
    n_atoms = 0
    #list indices correspond to atom index
    mv_index_list = []
    edge_index_list = []#each entry is a list of tuples containnig mv indices the mv connects to and the data partion name for the edge
    parent_index_list = [] #each entry is the ID of the parent vertex - used to obtain parent spike IDs
    for ear in range(n_ears):
        ome_index = len(mv_index_list)
        mv_index_list.append('ome')
        parent_index_list.append([])
        edge_index_list.append([])

        n_drnls = n_channels
        #generate MCACK tree
        mack_edge_list=[]
        mack_mv_list = []
        n_mack_tree_rows = int(np.ceil(math.log(n_channels, n_macks)))
        row_macks = range(n_drnls)
        for _ in row_macks:
            mack_mv_list.append('drnl')
            mack_edge_list.append([])

        mack_count=0
        for k in range(n_mack_tree_rows):
            parents_look_up = []
            parents = []
            for mack_index in row_macks:
                parent_index = int(mack_index / n_macks)
                if parent_index not in parents_look_up:
                    if len(row_macks) <= n_macks:
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
        #flip the order of the mack lists so it's bottom up (OME->DRNLs)
        for entry in mack_edge_list.reverse():
            edge_index_list.append(entry)
        for entry in mack_mv_list.reverse():
            mv_index_list.append(entry)
            parent_index_list.append([])

        #DRNLs should already be in the mv list so now we need to add the relevant ihc mvs
        drnl_indices = [i for i,j in enumerate(mv_index_list) if j=='drnl']

        for i in drnl_indices:
            # Add the data edges (OME->DRNLs) to the ome entry in the edge list
            edge_index_list[ome_index].append((i,data_partition_dict['ome']))
            parent_index_list[i]=ome_index
            #Generate the corresponding IHC mv indices
            for j in range(n_ihcs):
                mv_index_list.append('ihc')
                #add the IHC mv index to the DRNL edge list entries
                edge_index_list[i]

                #add the drnl parent index to the ihc


#TODO: find out how we can constrain OMEs to ethernet chips

class SpiNNakEarVertex(ApplicationVertex,
                       AbstractAcceptsIncomingSynapses,
                       ConstrainedObject,
                       AbstractGeneratesDataSpecification,
                       AbstractHasAssociatedBinary,
                       AbstractProvidesOutgoingPartitionConstraints,
                       SimplePopulationSettable,
                       ):
    def __init__(
            self, n_neurons, audio_input,fs,n_channels,n_ears,
            port, tag,   ip_address,board_address,
            max_on_chip_memory_usage_for_spikes_in_bytes,
            space_before_notification, constraints, label,
            spike_recorder_buffer_size, buffer_size_before_receive,
            max_atoms_per_core, model):
        self._model_name = "SpikeSourceSpiNNakEar"
        self._model = model
        self._n_atoms,self._mv_index_list,self._parent_index_list,self._edge_index_list = calculate_n_atoms(n_channels,n_ears)
        self._mv_list = []#append to each time create_machine_vertex is called

        config = globals_variables.get_simulator().config
        self._ip_address = ip_address
        if ip_address is None:
            self._ip_address = config.get("Buffers", "receive_buffer_host")
        self._port = port
        if port is None:
            self._port = helpful_functions.read_config_int(
                config, "Buffers", "receive_buffer_port")

    # **HACK** for Projection to connect a synapse type is required
    synapse_type = SpiNNakEarSynapseType()

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

    def clear_connection_cache(self):
        pass

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
        send_buffer_times = self._send_buffer_times
        if send_buffer_times is not None and len(send_buffer_times):
            if hasattr(send_buffer_times[0], "__len__"):
                send_buffer_times = send_buffer_times[
                                    vertex_slice.lo_atom:vertex_slice.hi_atom + 1]
        vertex = ReverseIPTagMulticastSourceMachineVertex(
            n_keys=vertex_slice.n_atoms,
            label=label, constraints=constraints,
            board_address=self._board_address,
            receive_port=self._receive_port,
            receive_sdp_port=self._receive_sdp_port,
            receive_tag=self._receive_tag,
            receive_rate=self._receive_rate,
            virtual_key=self._virtual_key, prefix=self._prefix,
            prefix_type=self._prefix_type, check_keys=self._check_keys,
            send_buffer_times=send_buffer_times,
            send_buffer_partition_id=self._send_buffer_partition_id,
            send_buffer_max_space=self._send_buffer_max_space,
            send_buffer_space_before_notify=(
                self._send_buffer_space_before_notify),
            buffer_notification_ip_address=(
                self._buffer_notification_ip_address),
            buffer_notification_port=self._buffer_notification_port,
            buffer_notification_tag=self._buffer_notification_tag,
            reserve_reverse_ip_tag=self._reserve_reverse_ip_tag)
        if self._record_buffer_size > 0:
            vertex.enable_recording(
                self._record_buffer_size,
                self._record_buffer_size_before_receive,
                self._record_time_between_requests)
        self._machine_vertices.append((vertex_slice, vertex))
        return vertex

    @property
    @overrides(ApplicationVertex.n_atoms)
    def n_atoms(self):

        return self._n_atoms

