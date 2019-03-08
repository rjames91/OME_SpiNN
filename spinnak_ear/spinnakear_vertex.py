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
import random
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
    _N_LSR_PER_IHC = 2#3
    _N_MSR_PER_IHC = 2#3
    _N_HSR_PER_IHC = 6#4
    _N_FIBRES_PER_IHC = _N_LSR_PER_IHC + _N_MSR_PER_IHC + _N_HSR_PER_IHC

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
        self._n_ihc = 5 #number of ihcs per parent drnl

        self._sdram_resource_bytes = audio_input.dtype.itemsize * audio_input.size
        self._todo_edges = []
        self._todo_mack_reg = []
        self._ear_index=ear_index
        self._n_group_tree_rows = int(np.ceil(math.log((n_channels*self._N_FIBRES_PER_IHC)/self._N_FIBRES_PER_IHCAN, self._MAX_N_ATOMS_PER_CORE)))
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
                                                                      n_ihcs=self._n_ihc)
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
            sdram_resource_bytes += 512 * 12#key mask tab
            sdram_resource_bytes += 8 * 8
            sdram_resource_bytes += 256 #max n bytes for conn_lut
            if self._is_recording_moc:
                sdram_resource_bytes += self._data_size_bytes
            reverse_iptags = None

        # elif vertex_label == "ihc":
        elif "ihc" in vertex_label:
            if self._is_recording_spikes:
                sdram_resource_bytes = 15*4 + 1 * self._KEY_ELEMENT_TYPE.size + self._N_FIBRES_PER_IHCAN * np.ceil(self.audio_input.size/8.) * 4
            else:
                sdram_resource_bytes = 15*4 + 1 * self._KEY_ELEMENT_TYPE.size
            reverse_iptags = None
        else:#angroup
            child_vertices = [self._mv_list[vertex_index] for vertex_index in self._parent_index_list[vertex_slice.lo_atom]]
            n_child_keys = len(child_vertices)
            sdram_resource_bytes = 5*4 + 12 * n_child_keys

        container = ResourceContainer(
            sdram=SDRAMResource(
                sdram_resource_bytes +
                front_end_common_constants.SYSTEM_BYTES_REQUIREMENT)+ 8,
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

        elif 'ihc' in mv_type:#mv_type == 'ihc':
            for parent in parent_mvs:
                n_lsr = int(mv_type[-3])
                n_msr = int(mv_type[-2])
                n_hsr = int(mv_type[-1])
                vertex = IHCANVertex(self._mv_list[parent], 1,
                                    self._ihc_seeds[self._seed_index:self._seed_index + 4],self._is_recording_spikes,
                                     ear_index=self._ear_index,bitfield=True, profile=False,n_lsr=n_lsr,n_msr=n_msr,n_hsr=n_hsr)
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

    def calculate_n_atoms(self,n_channels,n_group_tree_rows,n_ihcs=5):
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

        for _ in range(n_channels):
            drnl_index = len(mv_index_list)
            mv_index_list.append('drnl')
            parent_index_list.append([ome_index])
            edge_index_list.append([])
            #OME command
            # edge_index_list[ome_index].append((drnl_index, command_partition_dict['ome']))
            # Add the data edges (OME->DRNLs) to the ome entry in the edge list
            edge_index_list[ome_index].append((drnl_index, data_partition_dict['ome']))
            fibres = []
            for _ in range(self._N_HSR_PER_IHC):
                fibres.append(2)
            for _ in range(self._N_MSR_PER_IHC):
                fibres.append(1)
            for _ in range(self._N_LSR_PER_IHC):
                fibres.append(0)

            random.shuffle(fibres)

            for j in range(n_ihcs):
                ihc_index = len(mv_index_list)
                #randomly pick fibre types
                chosen_indices=[fibres.pop() for _ in range(self._N_FIBRES_PER_IHCAN)]
                mv_index_list.append('ihc{}{}{}'.format(chosen_indices.count(0),chosen_indices.count(1),chosen_indices.count(2)))
                #drnl data/command
                # add the IHC mv index to the DRNL edge list entries
                edge_index_list[drnl_index].append((ihc_index, data_partition_dict['drnl']))
                # add the drnl parent index to the ihc
                parent_index_list.append([drnl_index])
                edge_index_list.append([])

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
                ang_indices = [i for i, label in enumerate(mv_index_list) if "ihc" in label]
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
                # edge_index_list[ome_index].append((ang_index, command_partition_dict['ome']))
            n_angs = n_row_angs

        return len(mv_index_list), mv_index_list, parent_index_list, edge_index_list, ihc_seeds, ome_indices