from pacman.model.graphs.machine import MachineVertex
from pacman.model.resources.resource_container import ResourceContainer
from pacman.model.resources.dtcm_resource import DTCMResource
from pacman.model.resources.sdram_resource import SDRAMResource
from pacman.model.resources.cpu_cycles_per_tick_resource \
    import CPUCyclesPerTickResource
from pacman.model.decorators.overrides import overrides
from pacman.executor.injection_decorator import inject_items

from data_specification.enums.data_type import DataType

from spinn_front_end_common.abstract_models.abstract_has_associated_binary \
    import AbstractHasAssociatedBinary
from spinn_front_end_common.abstract_models\
    .abstract_generates_data_specification \
    import AbstractGeneratesDataSpecification
from spinn_front_end_common.interface.buffer_management.buffer_models\
    .abstract_receive_buffers_to_host import AbstractReceiveBuffersToHost
from spinn_front_end_common.utilities import helpful_functions
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities
from spinn_front_end_common.utilities.utility_objs.executable_start_type \
    import ExecutableStartType

from spinn_front_end_common.abstract_models\
    .abstract_provides_n_keys_for_partition \
    import AbstractProvidesNKeysForPartition

from enum import Enum
import numpy

from spinn_machine.utilities.progress_bar import ProgressBar


class OMEVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification,
        AbstractProvidesNKeysForPartition
        ):
    """ A vertex that runs the OME algorithm
    """

    # The number of bytes for the parameters
    _N_PARAMETER_BYTES = 6*4
    # The data type of each data element
    _DATA_ELEMENT_TYPE = DataType.FLOAT_32
    # The data type of the data count
    _DATA_COUNT_TYPE = DataType.UINT32
    # The numpy data type of each data element
    _NUMPY_DATA_ELEMENT_TYPE = numpy.single
    # The data type of the keys
    _KEY_ELEMENT_TYPE = DataType.UINT32
    # the data type of the coreID
    _COREID_TYPE = DataType.UINT32

    def __init__(self, data,fs,num_bfs,data_partition_name="OMEData",
            acknowledge_partition_name="OMEDataAck",command_partition_name="OMECommand"):#TODO:add Fs to params
        """

        :param coordinator: The coordinator vertex
        :param model: The model being simulated
        """

        MachineVertex.__init__(self, label="OME Node", constraints=None)
        AbstractProvidesNKeysForPartition.__init__(self)
        self._data = data
        self._data_partition_name = data_partition_name
        self._acknowledge_partition_name = acknowledge_partition_name
        self._command_partition_name = command_partition_name
        self._fs=fs
        self._num_bfs = num_bfs

        self._drnl_vertices = list()
        self._drnl_placements = list()
        self._data_receiver = dict()
        self._mack_vertices = list()

        self._data_size = (
            (len(self._data) * self._DATA_ELEMENT_TYPE.size) +
            self._DATA_COUNT_TYPE.size
        )
        self._sdram_usage = (
            self._N_PARAMETER_BYTES + self._data_size
        )

    @property
    def data_partition_name(self):
        return self._data_partition_name

    @property
    def acknowledge_partition_name(self):
        return self._acknowledge_partition_name

    @property
    def command_partition_name(self):
        return self._command_partition_name

    def register_processor(self, drnl_vertex):
        self._drnl_vertices.append(drnl_vertex)

    def register_mack_processor(self, mack_vertex):
        self._mack_vertices.append(mack_vertex)


    def get_acknowledge_key(self, placement, routing_info):
        key = routing_info.get_first_key_from_pre_vertex(
            placement.vertex, self._acknowledge_partition_name)
        return key

    def get_mask(self, placement, routing_info):
        mask = routing_info.get_routing_info_from_pre_vertex(
            self, self._data_partition_name).first_mask
        return ~mask & 0xFFFFFFFF

    @property
    def n_data_points(self):
        return len(self._data)
    @property
    def fs(self):
        return self._fs

    @property
    def get_num_bfs(self):
        return self._num_bfs

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        sdram = self._N_PARAMETER_BYTES + self._data_size
        sdram += len(self._drnl_vertices) * self._KEY_ELEMENT_TYPE.size

        resources = ResourceContainer(
            dtcm=DTCMResource(0),
            sdram=SDRAMResource(sdram),
            cpu_cycles=CPUCyclesPerTickResource(0),
            iptags=[], reverse_iptags=[])
        return resources

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "OME.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableStartType.SYNC

    @overrides(AbstractProvidesNKeysForPartition.get_n_keys_for_partition)
    def get_n_keys_for_partition(self, partition, graph_mapper):
        return 2  # 2 for control IDs


    @inject_items({
        "routing_info": "MemoryRoutingInfos",
        "tags": "MemoryTags",
        "placements": "MemoryPlacements"
    })
    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification,
        additional_arguments=["routing_info", "tags", "placements"])
    def generate_data_specification(
            self, spec, placement, routing_info, tags, placements):

        # Reserve and write the parameters region
        region_size = self._N_PARAMETER_BYTES + self._data_size
        region_size += len(self._drnl_vertices) * self._KEY_ELEMENT_TYPE.size
        spec.reserve_memory_region(0, region_size)
        spec.switch_write_focus(0)

        # Get the placement of the vertices and find out how many chips
        # are needed
        keys = list()
        for vertex in self._drnl_vertices:
            drnl_placement = placements.get_placement_of_vertex(vertex)
            self._drnl_placements.append(drnl_placement)
            key = routing_info.get_first_key_from_pre_vertex(
                vertex, self._acknowledge_partition_name)
            keys.append(key)
        keys.sort()

        # Write the data size in words
        spec.write_value(
            len(self._data) * (float(self._DATA_ELEMENT_TYPE.size) / 4.0),
            data_type=self._DATA_COUNT_TYPE)

        # Write the CoreID
        spec.write_value(
            placement.p, data_type=self._COREID_TYPE)

     #   print "OME placement=",placement.p

        # Write number of drnls
        spec.write_value(
            len(self._drnl_vertices), data_type=self._COREID_TYPE)

        # Write number of macks
        spec.write_value(
            len(self._mack_vertices), data_type = self._COREID_TYPE)

        # Write the sampling frequency
        spec.write_value(
            self._fs, data_type=DataType.UINT32)

        spec.write_value(
            self._num_bfs, data_type=DataType.UINT32)

        # Write the key
        if len(keys)>0:
            spec.write_value(routing_info.get_routing_info_from_pre_vertex(
                self, self._data_partition_name).first_key, data_type=DataType.UINT32)
        else:
            spec.write_value(0, data_type=DataType.UINT32)

        # write the command key
        spec.write_value(routing_info.get_routing_info_from_pre_vertex(
            self,self._command_partition_name).first_key, data_type=DataType.UINT32)


        # Write the data - Arrays must be 32-bit values, so convert
        data = numpy.array(self._data, dtype=self._NUMPY_DATA_ELEMENT_TYPE)
        spec.write_array(data.view(numpy.uint32))

        # End the specification
        spec.end_specification()

    @overrides(AbstractProvidesNKeysForPartition.get_n_keys_for_partition)
    def get_n_keys_for_partition(self, partition, graph_mapper):
        return len(self._drnl_vertices)

    def read_samples(self, buffer_manager):
        """ Read back the samples
        """
        progress = ProgressBar(len(self._drnl_placements), "Reading results")
        samples = list()
        for placement in self._drnl_placements:

            # Read the data recorded
            samples.append(
                placement.vertex.read_samples(buffer_manager))
            progress.update()
        progress.end()

        # Merge all the arrays
        return numpy.hstack(samples)