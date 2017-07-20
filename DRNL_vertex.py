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

from enum import Enum
import numpy


class DRNLVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification,
        ):
    """ A vertex that runs the DRNL algorithm
    """
    # The number of bytes for the parameters
    _N_PARAMETER_BYTES = 8*4
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

    def __init__(self, ome,CF,data_partition_name="DRNLData",
            acknowledge_partition_name="DRNLDataAck"):#TODO:add Fs to params
        """

        :param ome: The connected ome vertex
        """

        MachineVertex.__init__(self, label="DRNL Node", constraints=None)
        self._ome = ome
        self._ome.register_processor(self)
        self._CF=CF

        self._ihcan_vertices = list()
        self._ihcan_placements = list()

        self._num_data_points = ome.n_data_points
        self._data_size = (
            self._num_data_points * self._DATA_ELEMENT_TYPE.size +
            self._DATA_COUNT_TYPE.size
        )

        self._sdram_usage = (
            self._N_PARAMETER_BYTES + self._data_size
        )

        self._data_partition_name = data_partition_name
        self._acknowledge_partition_name = acknowledge_partition_name

    def register_processor(self, ihcan_vertex):
        self._ihcan_vertices.append(ihcan_vertex)


    def get_acknowledge_key(self, placement, routing_info):
        key = routing_info.get_first_key_from_pre_vertex(
            placement.vertex, self._acknowledge_partition_name)
        return key

    @property
    def n_data_points(self):
        return self._num_data_points

    @property
    def data_partition_name(self):
        return self._data_partition_name

    @property
    def acknowledge_partition_name(self):
        return self._acknowledge_partition_name


    def _get_model_parameters_array(self):
        parameters = self._model.get_parameters()
        numpy_format = list()
        numpy_values = list()
        for i, param in enumerate(parameters):
            numpy_format.append(('f{}'.format(i), param.data_type))
            numpy_values.append(param.value)
        return numpy.array(
            [tuple(numpy_values)], dtype=numpy_format).view("uint32")

    def _get_model_state_array(self):
        state = self._model.get_state_variables()
        numpy_format = list()
        numpy_values = list()
        for i, param in enumerate(state):
            numpy_format.append(('f{}'.format(i), param.data_type))
            numpy_values.append(param.initial_value)
        return numpy.array(
            [tuple(numpy_values)], dtype=numpy_format).view("uint32")

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        sdram = self._N_PARAMETER_BYTES + self._data_size
        sdram += len(self._ihcan_vertices) * self._KEY_ELEMENT_TYPE.size

        resources = ResourceContainer(
            dtcm=DTCMResource(0),
            sdram=SDRAMResource(sdram),
            cpu_cycles=CPUCyclesPerTickResource(0),
            iptags=[], reverse_iptags=[])
        return resources

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "DRNL.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableStartType.SYNC

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

        OME_placement=placements.get_placement_of_vertex(self._ome).p#how can this distinguish between multiple instances on chip?

        # Reserve and write the parameters region
        region_size = self._N_PARAMETER_BYTES
        region_size += (1 + len(self._ihcan_vertices)) * self._KEY_ELEMENT_TYPE.size
        spec.reserve_memory_region(0, region_size)
        spec.switch_write_focus(0)

        # Get the placement of the vertices and find out how many chips
        # are needed
        keys = list()
        for vertex in self._ihcan_vertices:
            ihcan_placement = placements.get_placement_of_vertex(vertex)
            self._ihcan_placements.append(ihcan_placement)
            key = routing_info.get_first_key_from_pre_vertex(
                vertex, self._acknowledge_partition_name)
            keys.append(key)
        keys.sort()

        # Write the data size in words
        spec.write_value(
            self._num_data_points * (float(self._DATA_ELEMENT_TYPE.size) / 4.0),
            data_type=self._DATA_COUNT_TYPE)

        # Write the OMECoreID
        spec.write_value(
            OME_placement, data_type=self._COREID_TYPE)

        # Write the CoreID
        spec.write_value(
            placement.p, data_type=self._COREID_TYPE)

        #Write the OMEAppID
        spec.write_value(
            0, data_type=self._COREID_TYPE)

        # Write the Acknowledge key
        spec.write_value(self._ome.get_acknowledge_key(
            placement, routing_info))

        # Write the key
        if len(keys)>0:
            routing_info = routing_info.get_routing_info_from_pre_vertex(
                self, self._data_partition_name)
            spec.write_value(routing_info.first_key, data_type=DataType.UINT32)
        else:
            spec.write_value(0, data_type=DataType.UINT32)

        # Write number of ihcans
        spec.write_value(
            len(self._ihcan_vertices), data_type=self._COREID_TYPE)

        # Write the centre frequency
        spec.write_value(self._CF,data_type=DataType.UINT32)


        print "DRNL OME placement=",OME_placement

        print "DRNL placement=",placement.p

        # End the specification
        spec.end_specification()

