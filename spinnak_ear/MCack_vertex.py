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
from spinn_front_end_common.utilities.utility_objs import ExecutableType

import numpy

from spinn_front_end_common.abstract_models\
    .abstract_provides_n_keys_for_partition \
    import AbstractProvidesNKeysForPartition

class MCackVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification,
        AbstractProvidesNKeysForPartition
        ):
    """ A vertex that runs the multi-cast acknowledge algorithm
    """
    # The number of bytes for the parameters
    _N_PARAMETER_BYTES = 2*4
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

    def __init__(self, parent=None,command_partition_name="MCackData",
            acknowledge_partition_name="MCackDataAck"):
        """

        :param ome: The connected ome vertex
        """
        AbstractProvidesNKeysForPartition.__init__(self)
        MachineVertex.__init__(self, label="MCack Node", constraints=None)

        if parent is not None:
            self.register_parent_processor(parent_vertex=parent)
        self._child_vertices = list()
        self._child_placements = list()
        self._command_partition_name = command_partition_name
        self._acknowledge_partition_name = acknowledge_partition_name

    def register_mack_processor(self, child_vertex):
        self._child_vertices.append(child_vertex)

    def register_parent_processor(self,parent_vertex):
        self._parent = parent_vertex

    def get_acknowledge_key(self, placement, routing_info):
        key = routing_info.get_first_key_from_pre_vertex(
            placement.vertex, self._acknowledge_partition_name)

        return key

    @property
    def fs(self):
        return self._fs

    @property
    def command_partition_name(self):
        return self._command_partition_name

    @property
    def acknowledge_partition_name(self):
        return self._acknowledge_partition_name

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        sdram = self._N_PARAMETER_BYTES

        resources = ResourceContainer(
            dtcm=DTCMResource(0),
            sdram=SDRAMResource(sdram),
            cpu_cycles=CPUCyclesPerTickResource(0),
            iptags=[], reverse_iptags=[])
        return resources

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "SpiNNakEar_MCACK.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableType.SYNC

    @overrides(AbstractProvidesNKeysForPartition.get_n_keys_for_partition)
    def get_n_keys_for_partition(self, partition, graph_mapper):
        return 4#for control IDs

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
        region_size = self._N_PARAMETER_BYTES
        spec.reserve_memory_region(0, region_size)
        spec.switch_write_focus(0)

        # Write the Parent key
        spec.write_value(self._parent.get_acknowledge_key(
            placement, routing_info))

        # Write number of child nodes
        spec.write_value(
            len(self._child_vertices), data_type=self._COREID_TYPE)

        # End the specification
        spec.end_specification()

