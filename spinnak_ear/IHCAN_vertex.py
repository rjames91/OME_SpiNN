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
from spinn_front_end_common.utilities import helpful_functions, constants
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities
from spinn_front_end_common.utilities.utility_objs import ExecutableType

from spinn_front_end_common.abstract_models\
    .abstract_provides_n_keys_for_partition \
    import AbstractProvidesNKeysForPartition

from spinn_front_end_common.interface.profiling.abstract_has_profile_data \
    import AbstractHasProfileData
from spinn_front_end_common.interface.profiling import profile_utils
from spinn_front_end_common.interface.simulation import simulation_utilities
from spinn_front_end_common.interface.profiling.profile_data \
    import ProfileData

from spinn_front_end_common.abstract_models import AbstractRecordable
from spinn_front_end_common.utilities import globals_variables

from enum import Enum
import numpy

class IHCANVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification,
        AbstractProvidesNKeysForPartition,
        AbstractHasProfileData,
        AbstractRecordable
        ):
    """ A vertex that runs the DRNL algorithm
    """
    # The number of bytes for the parameters
    _N_PARAMETER_BYTES = 11*4#
    # The data type of each data element
    _DATA_ELEMENT_TYPE = DataType.FLOAT_32#DataType.FLOAT_64
    # The data type of the data count
    _DATA_COUNT_TYPE = DataType.UINT32
    # The numpy data type of each data element
    _NUMPY_DATA_ELEMENT_TYPE = numpy.single
    # The data type of the keys
    _KEY_ELEMENT_TYPE = DataType.UINT32
    # the data type of the coreID
    _COREID_TYPE = DataType.UINT32

    PROFILE_TAG_LABELS = {
        0: "TIMER",
        1: "DMA_READ",
        2: "INCOMING_SPIKE",
        3: "PROCESS_FIXED_SYNAPSES",
        4: "PROCESS_PLASTIC_SYNAPSES"}

    REGIONS = Enum(
        value="REGIONS",
        names=[('SYSTEM', 0),
               ('PARAMETERS', 1),
               ('RECORDING', 2),
               ('PROFILE', 3)])

    def __init__(self, drnl,resample_factor,seed,is_recording, minimum_buffer_sdram,
            buffered_sdram_per_timestep,spike_recorder,buffer_size_before_receive,
                 max_spikes_per_ts,maximum_sdram_for_buffering,ear_index=0,bitfield=True,profile=True):#TODO:add Fs to params
        """
        :param ome: The connected ome vertex    """

        MachineVertex.__init__(self, label="IHCAN Node", constraints=None)
        AbstractProvidesNKeysForPartition.__init__(self)
        AbstractRecordable.__init__(self)
        self._is_recording = is_recording
        self._minimum_buffer_sdram = minimum_buffer_sdram
        self._buffered_sdram_per_timestep = buffered_sdram_per_timestep
        self._spike_recorder = spike_recorder
        self._buffer_size_before_receive = buffer_size_before_receive
        self._max_spikes_per_ts = max_spikes_per_ts
        self._maximum_sdram_for_buffering = maximum_sdram_for_buffering
        self._ear_index = ear_index

        self._drnl = drnl
        self._drnl.register_processor(self)
        self._resample_factor=resample_factor
        self._fs=drnl.fs
        self._num_data_points = 2 * drnl.n_data_points # num of points is double previous calculations due to 2 fibre output of IHCAN model
        if bitfield:
            self._recording_size = numpy.ceil(float(self._num_data_points * self._DATA_ELEMENT_TYPE.size)/(self._DATA_ELEMENT_TYPE.size*8))
        else:
            self._recording_size = self._num_data_points * self._DATA_ELEMENT_TYPE.size
        self._seed = seed
        self._data_size = (
            self._num_data_points * self._DATA_ELEMENT_TYPE.size +
            self._DATA_COUNT_TYPE.size
        )
        self._sdram_usage = (
            self._N_PARAMETER_BYTES
        )
        # Set up for profiling
        self._n_profile_samples = 10000
        self._bitfield = bitfield
        self._profile = profile
        self._process_profile_times = None

        config = globals_variables.get_simulator().config

        self._time_between_requests = config.getint(
            "Buffers", "time_between_requests")
        self._receive_buffer_host = config.get(
            "Buffers", "receive_buffer_host")
        self._receive_buffer_port = helpful_functions.read_config_int(
            config, "Buffers", "receive_buffer_port")
        self._minimum_buffer_sdram = config.getint(
            "Buffers", "minimum_buffer_sdram")
        self._using_auto_pause_and_resume = config.getboolean(
            "Buffers", "use_auto_pause_and_resume")

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        sdram = self._N_PARAMETER_BYTES
        sdram += 1 * self._KEY_ELEMENT_TYPE.size
        # sdram += self._recording_size
        if self._profile:
            sdram += profile_utils.get_profile_region_size(self._n_profile_samples)

        resources = ResourceContainer(
            dtcm=DTCMResource(0),
            sdram=SDRAMResource(sdram),
            cpu_cycles=CPUCyclesPerTickResource(0),
            iptags=[], reverse_iptags=[])
        return resources

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "SpiNNakEar_IHCAN.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableType.USES_SIMULATION_INTERFACE

    @overrides(AbstractProvidesNKeysForPartition.get_n_keys_for_partition)
    def get_n_keys_for_partition(self, partition, graph_mapper):
        return 4#for control IDs

    @overrides(AbstractHasProfileData.get_profile_data)
    def get_profile_data(self, transceiver, placement):
        if self._profile:
            profiles = profile_utils.get_profiling_data(
                self.REGIONS.PROFILE.value,
                self.PROFILE_TAG_LABELS, transceiver, placement)
            self._process_profile_times = profiles._tags['TIMER'][1]
        else:
            profiles=ProfileData(self.PROFILE_TAG_LABELS)
        return profiles

    @inject_items({
        "routing_info": "MemoryRoutingInfos",
        "tags": "MemoryTags",
        "placements": "MemoryPlacements",
        "machine_graph":"MemoryMachineGraph",
        "n_machine_time_steps": "TotalMachineTimeSteps",
        "machine_time_step": "MachineTimeStep",
    })
    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification,
        additional_arguments=["routing_info", "tags", "placements","machine_graph"
                              ,"n_machine_time_steps","graph_mapper","machine_time_step"])
    def generate_data_specification(
            self, spec, placement, routing_info, tags, placements,machine_graph,
            n_machine_time_steps,machine_time_step):

        DRNL_placement=placements.get_placement_of_vertex(self._drnl).p

        # Setup words + 1 for flags + 1 for recording size
        setup_size = constants.SYSTEM_BYTES_REQUIREMENT + 8
        # reserve system region
        spec.reserve_memory_region(
            region=self.REGIONS.SYSTEM.value,
            size=setup_size, label='systemInfo')
        # Reserve and write the parameters region
        region_size = self._N_PARAMETER_BYTES
        region_size += 1 * self._KEY_ELEMENT_TYPE.size
        spec.reserve_memory_region(self.REGIONS.PARAMETERS.value, region_size)

        #reserve recording region
        spec.reserve_memory_region(
            self.REGIONS.RECORDING.value,
            recording_utilities.get_recording_header_size(1))
        if self._profile:
            #reserve profile region
            profile_utils.reserve_profile_region(
                spec, self.REGIONS.PROFILE.value,
                self._n_profile_samples)

        # simulation.c requirements
        spec.switch_write_focus(self.REGIONS.SYSTEM.value)
        spec.write_array(simulation_utilities.get_simulation_header_array(
            self.get_binary_file_name(), 1,
            1))

        #write parameters
        spec.switch_write_focus(self.REGIONS.PARAMETERS.value)

        # Write the data size in words
        spec.write_value(
            self._num_data_points * (float(self._DATA_ELEMENT_TYPE.size) / 4.0),
            data_type=self._DATA_COUNT_TYPE)

        # Write the DRNLCoreID
        spec.write_value(
            DRNL_placement, data_type=self._COREID_TYPE)

        # Write the CoreID
        spec.write_value(
            placement.p, data_type=self._COREID_TYPE)

        #Write the DRNLAppID
        spec.write_value(
            0, data_type=self._COREID_TYPE)

        # Write the Acknowledge key
        spec.write_value(self._drnl.get_acknowledge_key(
            placement, routing_info))

        #Write the spike resample factor
        spec.write_value(
            self._resample_factor, data_type=self._COREID_TYPE)

        #Write the sampling frequency
        spec.write_value(
            self._fs, data_type=self._COREID_TYPE)

        #Write the routing key
        partitions = machine_graph \
            .get_outgoing_edge_partitions_starting_at_vertex(self)
        for partition in partitions:
            if partition.identifier == 'SPIKE':
                rinfo = routing_info.get_routing_info_from_partition(
                    partition)
                key = rinfo.first_key
                mask = rinfo.first_mask
                spec.write_value(
                   key, data_type=self._COREID_TYPE)

        #Write the seed
        data = numpy.array(self._seed, dtype=numpy.uint32)
        spec.write_array(data.view(numpy.uint32))

        # Write the recording regions
        # spec.switch_write_focus(self.REGIONS.RECORDING.value)
        # ip_tags = tags.get_ip_tags_for_vertex(self) or []
        # spec.write_array(recording_utilities.get_recording_header_array(
        #     [self._recording_size], ip_tags=ip_tags))
        # write recording data
        ip_tags = tags.get_ip_tags_for_vertex(self) or []
        spec.switch_write_focus(self.REGIONS.RECORDING.value)
        recorded_region_sizes = recording_utilities.get_recorded_region_sizes(
            [self._spike_recorder.get_sdram_usage_in_bytes(
                2, self._max_spikes_per_ts(
                    n_machine_time_steps, machine_time_step),
                n_machine_time_steps)],
            self._maximum_sdram_for_buffering)
        spec.write_array(recording_utilities.get_recording_header_array(
            recorded_region_sizes, self._time_between_requests,
            self._buffer_size_before_receive, ip_tags))

        #Write profile regions
        if self._profile:
            profile_utils.write_profile_region_data(
                spec, self.REGIONS.PROFILE.value,
                self._n_profile_samples)

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        """ Read back the spikes """

        # Read the data recorded
        data_values, _ = buffer_manager.get_data_for_vertex(placement, 0)
        data = data_values.read_all()

        numpy_format=list()
        numpy_format.append(("AN",numpy.float32))
        if self._bitfield:
            formatted_data = numpy.array(data, dtype=numpy.uint8)
            unpacked = numpy.unpackbits(formatted_data)
            output_data = unpacked.astype(numpy.float32)
        else:
            formatted_data = numpy.array(data, dtype=numpy.uint8,copy=True).view(numpy_format)
            output_data = formatted_data.copy()


        #check all expected data has been recorded
        if len(output_data) != self._num_data_points:
            #if not set output to zeros of correct length, this will cause an error flag in run_ear.py
            #raise Warning
            print("recording not complete, reduce Fs or disable RT!\n"
                            "recorded output length:{}, expected length:{} "
                            "at placement:{},{},{}".format(len(output_data),
                            self._num_data_points,placement.x,placement.y,placement.p))

            # output_data = numpy.zeros(self._num_data_points)
            output_data.resize(self._num_data_points,refcheck=False)
        #return formatted_data
        return output_data

    # def get_n_timesteps_in_buffer_space(self, buffer_space, machine_time_step):
    #     return recording_utilities.get_n_timesteps_in_buffer_space(
    #         buffer_space, 4)
    #
    # def get_recorded_region_ids(self):
    #     return [0]

    @overrides(AbstractReceiveBuffersToHost.get_recording_region_base_address)
    def get_recording_region_base_address(self, txrx, placement):
        return helpful_functions.locate_memory_region_for_placement(
            placement, self.REGIONS.RECORDING.value, txrx)

    @overrides(AbstractRecordable.is_recording)
    def is_recording(self):
        return self._is_recording

    @overrides(AbstractReceiveBuffersToHost.get_minimum_buffer_sdram_usage)
    def get_minimum_buffer_sdram_usage(self):
        return self._minimum_buffer_sdram

    @overrides(AbstractReceiveBuffersToHost.get_n_timesteps_in_buffer_space)
    def get_n_timesteps_in_buffer_space(self, buffer_space, machine_time_step):
        return recording_utilities.get_n_timesteps_in_buffer_space(
            buffer_space, [self._buffered_sdram_per_timestep])

    @overrides(AbstractReceiveBuffersToHost.get_recorded_region_ids)
    def get_recorded_region_ids(self):
        if self._is_recording:
            return [0]
        return []