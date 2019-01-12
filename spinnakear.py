from spynnaker.pyNN.models.abstract_pynn_model import AbstractPyNNModel
from spinn_utilities.overrides import overrides
from spynnaker.pyNN.utilities import constants
from spinnakear_vertex import SpiNNakEarVertex
import numpy as np

_population_parameters = {
        'port': None, 'tag': None, 'ip_address': None, 'board_address': None,
        'max_on_chip_memory_usage_for_spikes_in_bytes': (
            constants.SPIKE_BUFFER_SIZE_BUFFERING_IN),
        'space_before_notification': 640,
        'spike_recorder_buffer_size': (
            constants.EIEIO_SPIKE_BUFFER_SIZE_BUFFERING_OUT),
        'buffer_size_before_receive': (
            constants.EIEIO_BUFFER_SIZE_BEFORE_RECEIVE)
    }
#audio segment size
SEG_SIZE = 96

class SpiNNakEar(AbstractPyNNModel):
    default_population_parameters = _population_parameters

    def __init__(self, audio_input=np.asarray([]),fs=22050.,n_channels=3000):
        if isinstance(audio_input,list):
            audio_input = np.asarray(audio_input)
        self._n_ears = audio_input.shape[0]
        self._audio_input = []
        for ear in range(self._n_ears):
            self._audio_input.append(audio_input[ear][0:int(np.floor(len(audio_input[ear])/SEG_SIZE)*SEG_SIZE)])
        self._audio_input=np.asarray(self._audio_input)
        self._fs = fs
        self._n_channels = n_channels

    @overrides(AbstractPyNNModel.create_vertex,
               additional_arguments=_population_parameters.keys())
    def create_vertex(
            self, n_neurons, label, constraints, port, tag, ip_address,
            board_address, max_on_chip_memory_usage_for_spikes_in_bytes,
            space_before_notification, spike_recorder_buffer_size,
            buffer_size_before_receive):
        max_atoms = 1#self.get_max_atoms_per_core()

        return SpiNNakEarVertex(
            n_neurons,  self._audio_input,self._fs,self._n_channels,self._n_ears,
            port, tag, ip_address, board_address,
            max_on_chip_memory_usage_for_spikes_in_bytes,
            space_before_notification, constraints, label,
            spike_recorder_buffer_size, buffer_size_before_receive, max_atoms,
            self)

