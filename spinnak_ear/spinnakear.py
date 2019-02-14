from spynnaker.pyNN.models.abstract_pynn_model import AbstractPyNNModel
from spinn_utilities.overrides import overrides
from spynnaker.pyNN.utilities import constants
from spinnak_ear.spinnakear_vertex import SpiNNakEarVertex
import numpy as np
from spinn_utilities.executable_finder import ExecutableFinder

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

    def __init__(self, audio_input=np.asarray([]),fs=22050.,n_channels=3000,pole_freqs=None,param_file=None,ear_index=0):
        if isinstance(audio_input,list):
            audio_input = np.asarray(audio_input)
        if len(audio_input.shape)>1:
            raise Exception("For binaural simulation please create separate SpiNNak-Ear populations (left/right)")
        self._audio_input = audio_input[0:int(np.floor(len(audio_input)/SEG_SIZE)*SEG_SIZE)]
        self._audio_input=np.asarray(self._audio_input)
        self._fs = fs
        self._n_channels = n_channels
        self._pole_freqs = pole_freqs
        self._param_file=param_file
        self._ear_index=ear_index

    @overrides(AbstractPyNNModel.create_vertex,
               additional_arguments=_population_parameters.keys())
    def create_vertex(
            self, n_neurons, label, constraints, port, tag, ip_address,
            board_address, max_on_chip_memory_usage_for_spikes_in_bytes,
            space_before_notification, spike_recorder_buffer_size,
            buffer_size_before_receive):
        max_atoms = 1

        self._vertex = SpiNNakEarVertex(
            n_neurons,  self._audio_input,self._fs,self._n_channels,self._pole_freqs,self._param_file,self._ear_index,
            port, tag, ip_address, board_address,
            max_on_chip_memory_usage_for_spikes_in_bytes,
            space_before_notification, constraints, label,
            spike_recorder_buffer_size, buffer_size_before_receive, max_atoms,
            self)

        return self._vertex

    def save_pre_gen_vars(self,filepath):
        self._vertex.save_pre_gen_vars(filepath=filepath)
