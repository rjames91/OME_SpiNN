from spynnaker.pyNN.models.abstract_pynn_model import AbstractPyNNModel
from spinn_utilities.overrides import overrides
from spynnaker.pyNN.utilities import constants
from spinnak_ear.spinnakear_vertex import SpiNNakEarVertex
import numpy as np
import math
from spinn_utilities.executable_finder import ExecutableFinder

#audio segment size
SEG_SIZE = 8

class SpiNNakEar(AbstractPyNNModel):
    default_population_parameters = {}

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

    @overrides(AbstractPyNNModel.create_vertex)
    def create_vertex(
            self, n_neurons, label, constraints):
        max_atoms = 1

        self._vertex = SpiNNakEarVertex(
            n_neurons,  self._audio_input,self._fs,self._n_channels,
            self._pole_freqs,self._param_file,self._ear_index,
            constraints, label,max_atoms,self)

        return self._vertex
    # def save_pre_gen_vars(self,filepath):
    #     self._vertex.save_pre_gen_vars(filepath=filepath)

#naive in that it assumes the mack or an group cores won't fit onto exisitng chips (an overestimate of n_chips)
def naive_n_chips_calc(n_channels,n_ears,neuron_pops=None,n_macks=4,n_an_group=2,n_cores_per_chip=15.):
    n_mack_tree_rows = int(np.ceil(math.log(n_channels,n_macks)))
    acc=0.
    for i in range(n_mack_tree_rows-1):
        acc += n_macks**i
    mack_chips = int(np.ceil(acc/n_cores_per_chip))
    ome_drnl_ihc_chips = 1+int(np.ceil(n_channels/2.))
    n_group_tree_rows = int(np.ceil(math.log((n_channels * 10) / 2., 2.)))
    max_n_atoms_per_group_tree_row = (2. ** np.arange(1,n_group_tree_rows + 1)) * 2.
    max_n_atoms_per_group_tree_row = max_n_atoms_per_group_tree_row[max_n_atoms_per_group_tree_row < 256]
    n_group_tree_rows = max_n_atoms_per_group_tree_row.size
    for i in range(n_group_tree_rows-1):
        acc += n_an_group**i
    an_group_chips = int(np.ceil(acc/n_cores_per_chip))

    n_chips_total = n_ears*(ome_drnl_ihc_chips+mack_chips+an_group_chips)

    if neuron_pops is not None:
        neuron_cores = 0
        for (n_atoms,n_per_core) in neuron_pops:
            neuron_cores += int(np.ceil(float(n_atoms)/n_per_core))
        n_chips_total += (n_ears * int(np.ceil(neuron_cores/n_cores_per_chip)))

    return n_chips_total*2# double requested number (seems to fix issues with larger sims)




