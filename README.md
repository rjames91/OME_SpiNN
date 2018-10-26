# OME_SpiNN

repository needed to launch SpiNNakEar simulation using https://github.com/rjames91/OME_SpiNN/blob/master/run_ear.py

n.b. C code of DRNL and IHCAN models can be found in https://github.com/rjames91/DRNL_SpiNN and https://github.com/rjames91/IHC_AN_SpiNN
To build SpiNNakEar binaries the "-ffrestanding" compiler flags must be removed from spinnaker_tools.mk

general spinnaker plotting functions etc used in these simulation scripts can be found in https://github.com/SpiNNakerManchester/PyNN8Examples/blob/rob_testing/signal_prep.py

N.B. When generating audio stimulus files please refrain from using long duration 0 value samples. This will lead to an increase in processing time and compromise real-time performance.
    The generate_signal function in signal_prep.py models silence as a low amplitude (-20dBSPL) background noise.