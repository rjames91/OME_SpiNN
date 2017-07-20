import spinnaker_graph_front_end as g

from OME_vertex import OMEVertex
from DRNL_vertex import DRNLVertex
from IHCAN_vertex import IHCANVertex
import model_binaries

from pacman.model.constraints.placer_constraints\
    .chip_and_core_constraint import ChipAndCoreConstraint
from pacman.model.graphs.machine import MachineEdge

from spinnman.model.enums.cpu_state import CPUState

from spinn_front_end_common.utilities import globals_variables

import numpy
import logging
import time

logger = logging.getLogger(__name__)


def run_model(
        data, n_chips=None,n_drnl=0,CF=[4000],n_ihcan=0,fs=44100,resample_factor=1):
    """ Executes an MCMC model, returning the received samples

    :param data: The audio input data

    :return: nothing (output checked from ybug dump)
      """

    # Set up the simulation
    g.setup(n_chips_required=n_chips, model_binary_module=model_binaries)

    # Get the number of cores available for use
    n_cores = 0
    machine = g.machine()

    # Create a OME for each chip
    omes = dict()
    drnls = dict()
    ihcans = dict()
    boards = dict()
    for chip in machine.chips:
        #create OME
        ome=OMEVertex(data)

        g.add_machine_vertex_instance(ome)

        # constrain placement to local chip
        ome.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
        omes[chip.x, chip.y] = ome
        boards[chip.x, chip.y] = chip.ip_address

        #create DRNLs
        for i in range(n_drnl):
            drnl=DRNLVertex(ome,CF[i])
            g.add_machine_vertex_instance(drnl)
            # constrain placement to local chip
            drnl.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
            drnls[chip.x, chip.y,i] = drnl

            for j in range(n_ihcan):
                ihcan=IHCANVertex(drnl,fs,resample_factor)
                g.add_machine_vertex_instance(ihcan)
                # constrain placement to local chip
                ihcan.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
                ihcans[chip.x, chip.y,j] = ihcan
                # Add an edge from the DRNL to the IHCAN, to send the data
                g.add_machine_edge_instance(
                    MachineEdge(drnl, ihcan),
                    drnl.data_partition_name)

                # Add an edge from the IHCAN to the DRNL,
                # to send acknowledgement
                g.add_machine_edge_instance(
                    MachineEdge(ihcan, drnl),
                    drnl.acknowledge_partition_name)

            # Add an edge from the OME to the DRNL, to send the data
            g.add_machine_edge_instance(
                MachineEdge(ome, drnl),
                ome.data_partition_name)

            # Add an edge from the DRNL to the OME,
            # to send acknowledgement
            g.add_machine_edge_instance(
                MachineEdge(drnl, ome),
                ome.acknowledge_partition_name)



# Run the simulation
    g.run(None)

    # Wait for the application to finish
    txrx = g.transceiver()
    app_id = globals_variables.get_simulator()._app_id
    #logger.info("Running {} worker cores".format(n_workers))
    logger.info("Waiting for application to finish...")
    running = txrx.get_core_state_count(app_id, CPUState.RUNNING)
    while running > 0:
        time.sleep(0.5)
        error = txrx.get_core_state_count(app_id, CPUState.RUN_TIME_EXCEPTION)
        watchdog = txrx.get_core_state_count(app_id, CPUState.WATCHDOG)
        if error > 0 or watchdog > 0:
            error_msg = "Some cores have failed ({} RTE, {} WDOG)".format(
                error, watchdog)
            raise Exception(error_msg)
        running = txrx.get_core_state_count(30, CPUState.RUNNING)

    # Get the data back
    samples = list()
    '''for coordinator in coordinators.itervalues():
        samples.append(coordinator.read_samples(g.buffer_manager()))
    samples = numpy.hstack(samples)'''

    # Close the machine
    g.stop()

    return samples
