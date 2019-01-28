// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/SpiNNakEar_MCack_node.c

/*
 ============================================================================
 Name        : SpiNNakEar_MCack_node.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : receives acknowledgement signals from an number of child nodes
               and passes an acknowledgement onto a parent node
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdfix.h>
#include "spin1_api.h"
#include "math.h"
#include "complex.h"
#include "random.h"
#include "stdfix-exp.h"
#include "log.h"
#include <data_specification.h>
#include <debug.h>
//=========GLOBAL VARIABLES============//
uint coreID;
uint chipID;

// The parameters to be read from memory
enum params {
    PARENT_KEY,
    CHILD_KEY,
    NUM_CHILDREN
    };

uint parent_key;
uint child_key;
uint num_children;
uint final_ack;
uint mask;
uint sync_count;

void app_init(void)
{
	/* say hello */
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);

	//obtain data spec
	address_t data_address = data_specification_get_data_address();
    address_t params = data_specification_get_region(0, data_address);

    //get parent key
    parent_key=params[PARENT_KEY];
    //get child key
    child_key=params[CHILD_KEY];
    //get number of children
    num_children=params[NUM_CHILDREN];

    log_mini_info("%u%d", 8040,num_children);  /* "num_child nodes:%d"*/
    log_mini_info("%u%d", 8041,parent_key);  /* "parent key:%d"*/
    log_mini_info("%u%d", 8042,child_key);  /* "child key:%d"*/

    mask=3;
    final_ack=0;
    sync_count=0;
}

void app_end()
{
    spin1_exit(0);
    io_printf (IO_BUF, "spinn_exit\n");
}

void sync_check(uint mc_key, uint null)
{
    sync_count++;
    log_mini_info("%u%d%d", 8043,mc_key-2,sync_count);  /* "ack mcpacket recieved from key=%d sync_count=%d"*/

    if(sync_count>=num_children)
    {
        //send acknowledgement to parent node
        while (!spin1_send_mc_packet(parent_key|2, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
        //first wave of acks
        if(!final_ack)
        {
            sync_count=0;
            final_ack=1;
        }
        //second wave of acks
        else
        {
            app_end();
        }
    }

}

void r2s_forward(uint mc_key, uint null)
{
    while (!spin1_send_mc_packet(child_key|1, 0, WITH_PAYLOAD))
    {
        spin1_delay_us(1);
    }
}

void app_done ()
{
  // report simulation time
  io_printf (IO_BUF, "[core %d] simulation lasted %d ticks\n", coreID,
             spin1_get_simulation_time());

  // say goodbye
  io_printf (IO_BUF, "[core %d] stopping simulation\n", coreID);
  io_printf (IO_BUF, "[core %d] -------------------\n", coreID);
}

void c_main()
{
  // Get core and chip IDs
  coreID = spin1_get_core_id ();
  chipID = spin1_get_chip_id ();

  app_init();

  //setup callbacks
  spin1_callback_on (MC_PACKET_RECEIVED,sync_check,-1);
  spin1_callback_on (MCPL_PACKET_RECEIVED,r2s_forward,-1);

  spin1_start (SYNC_WAIT);
  
  app_done ();

}

