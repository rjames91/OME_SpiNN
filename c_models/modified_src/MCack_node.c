// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/MCack_node.c

/*
 ============================================================================
 Name        : MCack_node.c
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

//#define TOTAL_TICKS 62//240//173//197
//#define PROFILE
//#define LOOP_PROFILE
//#define PRINT

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
//application initialisation
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

    log_mini_info("%u%d", 8014,num_children);  /* "num_child nodes:%d"*/
    log_mini_info("%u%d", 8015,parent_key);  /* "parent key:%d"*/
    log_mini_info("%u%d", 8016,child_key);  /* "child key:%d"*/

    mask=3;//1;
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
    uint command= mc_key & mask;

    //if (command)//ack received from child
    if (command==2)//ack received from child
    {
        sync_count++;
        log_mini_info("%u%d%d", 8017,mc_key-2,sync_count);  /* "ack mcpacket recieved from key=%d sync_count=%d"*/

        if(sync_count>=num_children)
        {
            //send acknowledgement to parent node
//            while (!spin1_send_mc_packet(parent_key+1, 0, NO_PAYLOAD)) {
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
    else if(command==1)//else if(command==0)r2s received
    {
       // log_info("r2s received from parent");
        //send r2s to child nodes
        //while (!spin1_send_mc_packet(child_key, 0, NO_PAYLOAD)) {
        while (!spin1_send_mc_packet(child_key|1, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
    }
}

void app_done ()
{
  // report simulation time
  io_printf (IO_BUF, "[core %d] simulation lasted %d ticks\n", coreID,
             spin1_get_simulation_time());

  //copy profile data
#ifdef PROFILE
  //io_printf (IO_BUF, "[core %d] saving profile data...\n", coreID);
	for (uint i=0;i<3*TOTAL_TICKS;i++)
	{
	//	profile_buffer[i]  = dtcm_profile_buffer[i];
	}
#endif
  
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

  spin1_start (SYNC_WAIT);//(SYNC_NOWAIT);//
  
  app_done ();

}

