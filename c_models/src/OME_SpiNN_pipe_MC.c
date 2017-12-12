/*
 ============================================================================
 Name        : IHC_AN_softfloat.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdfix.h>
#include "OME_SpiNN.h"
#include "spin1_api.h"
#include "math.h"
#include "complex.h"
#include "random.h"
#include "stdfix-exp.h"
#include "log.h"
#include <data_specification.h>
#include <profiler.h>
#include <profile_tags.h>


//#define TOTAL_TICKS 62//240//173//197
#define PROFILE
//#define LOOP_PROFILE
//#define PRINT

//=========GLOBAL VARIABLES============//
REAL Fs,dt,max_rate;
uint coreID;
uint chipID;
uint test_DMA;
uint seg_index;
uint cbuff_index;
uint cbuff_numseg;
uint read_switch;
uint write_switch;
uint processing;
uint index_x;
uint index_y;
uint TOTAL_TICKS;
uint TIMER_TICK_PERIOD;
uint final_ack;

uint_float_union MC_union;

REAL conchaL,conchaH,conchaG,earCanalL,earCanalH,earCanalG,stapesH,stapesL,stapesScalar,
	ARtau,ARdelay,ARrateThreshold,rateToAttenuationFactor,BFlength;

REAL concha_q,concha_j,concha_k,concha_l,conchaGainScalar,recip_conchaFilter_a0,
	earCanal_q,earCanal_j,earCanal_k,earCanal_l,earCanalGainScalar,recip_earCanalFilter_a0,
	ARatt,Wn,stapesHP_order,sf,stapesLP_b,stapes_tau,past_stapesDisp;

REAL conchaFilter_b[3],conchaFilter_a[3],earCanalFilter_b[3],earCanalFilter_a[3],stapesHP_b[3],stapesHP_a[3],stapesLP_a[2],past_input[2],past_concha[2],past_earCanalInput[2],past_earCanal[2],past_stapesInput[2]
,past_stapes[2];


//uint seed_selection[SEED_SEL_SIZE];//TODO:this needs to be moved to SDRAM

int start_count_process;
int end_count_process;
int start_count_read;
int end_count_read;
int start_count_write;
int end_count_write;
int start_count_full;
int end_count_full;

uint sync_count=0;
uint read_ticks;

REAL *dtcm_buffer_a;
REAL *dtcm_buffer_b;
REAL *dtcm_buffer_x;
REAL *dtcm_buffer_y;
REAL *dtcm_profile_buffer;

REAL *sdramin_buffer;
REAL *sdramout_buffer;
REAL *profile_buffer;

// The parameters to be read from memory
enum params {
    DATA_SIZE = 0,
    COREID,
    NUM_DRNL,
    NUM_MACK,
    FS,
    NUM_BFS,
    KEY,
    R2S_KEY,
    DATA
    };

// The size of the remaining data to be sent
uint data_size;
// Pointer to the start of the data still to be sent and confirmed
//REAL *data;
// The mask which indicates the sequence number
uint sequence_mask;
// the core ID given by the placement software
uint placement_coreID;

uint key;

uint r2s_key;

uint num_drnls;

uint num_macks;

uint sampling_frequency;

//application initialisation
void app_init(void)
{

	seg_index=0;
	cbuff_index=0;
	read_switch=0;
	write_switch=0;
	cbuff_numseg=3;
	
	/* say hello */
	
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);

	//obtain data spec

	address_t data_address = data_specification_get_data_address();
    address_t params = data_specification_get_region(0, data_address);

    // Get the size of the data in words
    data_size = params[DATA_SIZE];
    read_ticks=0;
    TOTAL_TICKS= data_size/SEGSIZE;
    log_info("data_size=%d",data_size);
    log_info("TOTAL_TICKS=%d",TOTAL_TICKS);

    // Get a pointer to the data - not worth copying at present
    sdramin_buffer = (REAL *) &(params[DATA]);

    //obtain this core ID from the host placement perspective
    placement_coreID = params[COREID];

    // Get the key to send the data with
    key = params[KEY];
    io_printf (IO_BUF, "OME-->DRNL key=%d\n",key);

    r2s_key=params[R2S_KEY];

    //Get number of child DRNL vertices
    num_drnls=params[NUM_DRNL];
    io_printf (IO_BUF, "num drnls=%d\n",num_drnls);

    //get number of macks
    num_macks=params[NUM_MACK];
    io_printf (IO_BUF, "num macks=%d\n",num_macks);

    //Get sampling frequency
    sampling_frequency = params[FS];

    Fs= (REAL)sampling_frequency;
    dt=(1.0/Fs);

    //real-time timer period calculation (us)
    if(Fs<=22050.0f)TIMER_TICK_PERIOD = (uint)(0.95*1e6 * ((REAL)SEGSIZE/Fs));//RT
    else TIMER_TICK_PERIOD = 5000;//Not RT

    log_info("timer period=%d",(uint)TIMER_TICK_PERIOD);
    log_info("Fs=%d",sampling_frequency);
    final_ack=0;

	// Allocate buffers somewhere in SDRAM

    //hack for smaller SDRAM intermediate buffers
	//data_size=cbuff_numseg*SEGSIZE;

	dtcm_buffer_a = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_b = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
/*	dtcm_buffer_x = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_y = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_profile_buffer = (REAL *) sark_alloc (3*TOTAL_TICKS, sizeof(REAL));*/
	
	//if (dtcm_buffer_a == NULL ||dtcm_buffer_b == NULL ||dtcm_buffer_x == NULL ||dtcm_buffer_y == NULL
	//		|| sdramin_buffer == NULL || dtcm_profile_buffer == NULL)
	if(dtcm_buffer_a == NULL ||dtcm_buffer_b == NULL ||sdramin_buffer == NULL)
	{
		test_DMA = FALSE;
		io_printf (IO_BUF, "[core %d] error - cannot allocate buffer\n", coreID);
	}
	else
	{
		test_DMA = TRUE;
		// initialize sections of DTCM, system RAM and SDRAM
		for (uint i = 0; i < SEGSIZE; i++)
		{
			dtcm_buffer_a[i]   = 0;
			dtcm_buffer_b[i]   = 0;
		}
	/*	for (uint i = 0; i < SEGSIZE; i++)
		{
			dtcm_buffer_x[i]   = 0;
			dtcm_buffer_y[i]   = 0;
		}

		
		for (uint i=0;i<3 * TOTAL_TICKS;i++)
		{
			dtcm_profile_buffer[i]  = 0;
			//profile_buffer[i]  = 0;
		}*/

        io_printf (IO_BUF, "[core %d] data spec output buffer tag= %d\n", coreID,
           (uint) placement_coreID);
		
	/*	io_printf (IO_BUF, "[core %d] dtcm buffer a @ 0x%08x\n", coreID,
				   (uint) dtcm_buffer_a);
		io_printf (IO_BUF, "[core %d] profile buffer @ 0x%08x\n", coreID,
				   (uint) profile_buffer);	*/

	
	//============MODEL INITIALISATION================//
    BFlength = params[NUM_BFS];
	conchaL=1500.0;
	conchaH=3000.0;
	conchaG=5.0;
	earCanalL=3000.0;
	earCanalH=3800.0;
	earCanalG=5.0;
	stapesH=700.0;
	stapesL=10.0;
	stapesScalar=5e-7;
	ARtau=0.2;
	ARdelay=0.0085;
	ARrateThreshold=100.0;
	rateToAttenuationFactor=0.1/BFlength;

	concha_q= (REAL)PI * dt * (conchaH - conchaL);
	concha_j= 1.0/(1.0+ (1.0/tan(concha_q)));
	concha_k= (2.0 * cos((REAL)PI * dt * (conchaH + conchaL)))
	            / ((1.0 + tan(concha_q)) * cos(concha_q));
	concha_l= (tan(concha_q) - 1.0)/(tan(concha_q) + 1.0);
	conchaGainScalar=pow(10.0,conchaG/20.0);

	conchaFilter_b[0]=concha_j;
	conchaFilter_b[1]=0.0;
	conchaFilter_b[2]=-1.0*concha_j;
	conchaFilter_a[0]=1.0;
	conchaFilter_a[1]=-1.0*concha_k;
	conchaFilter_a[2]=-1.0*concha_l;
	recip_conchaFilter_a0=1.0/conchaFilter_a[0];

	earCanal_q= PI * dt * (earCanalH - earCanalL);
	earCanal_j= 1.0/(1.0+ (1.0/tan(earCanal_q)));
	earCanal_k= (2.0 * cos(PI * dt * (earCanalH + earCanalL))) / ((1.0 + tan(earCanal_q)) * cos(earCanal_q));
	earCanal_l= (tan(earCanal_q) - 1.0)/(tan(earCanal_q) + 1.0);
	earCanalGainScalar=pow(10.0,earCanalG/20.0);

	earCanalFilter_b[0]=earCanal_j;
	earCanalFilter_b[1]=0.0;
	earCanalFilter_b[2]=-earCanal_j;
	earCanalFilter_a[0]=1.0;
	earCanalFilter_a[1]=-earCanal_k;
	earCanalFilter_a[2]=-earCanal_l;
	recip_earCanalFilter_a0=1.0/earCanalFilter_a[0];

	//stapes filter coeffs hard coded due to butterworth calc code overflow
	//N.B. these will need to be altered if Fs is changed!!!
	stapesHP_b[0] = 0.868412544196881;//0.931905;
	stapesHP_b[1] = -1.736825088393761;//-1.863809;
	stapesHP_b[2] = 0.868412544196881;//0.931905;
	stapesHP_a[0] = 1.0;
	stapesHP_a[1] = -1.719434219286951;//-1.859167;
	stapesHP_a[2] = 0.754215957500572;//0.868451;

	stapes_tau= 1.0/ (2 * PI * stapesL);
	stapesLP_a[0]= 1.0;
	stapesLP_a[1]= dt/stapes_tau -1.0;
	stapesLP_b= 1.0 + stapesLP_a[1];

	past_input[0]=0.0;
	past_input[1]=0.0;
	past_concha[0]=0.0;
	past_concha[1]=0.0;

	past_earCanalInput[0]=0.0;
	past_earCanalInput[1]=0.0;
	past_earCanal[0]=0.0;
	past_earCanal[1]=0.0;

	ARatt=1.0; //TODO: change this to be determined by spiking input

	past_stapesInput[0]=0.0;
	past_stapesInput[1]=0.0;
	past_stapes[0]=0.0;
	past_stapes[1]=0.0;

	past_stapesDisp=0.0;		

#ifdef PROFILE
    // configure timer 2 for profiling
    // enabled, free running, interrupt disabled, no pre-scale, 32 bit, free-running mode
    //tc[T2_CONTROL] = TIMER2_CONF;
        // Setup profiler
    profiler_init(
        data_specification_get_region(1, data_address));
#endif
	}
}
void data_write(uint null_a, uint null_b)
{
	REAL *dtcm_buffer_out;
	uint out_index;
	
	if(test_DMA == TRUE)
	{
		if(!write_switch)
		{
			out_index=index_x;
			dtcm_buffer_out=dtcm_buffer_x;
#ifdef PRINT	
			io_printf (IO_BUF, "buff_x write\n");
#endif
		}
		else
		{
			out_index=index_y;
			dtcm_buffer_out=dtcm_buffer_y;
#ifdef PRINT
			io_printf (IO_BUF, "buff_y write\n");
#endif
		}
		spin1_dma_transfer(DMA_WRITE,&sdramout_buffer[out_index],dtcm_buffer_out,DMA_WRITE,
		  						SEGSIZE*sizeof(REAL));
#ifdef PRINT
        log_info("[core %d] segment %d written to @ 0x%08x\n", coreID,seg_index,
              (uint) &sdramout_buffer[out_index]);

		io_printf (IO_BUF, "[core %d] segment %d written to @ 0x%08x - 0x%08x\n", coreID,seg_index,
							  (uint) &sdramout_buffer[out_index],(uint) &sdramout_buffer[out_index+(NUMFIBRES-1)*SEGSIZE+SEGSIZE-1]);
#endif
	}
}

void app_end(uint null_a,uint null_b)
{
    if(final_ack==1)
    {
        spin1_exit(0);
        io_printf (IO_BUF, "spinn_exit\n");
    }

    else
    {
        log_info("All data has been sent seg_index=%d",seg_index);
        //while (!spin1_send_mc_packet(r2s_key, 0, NO_PAYLOAD)) {
        while (!spin1_send_mc_packet(r2s_key|1, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
        final_ack=1;
    }

}
//DMA read
void data_read(uint null_a, uint null_b)
{
#ifdef PROFILE
    profiler_write_entry_disable_irq_fiq(PROFILER_ENTER | PROFILER_TIMER);
#endif

	REAL *dtcm_buffer_in;
	if(test_DMA == TRUE && sync_count<num_macks && seg_index==0)
	{
	    //log_info("sending r2s");
        //send ready to send MC packet
	    //while (!spin1_send_mc_packet(r2s_key, 0, NO_PAYLOAD)) {
	    while (!spin1_send_mc_packet(r2s_key|1, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
       // spin1_delay_us(1000);//to prevent unnecessarily frequent TX of r2s

	}
	//read from DMA and copy into DTCM
	else if(read_ticks<TOTAL_TICKS && test_DMA == TRUE && sync_count==num_macks)
	{
	    #ifdef PROFILE
	    if(seg_index==0)profiler_write_entry_disable_irq_fiq(PROFILER_ENTER | PROFILER_DMA_READ);
        #endif
	    //if(seg_index==0)start_count_full = tc[T2_COUNT];

	    //log_info("read_ticks=%d",read_ticks);
	    read_ticks++;
		//assign recieve buffer
		if(!read_switch)	
		{
			dtcm_buffer_in=dtcm_buffer_a;
			read_switch=1;
		}
		else
		{
			dtcm_buffer_in=dtcm_buffer_b;
			read_switch=0;
		}

		spin1_dma_transfer(DMA_READ,&sdramin_buffer[seg_index*SEGSIZE], dtcm_buffer_in, DMA_READ,
			   SEGSIZE*sizeof(REAL));
	}
	// stop if desired number of ticks reached
	else if (read_ticks >= TOTAL_TICKS && sync_count==num_macks)
	{
	    if(final_ack==1)//2nd time in this call therefore
	                    //child models have replied with finished processing acks
	    {
           // end_count_full = tc[T2_COUNT];
           // log_info("[Chip%d] total processing time=%d ticks",chipID,start_count_full-end_count_full);
            spin1_schedule_callback(app_end,NULL,NULL,2);
	    }
	    else
	    {
            sync_count=0;
	        //send finished processing MC packet
            spin1_schedule_callback(app_end,NULL,NULL,2);
	    }

 /*       log_info("All data has been sent and confirmed seg_index=%d",seg_index);
        while (!spin1_send_mc_packet(key, 1, WITH_PAYLOAD)) {
            spin1_delay_us(1);
        }
        spin1_exit(0);
        io_printf (IO_BUF, "spinn_exit\n");*/
	}
	#ifdef PROFILE
    profiler_write_entry_disable_irq_fiq(PROFILER_EXIT | PROFILER_TIMER);
    #endif
}

void process_chan(REAL *in_buffer)
{
	uint i,j,k;
	uint si=0;
	REAL concha,earCanalInput,earCanalRes,earCanalOutput,ARoutput,stapesVelocity,stapesDisplacement;
		
#ifdef PRINT
	io_printf (IO_BUF, "[core %d] segment %d (offset=%d) starting processing\n", coreID,seg_index,segment_offset);
#endif
	
	for(i=0;i<SEGSIZE;i++)
	{
		//concha (invert input gives rarefaction effect for positive pressure (?!))
		concha= (conchaFilter_b[0]*-in_buffer[i]
			+ conchaFilter_b[1]*past_input[0]
			+ conchaFilter_b[2]*past_input[1]
			- conchaFilter_a[1]*past_concha[0]
			- conchaFilter_a[2]*past_concha[1]) * recip_conchaFilter_a0;
		//update vars		
		past_input[1]=past_input[0];
		past_input[0]=-in_buffer[i];

		past_concha[1]=past_concha[0];
		past_concha[0]=concha;

		earCanalInput= conchaGainScalar* concha - in_buffer[i];

		//ear canal
		earCanalRes= (earCanalFilter_b[0]*earCanalInput
				+ earCanalFilter_b[1]*past_earCanalInput[0]
				+ earCanalFilter_b[2]*past_earCanalInput[1]
				- earCanalFilter_a[1]*past_earCanal[0]
				- earCanalFilter_a[2]*past_earCanal[1]) * recip_earCanalFilter_a0;
		//update vars	
		past_earCanalInput[1]=past_earCanalInput[0];
		past_earCanalInput[0]=earCanalInput;
	
		past_earCanal[1]=past_earCanal[0];
		past_earCanal[0]=earCanalRes;

		earCanalOutput= earCanalGainScalar * earCanalRes + earCanalInput;

		//AR
		ARoutput= ARatt * stapesScalar * earCanalOutput;

		//stapes velocity
		stapesVelocity= stapesHP_b[0] * ARoutput + 
				stapesHP_b[1] * past_stapesInput[0] + 
				stapesHP_b[2] * past_stapesInput[1] 
				- stapesHP_a[1] * past_stapes[0]
				- stapesHP_a[2] * past_stapes[1];
		//update vars
		past_stapesInput[1]= past_stapesInput[0];
		past_stapesInput[0]= ARoutput;
		past_stapes[1]= past_stapes[0];
		past_stapes[0]= stapesVelocity;

		//stapes displacement
		stapesDisplacement= stapesLP_b * stapesVelocity - stapesLP_a[1] * past_stapesDisp;
		//update vars
		past_stapesDisp=stapesDisplacement;

		MC_union.f = stapesDisplacement;
		//MC_union.f = concha;//in_buffer[i];//recip_conchaFilter_a0;//in_buffer[i];//

       // log_info("payload = %d",MC_union.u);
		//save to buffer
		//out_buffer[i]=stapesDisplacement;
		//out_buffer[i]=stapesVelocity;//nlin_b0;//nonlinout1a;//nonlinout2b;//linout1;//nonlinout1a;
        while (!spin1_send_mc_packet(key,MC_union.u, WITH_PAYLOAD)) {
            spin1_delay_us(1);
        }
	}
    //log_info("processing complete %d",seg_index);
}

void transfer_handler(uint tid, uint ttag)
{
	if (ttag==DMA_READ)
	{
#ifdef PROFILE
  //end_count_read = tc[T2_COUNT];
  //dtcm_profile_buffer[seg_index*3]=(REAL)(start_count_read-end_count_read);
#endif
  //io_printf (IO_BUF, "read complete in %d ticks\n",start_count_read-end_count_read);

		//increment segment index
		seg_index++;
		//check circular buffer
/*		if(cbuff_index<cbuff_numseg)
		{    //increment circular buffer index
		    cbuff_index++;
		}
		else
		{
		    cbuff_index=1;
		}*/
		#ifdef PROFILE
		 // start_count_process = tc[T2_COUNT];
		#endif
		
		//choose current buffers
		if(!read_switch)
		{
			process_chan(dtcm_buffer_b);
		}
		else
		{
			process_chan(dtcm_buffer_a);
		}		  


	//	spin1_trigger_user_event(NULL,NULL);
	}
	else
	{
		io_printf(IO_BUF,"[core %d] invalid %d DMA tag!\n",coreID,ttag);
	}

}

void sync_check(uint mc_key, uint null)
{
    sync_count++;
    log_info("ack mcpacket recieved from key=%d sync_count=%d seg_index=%d\n",mc_key,sync_count,seg_index);
}

void app_done ()
{
#ifdef PROFILE
  //profiler_write_entry_disable_irq_fiq(PROFILER_EXIT | PROFILER_TIMER);
  profiler_write_entry_disable_irq_fiq(PROFILER_EXIT | PROFILER_DMA_READ);
  profiler_finalise();
#endif

  // report simulation time
  io_printf (IO_BUF, "[core %d] simulation lasted %d ticks\n", coreID,
             spin1_get_simulation_time());

  //copy profile data
#ifdef PROFILE
  //io_printf (IO_BUF, "[core %d] saving profile data...\n", coreID);
	/*for (uint i=0;i<3*TOTAL_TICKS;i++)
	{
	//	profile_buffer[i]  = dtcm_profile_buffer[i];
	}*/
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

  //set timer tick
  spin1_set_timer_tick (TIMER_TICK_PERIOD);
  //setup callbacks
  //process channel once data input has been read to DTCM
  spin1_callback_on (DMA_TRANSFER_DONE,transfer_handler,1);
  //reads from DMA to DTCM every tick
  spin1_callback_on (TIMER_TICK,data_read,0);
 // spin1_callback_on (USER_EVENT,data_write,1);
  spin1_callback_on (MC_PACKET_RECEIVED,sync_check,-1);

  spin1_start (SYNC_WAIT);//(SYNC_NOWAIT);//

  
  app_done ();

}

