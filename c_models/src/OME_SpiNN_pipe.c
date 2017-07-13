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
#define TIMER_TICK_PERIOD 2600//2300//REALTIME (2.3ms to process 100 44100Hz samples

#define TOTAL_TICKS 240//173//197       
#define PROFILE
//#define LOOP_PROFILE
//#define PRINT

//=========GLOBAL VARIABLES============//
REAL Fs,dt,max_rate;
uint coreID;
uint chipID;
uint test_DMA;
uint seg_index;
uint read_switch;
uint write_switch;
uint processing;
uint index_x;
uint index_y;

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
   // SEQUENCE_MASK,
    COREID,
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

//application initialisation
void app_init(void)
{
	Fs=SAMPLING_FREQUENCY;
	dt=(1.0/Fs);
	seg_index=0;
	read_switch=0;
	write_switch=0;
	
	/* say hello */
	
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);

	//obtain data spec

	address_t data_address = data_specification_get_data_address();
    address_t params = data_specification_get_region(0, data_address);

    // Get the size of the data in words
    data_size = params[DATA_SIZE];
    // Get a pointer to the data - not worth copying at present
    sdramin_buffer = (REAL *) &(params[DATA]);

    //obtain this core ID from the host placement perspective
    placement_coreID = params[COREID];

	// Allocate buffers somewhere in SDRAM

	
	//output results buffer
	sdramout_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					 data_size * sizeof(REAL),
					 placement_coreID,
					 ALLOC_LOCK);

	/*sdramin_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					MAX_SIGNAL_S*(uint)44100. *sizeof(REAL),
					coreID|32,
					ALLOC_LOCK);*/
	
	profile_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					3 * ((uint)44100./SEGSIZE) *sizeof(REAL),
					coreID|64,
					ALLOC_LOCK);
	
	// and a buffer in DTCM
	
	dtcm_buffer_a = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_b = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_x = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_y = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_profile_buffer = (REAL *) sark_alloc (3*TOTAL_TICKS, sizeof(REAL));
	
	if (dtcm_buffer_a == NULL ||dtcm_buffer_b == NULL ||dtcm_buffer_x == NULL ||dtcm_buffer_y == NULL 
			||  sdramout_buffer == NULL || sdramin_buffer == NULL || profile_buffer == NULL || dtcm_profile_buffer == NULL)
	/*if (sdramout_buffer == NULL || sdramin_buffer == NULL || dtcm_buffer_y == NULL 
			|| dtcm_buffer_a == NULL || dtcm_buffer_b == NULL || dtcm_profile_buffer == NULL ||dtcm_buffer_x == NULL)*/
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
		for (uint i = 0; i < SEGSIZE; i++)
		{
			dtcm_buffer_x[i]   = 0;
			dtcm_buffer_y[i]   = 0;
		}
	
		for (uint i=0;i<MAX_SIGNAL_S* ((uint)44100.);i++)
		{
			sdramout_buffer[i]  = 0;
		}
		for (uint i=0;i<MAX_SIGNAL_S * (uint)44100.;i++)
		{
			sdramin_buffer[i]  = 0;
		}
		
		for (uint i=0;i<3 * TOTAL_TICKS;i++)
		{
			dtcm_profile_buffer[i]  = 0;
			profile_buffer[i]  = 0;
		}
		
		io_printf (IO_BUF, "[core %d] dtcm buffer a @ 0x%08x\n", coreID,
				   (uint) dtcm_buffer_a);
		io_printf (IO_BUF, "[core %d] sdram out buffer @ 0x%08x\n", coreID,
				   (uint) sdramout_buffer);
		io_printf (IO_BUF, "[core %d] sdram in buffer @ 0x%08x\n", coreID,
				   (uint) sdramin_buffer);	
		io_printf (IO_BUF, "[core %d] profile buffer @ 0x%08x\n", coreID,
				   (uint) profile_buffer);	
	}
	
	//============MODEL INITIALISATION================//
	BFlength=1.0;//TODO change this to input parameter

	conchaL=1500.0;
	conchaH=3000.0;
	conchaG=5.0;
	earCanalL=3000.0;
	earCanalH=3800;
	earCanalG=5.0;
	stapesH=700.0;
	stapesL=10.0;
	stapesScalar=5e-7;
	ARtau=0.2;
	ARdelay=0.0085;
	ARrateThreshold=100.0;
	rateToAttenuationFactor=0.1/BFlength;

	concha_q= PI * dt * (conchaH - conchaL);
	concha_j= 1.0/(1.0+ (1.0/tan(concha_q)));
	concha_k= (2.0 * cos(PI * dt * (conchaH + conchaL))) / ((1.0 + tan(concha_q)) * cos(concha_q));
	concha_l= (tan(concha_q) - 1.0)/(tan(concha_q) + 1.0);
	conchaGainScalar=pow(10.0,conchaG/20.0);

	conchaFilter_b[0]=concha_j;
	conchaFilter_b[1]=0.0;
	conchaFilter_b[2]=-concha_j;
	conchaFilter_a[0]=1.0;
	conchaFilter_a[1]=-concha_k;
	conchaFilter_a[2]=-concha_l;
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
	stapesHP_b[0] = 0.931905;
	stapesHP_b[1] = -1.863809;
	stapesHP_b[2] = 0.931905;
	stapesHP_a[0] = 1.0;
	stapesHP_a[1] = -1.859167;
	stapesHP_a[2] = 0.868451;

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
    tc[T2_CONTROL] = TIMER2_CONF;
#endif
    
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
#ifdef PROFILE
  start_count_write = tc[T2_COUNT];
#endif
		spin1_dma_transfer(DMA_WRITE,&sdramout_buffer[out_index],dtcm_buffer_out,DMA_WRITE,
		  						SEGSIZE*sizeof(REAL));
#ifdef PRINT
		io_printf (IO_BUF, "[core %d] segment %d written to @ 0x%08x - 0x%08x\n", coreID,seg_index,
							  (uint) &sdramout_buffer[out_index],(uint) &sdramout_buffer[out_index+(NUMFIBRES-1)*SEGSIZE+SEGSIZE-1]);
#endif
	}
}

//DMA read
void data_read(uint ticks, uint null)
{
#ifdef PROFILE
  start_count_read = tc[T2_COUNT];
#endif

	REAL *dtcm_buffer_in;

	//read from DMA and copy into DTCM
	if(test_DMA == TRUE)
	{
		//assign recieve buffer
		if(!read_switch)	
		{
			dtcm_buffer_in=dtcm_buffer_a;
			read_switch=1;
#ifdef PRINT
			io_printf (IO_BUF, "buff_a read\n");
#endif
		}
		else
		{
			dtcm_buffer_in=dtcm_buffer_b;
			read_switch=0;
#ifdef PRINT
			io_printf (IO_BUF, "buff_b read\n");
#endif
		}

#ifdef PRINT
		io_printf (IO_BUF, "[core %d] sdram DMA read @ 0x%08x (segment %d)\n", coreID,
					  (uint) &sdramin_buffer[(seg_index)*SEGSIZE],seg_index+1);
#endif
		
		spin1_dma_transfer(DMA_READ,&sdramin_buffer[seg_index*SEGSIZE], dtcm_buffer_in, DMA_READ,
			   SEGSIZE*sizeof(REAL));
	}
	
	// stop if desired number of ticks reached
	if (ticks > TOTAL_TICKS) 
	{
		io_printf (IO_BUF, "spinn_exit\n");
		spin1_exit (0); 
	}
	
}


uint process_chan(REAL *out_buffer,REAL *in_buffer) 
{  
	uint segment_offset=SEGSIZE*(seg_index-1);
	uint i,j,k;
		
	uint si=0;
	REAL concha,earCanalInput,earCanalRes,earCanalOutput,ARoutput,stapesVelocity,stapesDisplacement;
		
#ifdef PRINT
	io_printf (IO_BUF, "[core %d] segment %d (offset=%d) starting processing\n", coreID,seg_index,segment_offset);
#endif
	
	for(i=0;i<SEGSIZE;i++)
	{
		//concha
		concha= (conchaFilter_b[0]*in_buffer[i]
			+ conchaFilter_b[1]*past_input[0]
			+ conchaFilter_b[2]*past_input[1]
			- conchaFilter_a[1]*past_concha[0]
			- conchaFilter_a[2]*past_concha[1]) * recip_conchaFilter_a0;
		//update vars		
		past_input[1]=past_input[0];
		past_input[0]=in_buffer[i];

		past_concha[1]=past_concha[0];
		past_concha[0]=concha;

		earCanalInput= conchaGainScalar* concha + in_buffer[i];

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
	
		//save to buffer
		out_buffer[i]=stapesDisplacement;
		//out_buffer[i]=stapesVelocity;//nlin_b0;//nonlinout1a;//nonlinout2b;//linout1;//nonlinout1a;
	}
		
	return segment_offset;
}

void transfer_handler(uint tid, uint ttag)
{
	if (ttag==DMA_READ)
	{
#ifdef PROFILE
  end_count_read = tc[T2_COUNT];
  dtcm_profile_buffer[seg_index*3]=(REAL)(start_count_read-end_count_read);
#ifdef PRINT
  io_printf (IO_BUF, "read complete in %d ticks\n",start_count_read-end_count_read);
#endif
#endif
		//increment segment index
		seg_index++;
		
		#ifdef PROFILE
		  start_count_process = tc[T2_COUNT];
		#endif
		
		//choose current buffers
		if(!read_switch && !write_switch)
		{
#ifdef PRINT 
io_printf (IO_BUF, "buff_b-->buff_x\n");
#endif
			index_x=process_chan(dtcm_buffer_x,dtcm_buffer_b);
		}
		else if(!read_switch && write_switch)
		{
#ifdef PRINT
io_printf (IO_BUF, "buff_b-->buff_y\n"); 
#endif
			index_y=process_chan(dtcm_buffer_y,dtcm_buffer_b);
		}
		else if(read_switch && !write_switch)
		{
#ifdef PRINT
io_printf (IO_BUF, "buff_a-->buff_x\n");
#endif	
			index_x=process_chan(dtcm_buffer_x,dtcm_buffer_a);

		}
		else
		{
#ifdef PRINT
io_printf (IO_BUF, "buff_a-->buff_y\n");
#endif
			index_y=process_chan(dtcm_buffer_y,dtcm_buffer_a);
		}		  
			
		#ifdef PROFILE
			  end_count_process = tc[T2_COUNT];
			  dtcm_profile_buffer[1+((seg_index-1)*3)]=start_count_process-end_count_process;
		#ifdef PRINT 
			io_printf (IO_BUF, "process complete in %d ticks (segment %d)\n",start_count_process-end_count_process,seg_index);
		#endif	
		#endif			
		
		spin1_trigger_user_event(NULL,NULL);
	}
	else if (ttag==DMA_WRITE)
	{
#ifdef PROFILE
  end_count_write = tc[T2_COUNT];
  dtcm_profile_buffer[2+((seg_index-1)*3)]=start_count_write-end_count_write;
#ifdef PRINT 
  io_printf (IO_BUF, "write complete in %d ticks\n",start_count_write-end_count_write);
#endif
#endif
		//flip write buffers
		write_switch=!write_switch;
	}
	else
	{
		#ifdef PRINT
		io_printf(IO_BUF,"[core %d] invalid %d DMA tag!\n",coreID,ttag);
		#endif
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
		profile_buffer[i]  = dtcm_profile_buffer[i];
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

  //set timer tick
  spin1_set_timer_tick (TIMER_TICK_PERIOD);

  //setup callbacks
  //process channel once data input has been read to DTCM
  spin1_callback_on (DMA_TRANSFER_DONE,transfer_handler,0);
  //reads from DMA to DTCM every tick
  spin1_callback_on (TIMER_TICK,data_read,-1);
  spin1_callback_on (USER_EVENT,data_write,0);

  app_init();

  spin1_start (SYNC_WAIT);
  
  app_done ();

}

