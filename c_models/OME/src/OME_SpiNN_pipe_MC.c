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

//#define PROFILE

//=========GLOBAL VARIABLES============//
double m_pi;
double Fs,dt,max_rate;
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

double conchaL,conchaH,conchaG,earCanalL,earCanalH,earCanalG,stapesH,stapesL,stapesScalar,
	ARtau,ARdelay,ARrateThreshold,rateToAttenuationFactor,BFlength;

double concha_q,concha_j,concha_k,concha_l;

double conchaGainScalar,recip_conchaFilter_a0,
	earCanal_q,earCanal_j,earCanal_k,earCanal_l,earCanalGainScalar,recip_earCanalFilter_a0,
	ARatt,Wn,stapesHP_order,sf,stapesLP_b,stapes_tau,past_stapesDisp;
/*REAL concha_q,concha_j,concha_k,concha_l,conchaGainScalar,recip_conchaFilter_a0,
	earCanal_q,earCanal_j,earCanal_k,earCanal_l,earCanalGainScalar,recip_earCanalFilter_a0,
	ARatt,Wn,stapesHP_order,sf,stapesLP_b,stapes_tau,past_stapesDisp;*/

double conchaFilter_b[3],conchaFilter_a[3],earCanalFilter_b[3],earCanalFilter_a[3],stapesHP_b[3],stapesHP_a[3],stapesLP_a[2],past_input[2],past_concha[2],past_earCanalInput[2],past_earCanal[2],past_stapesInput[2]
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
uint *dtcm_buffer_x;
REAL *dtcm_buffer_y;
REAL *dtcm_profile_buffer;

REAL *sdramin_buffer;
REAL *sdramout_buffer;
REAL *profile_buffer;

// The parameters to be read from memory
enum params_enum {
    DATA_SIZE = 0,
    COREID,
    NUM_DRNL,
    NUM_MACK,
    FS,
    NUM_BFS,
    KEY,
    R2S_KEY,
    RT,
    SHB1,
    SHB2,
    SHB3,
    SHA1,
    SHA2,
    SHA3,
    DATA
    };

struct parameters {
    uint DATA_SIZE;
    uint COREID;
    uint NUM_DRNL;
    uint NUM_MACK;
    uint FS;
    uint NUM_BFS;
    uint KEY;
    uint R2S_KEY;
    uint RT;
    REAL SHB1;
    REAL SHB2;
    REAL SHB3;
    REAL SHA1;
    REAL SHA2;
    REAL SHA3;
    uint DATA;
    };

// The general parameters
struct parameters params;

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

uint rt;

uint num_drnls;

uint num_macks;

uint sampling_frequency;

//application initialisation
void app_init(void)
{
    //m_pi = 4.0*atan(1.0);//acos(-1.0);
    m_pi = PI;//M_PI;
	seg_index=0;
	cbuff_index=0;
	read_switch=0;
	write_switch=0;
	cbuff_numseg=3;

	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);
	//obtain data spec
	address_t data_address = data_specification_get_data_address();
    address_t parameters_address = data_specification_get_region(0, data_address);
    address_t params_enum = data_specification_get_region(0, data_address);

    struct parameters *ome_params = (struct parameters *) parameters_address;
        spin1_memcpy(&params, ome_params, sizeof(struct parameters));

    // Get the size of the data in words
    data_size = params.DATA_SIZE;//params_enum[DATA_SIZE];//
    read_ticks=0;
    TOTAL_TICKS= data_size/SEGSIZE;
    log_info("data_size=%d",data_size);
    log_info("TOTAL_TICKS=%d",TOTAL_TICKS);

    // Get a pointer to the data - not worth copying at present
    //sdramin_buffer = (REAL *) &(params_enum[DATA]);//(REAL *) &params.DATA;//
    sdramin_buffer = (REAL *) &(params_enum[DATA+7]);//to account for double types filter params
    log_info("sdramin_buffer=0x%08x",(uint)sdramin_buffer);
    //obtain this core ID from the host placement perspective
    placement_coreID = params.COREID;//params_enum[COREID];//

    // Get the key to send the data with
    key = params.KEY;//params_enum[KEY];//
    io_printf (IO_BUF, "OME-->DRNL key=%d\n",key);

    r2s_key=params.R2S_KEY;//params[R2S_KEY];//
    log_info("r2s key=%d",r2s_key);

    rt=params.RT;
    log_info("rt=%d",rt);

    //Get number of child DRNL vertices
    num_drnls=params.NUM_DRNL;//params[NUM_DRNL];//
    io_printf (IO_BUF, "num drnls=%d\n",num_drnls);

    //get number of macks
    num_macks=params.NUM_MACK;//params[NUM_MACK];//
    io_printf (IO_BUF, "num macks=%d\n",num_macks);

    //Get sampling frequency
    sampling_frequency = params.FS;//params[FS];//

    Fs= (double)sampling_frequency;
    dt=(1.0/Fs);
    //real-time timer period calculation (us)
    //if(Fs<=22050.0f)TIMER_TICK_PERIOD = (uint)(0.95*1e6 * ((REAL)SEGSIZE/Fs));//RT
    //else TIMER_TICK_PERIOD = 5000;//Not RT
    //TIMER_TICK_PERIOD = (uint)(0.95*1e6 * ((REAL)SEGSIZE/Fs));//RT
    if(rt)TIMER_TICK_PERIOD = (uint)(1e6 * ((REAL)SEGSIZE/Fs));//RT
    else TIMER_TICK_PERIOD = 30000;//6000;//Not RT
     //TIMER_TICK_PERIOD = (uint)(0.8*1e6 * ((REAL)SEGSIZE/Fs));//RT
    //TODO: need to reconsider the running speed needed for total RT performance as the MC TX overhead is large
    //TIMER_TICK_PERIOD = (uint)(1e6 * ((REAL)SEGSIZE/48000.0));//top speed

    log_info("timer period=%d",(uint)TIMER_TICK_PERIOD);
    log_info("Fs=%d",sampling_frequency);
    final_ack=0;

	// Allocate buffers somewhere in SDRAM
    //hack for smaller SDRAM intermediate buffers
	//data_size=cbuff_numseg*SEGSIZE;
	dtcm_buffer_a = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_b = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_x = (uint *) sark_alloc (SEGSIZE, sizeof(uint));
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
		for (uint i = 0; i < SEGSIZE; i++)
		{
			dtcm_buffer_x[i]   = 0;
			//dtcm_buffer_y[i]   = 0;
		}
/*
		
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
    BFlength = params.NUM_BFS;//params[NUM_BFS];//
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

	concha_q= (double)m_pi * (double)dt * (double)(conchaH - conchaL);
	concha_j= 1.0/(1.0+ (1.0/tan(concha_q)));
	concha_k= (2.0 * cos((double)m_pi * (double)dt * (double)(conchaH + conchaL)))
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

	earCanal_q= m_pi * dt * (earCanalH - earCanalL);
	earCanal_j= 1.0/(1.0+ (1.0/tan(earCanal_q)));
	earCanal_k= (2.0 * cos(m_pi * dt * (earCanalH + earCanalL))) / ((1.0 + tan(earCanal_q)) * cos(earCanal_q));
	earCanal_l= (tan(earCanal_q) - 1.0)/(tan(earCanal_q) + 1.0);
	earCanalGainScalar=pow(10.0,earCanalG/20.0);

	earCanalFilter_b[0]=earCanal_j;
	earCanalFilter_b[1]=0.0;
	earCanalFilter_b[2]=-earCanal_j;
	earCanalFilter_a[0]=1.0;
	earCanalFilter_a[1]=-earCanal_k;
	earCanalFilter_a[2]=-earCanal_l;
	recip_earCanalFilter_a0=1.0/earCanalFilter_a[0];

	//stapes filter coeffs pre generated due to butterworth calc code overflow
	stapesHP_b[0] = 0.96937834471788530965596919486415572464466094970703125;//params.SHB1;//0;//0.868412544196881;//0.931905;
	stapesHP_b[1] = -1.9387566894357706193119383897283114492893218994140625;//params.SHB2;//0;//-1.736825088393761;//-1.863809;
	stapesHP_b[2] = 0.96937834471788530965596919486415572464466094970703125;//params.SHB3;//0;//0.868412544196881;//0.931905;
	stapesHP_a[0] = 1.0;//params.SHA1;//0;//1.0;
	stapesHP_a[1] = -1.937818783746783513066702653304673731327056884765625;//params.SHA2;//0;//-1.719434219286951;//-1.859167;
	stapesHP_a[2] = 0.93969459512475761453487166363629512488842010498046875;//params.SHA3;//0;//0.754215957500572;//0.868451;

    log_info("shp_b0:%k",(accum)stapesHP_b[0]);
    log_info("shp_b1:%k",(accum)stapesHP_b[1]);
    log_info("shp_b2:%k",(accum)stapesHP_b[2]);
	log_info("shp_a0:%k",(accum)stapesHP_a[0]);
	log_info("shp_a1:%k",(accum)stapesHP_a[1]);
    log_info("shp_a2:%k",(accum)stapesHP_a[2]);

	stapes_tau= 1.0/ (2 * (double)m_pi * stapesL);
	stapesLP_a[0]= 1.0;
	stapesLP_a[1]= -0.99937168146928201384326939660240896046161651611328125;//dt/stapes_tau -1.0;
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
	REAL *dtcm_buffer_in;
	if(test_DMA == TRUE && sync_count<num_macks && seg_index==0)
	{
	    //log_info("sending r2s");
        //send ready to send MC packet
	    //while (!spin1_send_mc_packet(r2s_key, 0, NO_PAYLOAD)) {
	    while (!spin1_send_mc_packet(r2s_key|1, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
        //spin1_delay_us(10000);//to prevent unnecessarily frequent TX of r2s

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
	else if (read_ticks >= (TOTAL_TICKS+3) && sync_count==num_macks)//+ 3 ticks to account for OME+DRNL+IHCAN processing latency
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
    else if (read_ticks>=TOTAL_TICKS && sync_count==num_macks) read_ticks++;//additional latency wait

}

void process_chan(REAL *in_buffer)
{
    #ifdef PROFILE
        //log_info("enter");
        profiler_write_entry_disable_irq_fiq(PROFILER_ENTER | PROFILER_TIMER);
    #endif
	uint i,j,k;
	uint si=0;
	double concha,earCanalInput,earCanalRes,earCanalOutput,ARoutput,stapesVelocity,stapesDisplacement;
	double filter_1,diff;
	double sub;
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

        /*filter_1 = conchaFilter_b[0]*in_buffer[i]
                    + conchaFilter_b[1]*past_input[0]
                    + conchaFilter_b[2]*past_input[1];
        concha = conchaFilter_a[0]*filter_1
                    - conchaFilter_a[1]*past_concha[0]
                    - conchaFilter_a[2]*past_concha[1];*/
        //filter_3 = conchaFilter_b[2]*past_input[1]+filter_2;
        //filter_4 = conchaFilter_a[1]*past_concha[0];
        //filter_5 = filter_4 - conchaFilter_a[2]*past_concha[1];
        //concha =   filter_5;

		//update vars
		past_input[1]=past_input[0];
		past_input[0]=in_buffer[i];

		past_concha[1]=past_concha[0];
		past_concha[0]=concha;

		earCanalInput= conchaGainScalar * concha + in_buffer[i];

		//ear canal
		/*earCanalRes= (earCanalFilter_b[0]*earCanalInput
				+ earCanalFilter_b[1]*past_earCanalInput[0]
				+ earCanalFilter_b[2]*past_earCanalInput[1]
				- earCanalFilter_a[1]*past_earCanal[0]
				- earCanalFilter_a[2]*past_earCanal[1]) * recip_earCanalFilter_a0;*/
		filter_1 = earCanalFilter_b[0]*earCanalInput
				    + earCanalFilter_b[1]*past_earCanalInput[0]
				    + earCanalFilter_b[2]*past_earCanalInput[1];
        earCanalRes = earCanalFilter_a[0] * filter_1
                      - earCanalFilter_a[1]*past_earCanal[0]
                      - earCanalFilter_a[2]*past_earCanal[1];
		//update vars
		past_earCanalInput[1]=past_earCanalInput[0];
		past_earCanalInput[0]=earCanalInput;
	
		past_earCanal[1]=past_earCanal[0];
		past_earCanal[0]=earCanalRes;

		earCanalOutput= earCanalGainScalar * earCanalRes + earCanalInput;

		//AR
		ARoutput= ARatt * stapesScalar * earCanalOutput;
		//stapes velocity
		/*stapesVelocity= stapesHP_b[0] * ARoutput +
				stapesHP_b[1] * past_stapesInput[0] + 
				stapesHP_b[2] * past_stapesInput[1] 
				- stapesHP_a[1] * past_stapes[0]
				- stapesHP_a[2] * past_stapes[1];*/
        filter_1 = stapesHP_b[0] * ARoutput +
				    stapesHP_b[1] * past_stapesInput[0] +
				    stapesHP_b[2] * past_stapesInput[1];
        stapesVelocity = stapesHP_a[0] * filter_1
                         - stapesHP_a[1] * past_stapes[0]
				         - stapesHP_a[2] * past_stapes[1];
		//update vars
		past_stapesInput[1]= past_stapesInput[0];
		past_stapesInput[0]= ARoutput;
		past_stapes[1]= past_stapes[0];
		past_stapes[0]= stapesVelocity;

		//stapes displacement
		filter_1 = stapesLP_b * stapesVelocity;
		//stapesDisplacement = stapesLP_a[0] * filter_1 - stapesLP_a[1] * past_stapesDisp;
		//stapesDisplacement= stapesLP_b * stapesVelocity - stapesLP_a[1] * past_stapesDisp;
		sub = stapesLP_a[1] * past_stapesDisp;
        //diff = past_stapesDisp - -5.2230809704382550137617236885060008499323186480049752145049524187925271689891815185546875000000000000e-17;
        //diff = past_stapesDisp - -5.1409406724764351467909566793607470131088916641837888166222114705306012183427810668945312500000000000e-17;
        //diff = in_buffer[i]- -1.1120786484703880892813372476771865970945896151533816009759902954101562500000000000000000000000000000e-08;
		//filter_1 = stapesLP_b * -1.3587128514911903004647232071887055606064433500163035617447349068243056535720825195312500000000000000e-15;
		//sub = stapesLP_a[1] *  -5.2230809704382550137617236885060008499323186480049752145049524187925271689891815185546875000000000000e-17;
		stapesDisplacement = filter_1 - sub;
		//update vars
		past_stapesDisp=stapesDisplacement;


		MC_union.f = stapesDisplacement;// - 6.2570166283915330472350035599279830890187798864230522832841074887255672365427017211914062500000000000e-18;
		//MC_union.f = in_buffer[i];//stapesVelocity;//diff;//stapesLP_a[1];//ARoutput;//earCanalOutput;//earCanalInput;//concha;//earCanalRes;//recip_conchaFilter_a0;//

       // log_info("payload = %d",MC_union.u);
		//save to buffer
		//out_buffer[i]=stapesDisplacement;
		//out_buffer[i]=stapesVelocity;//nlin_b0;//nonlinout1a;//nonlinout2b;//linout1;//nonlinout1a;
		//dtcm_buffer_x[i] = MC_union.u;
       /* while (!spin1_send_mc_packet(key,MC_union.u, WITH_PAYLOAD)) {
            spin1_delay_us(1);
        }*/
        spin1_send_mc_packet(key,MC_union.u, WITH_PAYLOAD);
	}
	/*for (i=0;i<SEGSIZE;i++)
	{
        while (!spin1_send_mc_packet(key,dtcm_buffer_x[i], WITH_PAYLOAD)) {
                spin1_delay_us(1);
            }
    }*/

    #ifdef PROFILE
    //log_info("enter");
    profiler_write_entry_disable_irq_fiq(PROFILER_EXIT | PROFILER_TIMER);
    #endif
    //log_info("processing complete %d",seg_index);
}

void transfer_handler(uint tid, uint ttag)
{
	if (ttag==DMA_READ)
	{
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
    log_info("ack mcpacket recieved from key=%d sync_count=%d seg_index=%d\n",mc_key-2,sync_count,seg_index);
}

void app_done ()
{
#ifdef PROFILE
  //profiler_write_entry_disable_irq_fiq(PROFILER_EXIT | PROFILER_TIMER);
  profiler_write_entry_disable_irq_fiq(PROFILER_EXIT | PROFILER_DMA_READ);
  profiler_finalise();
#endif
	log_info("b0:%k",(accum)conchaFilter_b[0]);
    log_info("b1:%k",(accum)conchaFilter_b[1]);
    log_info("b2:%k",(accum)conchaFilter_b[2]);
	log_info("a0:%k",(accum)conchaFilter_a[0]);
	log_info("a1:%k",(accum)conchaFilter_a[1]);
    log_info("a2:%k",(accum)conchaFilter_a[2]);
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
  //io_printf (IO_BUF, "[core %d] -------------------\n", coreID);
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

