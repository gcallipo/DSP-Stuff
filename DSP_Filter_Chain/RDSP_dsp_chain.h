/**
  ******************************************************************************
  * @file    
  * @author  Giuseppe Callipo - IK8YFW - ik8yfw@libero.it
  * @version V2.0.0
  * @date    05-01-2022
  * @brief    file
  *
  ******************************************************************************
  *
   */

#ifndef RDSP_DSP_CHAIN_H_INCLUDED
#define RDSP_DSP_CHAIN_H_INCLUDED

#include "RDSP_dsp_chain_exchange.h"
#include "RDSP_fir.h"
#include "RDSP_nr_spectral.h"
#include "RDSP_dynamicFilters.h"

// TODO: ahh here additional 
//#include "RDSP_nr_nlms.h"

//*************************************************************************************************************
//           THE DSP CHAIN USE THE SAMPLE RATE OF 44.100 KHz
//*************************************************************************************************************
#define SAMPLE_RATE ((float32_t)AUDIO_SAMPLE_RATE_EXACT)

//*************************************************************************************************************
//           CONNECTION WITH THE EXTERNAL MAIL PROGRAM AUDIO QUEUE- PLEASE DEFINE THESE VARIBLE OUTSIDE
//*************************************************************************************************************
extern AudioRecordQueue         Q_in_L;
extern AudioRecordQueue         Q_in_R;
extern AudioPlayQueue           Q_out_L;
extern AudioPlayQueue           Q_out_R;

//**************************************************************************************************************
//           DECIMATION / INTERPOLATION BASE SETTINGS
//**************************************************************************************************************

const float32_t n_samplerate = SAMPLE_RATE/1000.0; // samplerate before decimation
const float32_t n_desired_BW = 5.0; // desired max BW of the filters
const float32_t n_att = 90.0;
const float32_t n_fstop = ( (n_samplerate / DF) - n_desired_BW) / n_samplerate;
const float32_t n_fpass = n_desired_BW / n_samplerate;
const uint16_t  n_dec_taps = 1 + (uint16_t) (n_att / (22.0 * (n_fstop - n_fpass)));

// interpolate taps must be divisible by decimation factor - so round up.
const uint16_t  n_int_taps = ((uint16_t)((n_dec_taps + DF) / DF)) * (uint16_t)DF;
const uint16_t  n_int_states = (AUDIO_BLOCK_SAMPLES * N_BLOCKS) + (n_int_taps / DF) -1;

int16_t *inp, *inpR;
int16_t *sp_L, *sp_R;

arm_fir_decimate_instance_f32 FIR_dec;
float32_t FIR_dec_coeffs[(uint16_t)n_dec_taps];
float32_t FIR_dec_state [(int)(n_dec_taps + AUDIO_BLOCK_SAMPLES * N_BLOCKS - 1)];

float32_t FIR_int_coeffs[n_int_taps];
arm_fir_interpolate_instance_f32 FIR_int;
float32_t FIR_int_state [n_int_states];

const uint16_t  n_bpf_taps = 101;
const uint16_t  n_bpf_states = (AUDIO_BLOCK_SAMPLES * N_BLOCKS) + (n_bpf_taps / DF)  -1;

arm_fir_instance_f32 FIR_bpf;
float32_t FIR_bpf_coeffs[n_bpf_taps];
float32_t FIR_bpf_state [n_bpf_states];

//*************************************************************************************************************
//                   EXPORTED FUNCTIONS TO THE MAIN PROGRAM
//*************************************************************************************************************
extern void DSP_CHAIN_applyBFPfilter(double dFLoCut, double dFHiCut);  // Future use
extern void DSP_CHAIN_initialize();
extern void DSP_CHAIN_doProcessing(float iNRLevel, boolean bFilterEnabled, double dFLoCut, double dFHiCut);
//**************************************************************************************************************
//**************************************************************************************************************

//**************************************************************************************************************
//  Implemetation  ... yes must be placed in a cpp file, but now semplify the integration with radiodsp main program
//**************************************************************************************************************

// Call this when need to change the BPF settings
void DSP_CHAIN_applyBFPfilter(double dFLoCut, double dFHiCut){

  AudioNoInterrupts();

 
  /* Call FIR init function to initialize the instance structure. */
  audioFilter(FIR_bpf_coeffs, n_bpf_taps, ID_BANDPASS,  W_HAMMING, dFLoCut, dFHiCut, SAMPLE_RATE/DF);
  arm_fir_init_f32(&FIR_bpf, n_bpf_taps, FIR_bpf_coeffs, FIR_bpf_state, AUDIO_BLOCK_SAMPLES * N_BLOCKS / (uint32_t)DF);


  AudioInterrupts();

}

// Call this to initialize all the dsp chain
void DSP_CHAIN_initialize(){

  // Initialize sub modules
  // Init_LMS_NR(15);   // TODO: Add here the LMS noise reduction
  spectral_noise_reduction_init();
 
  // Initialize the slow path decimator
  calc_FIR_coeffs (FIR_dec_coeffs, n_dec_taps, (float32_t)(n_desired_BW * 1000.0), n_att, 0, 0.0, SAMPLE_RATE);
  if (arm_fir_decimate_init_f32(&FIR_dec, n_dec_taps, (uint32_t)DF , FIR_dec_coeffs, FIR_dec_state, AUDIO_BLOCK_SAMPLES * N_BLOCKS))
  {
    //Serial.print("DEC coeff fail");
    while(1);
  }

  calc_FIR_coeffs (FIR_int_coeffs, n_int_taps, (float32_t)(n_desired_BW * 1000.0), n_att, 0, 0.0, SAMPLE_RATE);
  if (arm_fir_interpolate_init_f32(&FIR_int, (uint8_t)DF, n_int_taps, FIR_int_coeffs, FIR_int_state, AUDIO_BLOCK_SAMPLES * N_BLOCKS / (uint32_t)DF))
  {
    //Serial.print("INT coeff fail");
    while(1);
  }

  // Initialize the filter mask
  
  /* Call FIR init function to initialize the instance structure. */
 audioFilter(FIR_bpf_coeffs, n_bpf_taps, ID_BANDPASS,  W_HAMMING, 200.0, 4500.0, SAMPLE_RATE/DF);
 arm_fir_init_f32(&FIR_bpf, n_bpf_taps, FIR_bpf_coeffs, FIR_bpf_state, AUDIO_BLOCK_SAMPLES * N_BLOCKS / (uint32_t)DF);

  /****************************************************************************************
     begin to queue the audio from the audio library
  ****************************************************************************************/
  Q_in_L.begin();
  Q_in_R.begin();
}


/*- Execute the main processing, at the moment only denoise spectral subtraction is active */
/*  iNRLevel = 0 - disabled
 *  iNRLevel > 0 - enabled
 *  bFilterEnabled - set on/off additional BPF filter (TODO)
 *  Low/High Hz of BPF filter (TODO)
 */
void DSP_CHAIN_doProcessing(float iNRLevel, boolean bFilterEnabled, double dFLoCut, double dFHiCut){

  // Read in N_BLOCKS at a time... we only care about the Left channel, as we are only mono mode right now.
  if (Q_in_L.available() > N_BLOCKS + 0 && Q_in_R.available() > N_BLOCKS + 0)
  {
    //Note when enough data became ready
    for (unsigned i = 0; i < N_BLOCKS; i++)
    {
      // We only process mono audio at the moment, even if the i2s is running in stereo mode..
      inp = Q_in_L.readBuffer();
      inpR = Q_in_R.readBuffer();
      
      arm_q15_to_float (inp, &float_buffer_L[i * AUDIO_BLOCK_SAMPLES], AUDIO_BLOCK_SAMPLES); // convert int_buffer to float 32bit
      arm_q15_to_float (inp, &float_buffer_R[i * AUDIO_BLOCK_SAMPLES], AUDIO_BLOCK_SAMPLES); // convert int_buffer to float 32bit
   
      Q_in_L.freeBuffer();
      Q_in_R.freeBuffer();
    
    }

    //Decimate the data down before we process
    // in-place does not seem to work for us?, so decimate into R, and copy back to L
    arm_fir_decimate_f32(&FIR_dec, float_buffer_L, float_buffer_R, AUDIO_BLOCK_SAMPLES * N_BLOCKS);
    memcpy(float_buffer_L, float_buffer_R, sizeof(float32_t) * BUFFER_SIZE * N_BLOCKS / (uint32_t)(DF));
  
    if (iNRLevel != 0) { // BYPASS
        // Reads input from L, leaves output in R and L
        spectral_noise_reduction();

    } else {   //full bypass mode
      // Just copy over then.
      // Ideally we would not even do the float convert in full bypass mode... but, we don't currently keep
      // the non-float data around for that.
      memcpy(float_buffer_R, float_buffer_L, sizeof(float32_t) * BUFFER_SIZE * N_BLOCKS / (uint32_t)(DF));
    }

    // BPF processing ...
    if (bFilterEnabled){
      arm_fir_f32(&FIR_bpf, float_buffer_R, float_buffer_L, (AUDIO_BLOCK_SAMPLES * N_BLOCKS) / (uint32_t)(DF));
    } else {   //full bypass mode
      memcpy(float_buffer_L, float_buffer_R, sizeof(float32_t) * BUFFER_SIZE * N_BLOCKS / (uint32_t)(DF));
    }
 
    //Interpolate the data back up before we play
    // R->L
    arm_fir_interpolate_f32(&FIR_int, float_buffer_L, float_buffer_R, (AUDIO_BLOCK_SAMPLES * N_BLOCKS) / (uint32_t)(DF));
    //And scale back up after interpolation. Hmm, should we be able to do this scale in the FIR filter itself ?
    //** L -> upscale after decimate -> R **/
    arm_scale_f32(float_buffer_R, DF, float_buffer_L, AUDIO_BLOCK_SAMPLES * N_BLOCKS);

   /**********************************************************************
      CONVERT TO INTEGER AND PLAY AUDIO - Push audio into I2S audio chain
   **********************************************************************/
    for (int i = 0; i < N_BLOCKS; i++)
    {
      sp_L = Q_out_L.getBuffer();    
      sp_R = Q_out_R.getBuffer();
      arm_float_to_q15 (&float_buffer_L[AUDIO_BLOCK_SAMPLES * i], sp_L, AUDIO_BLOCK_SAMPLES); 
      arm_float_to_q15 (&float_buffer_L[AUDIO_BLOCK_SAMPLES * i], sp_R, AUDIO_BLOCK_SAMPLES);
      Q_out_L.playBuffer(); // play it !  
      Q_out_R.playBuffer(); // play it !
    }

  } // end of processing an audio block set

}

//*****************************************
#endif /* RDSP_DSP_CHAIN_H_INCLUDED */

/**************************************END OF FILE****/
