/**
  ******************************************************************************
  * @file    
  * @author  Giuseppe Callipo - IK8YFW - ik8yfw@libero.it
  * @version V2.0.0
  * @date    29-01-2022
  * @brief    file
  *
  ******************************************************************************
  *
  *  NOTE:
  *  This is a rework of DSPham parts by  Graham Whaley
  *
  * spectral weighting noise reduction
  * based on:
  * Kim, H.-G. & D. Ruwisch (2002): Speech enhancement in non-stationary noise environments. 
  * – 7th International Conference on Spoken Language Processing [ICSLP 2002]. 
  * – ISCA Archive (http://www.isca-speech.org/archive)
  *
  ***********************************************************************************************************
  *    Noise reduction with spectral subtraction rule
  *    based on Romanin et al. 2009 & Schmitt et al. 2002
  *    and MATLAB voicebox
  *    and Gerkmann & Hendriks 2002
  *    and Yao et al. 2016
  *
  * STAND: UHSDR github 14.1.2018
  ************************************************************************************************************
  *
  * Frank DD4WH & Michael DL2FW, November 2017
  * NOISE REDUCTION BASED ON SPECTRAL SUBTRACTION
  * following Romanin et al. 2009 on the basis of Ephraim & Malah 1984 and Hu et al. 2001
  * detailed technical description of the implemented algorithm
  * can be found in our WIKI
  * https://github.com/df8oe/UHSDR/wiki/Noise-reduction
  *
  * half-overlapping input buffers (= overlap 50%)
  * sqrt von Hann window before FFT
  * sqrt von Hann window after inverse FFT
  * FFT256 - inverse FFT256
  * overlap-add
  * 
  *************************************************************************************************************/
  
 #include "RDSP_nr_spectral.h"
  
/*************************************************************************************************************/
/*                               Action SNR in dB
/*************************************************************************************************************/
float32_t asnr = 20;        // default active SNR in dB - Can be changed as: 30 or 40 dB 

/*************************************************************************************************************/
/*                               Base NR tuned parameters
/*************************************************************************************************************/
#define NR_FFT_L    256 // or FFT_length  // Fix the fft lenght for nr analisys
#define NR_N_frames 15
float32_t tinc = 0.00145;   // frame time 5.3333ms
const int16_t NR_width = 4;
const float32_t power_threshold = 0.4;
const float32_t tax = -tinc / log(0.8);  // noise output smoothing time constant = -tinc/ln(0.8)
const float32_t tap = -tinc / log(0.9);  // speech prob smoothing time constant = -tinc/ln(0.9) tinc = frame time (5.33ms)
const float32_t psthr = 0.99;            // threshold for smoothed speech probability [0.99]
const float32_t pnsaf = 0.01;            // noise probability safety value [0.01]
const float32_t psini = 0.5;             // initial speech probability [0.5]
const float32_t pspri = 0.5;             // prior speech probability [0.5]
const float32_t tavini = 0.064;

/*************************************************************************************************************/
/*                               Vorking variables and more ...
/*************************************************************************************************************/
const float32_t sqrtHann[256] = {0, 0.01231966, 0.024637449, 0.036951499, 0.049259941, 0.061560906, 0.073852527, 0.086132939, 
                                   0.098400278, 0.110652682, 0.122888291, 0.135105247, 0.147301698, 0.159475791, 0.171625679, 
                                   0.183749518, 0.195845467, 0.207911691, 0.219946358, 0.231947641, 0.24391372, 0.255842778,
                                   0.267733003, 0.279582593, 0.291389747, 0.303152674, 0.314869589, 0.326538713, 0.338158275,
                                   0.349726511, 0.361241666, 0.372701992, 0.384105749, 0.395451207, 0.406736643, 0.417960345,
                                   0.429120609, 0.440215741, 0.451244057, 0.462203884, 0.473093557, 0.483911424, 0.494655843,
                                   0.505325184, 0.515917826, 0.526432163, 0.536866598, 0.547219547, 0.557489439, 0.567674716, 
                                   0.577773831, 0.587785252, 0.597707459, 0.607538946, 0.617278221, 0.626923806, 0.636474236, 
                                   0.645928062, 0.65528385, 0.664540179, 0.673695644, 0.682748855, 0.691698439, 0.700543038, 
                                   0.709281308, 0.717911923, 0.726433574, 0.734844967, 0.743144825, 0.75133189, 0.759404917, 
                                   0.767362681, 0.775203976, 0.78292761, 0.790532412, 0.798017227, 0.805380919, 0.812622371,
                                   0.819740483, 0.826734175, 0.833602385, 0.840344072, 0.846958211, 0.853443799, 0.859799851, 
                                   0.866025404, 0.872119511, 0.878081248, 0.88390971, 0.889604013, 0.895163291, 0.900586702, 
                                   0.905873422, 0.911022649, 0.916033601, 0.920905518, 0.92563766, 0.930229309, 0.934679767, 
                                   0.938988361, 0.943154434, 0.947177357, 0.951056516, 0.954791325, 0.958381215, 0.961825643, 
                                   0.965124085, 0.968276041, 0.971281032, 0.974138602, 0.976848318, 0.979409768, 0.981822563, 
                                   0.984086337, 0.986200747, 0.988165472, 0.989980213, 0.991644696, 0.993158666, 0.994521895, 
                                   0.995734176, 0.996795325, 0.99770518, 0.998463604, 0.999070481, 0.99952572, 0.99982925, 
                                   0.999981027, 0.999981027, 0.99982925, 0.99952572, 0.999070481, 0.998463604, 0.99770518, 
                                   0.996795325, 0.995734176, 0.994521895, 0.993158666, 0.991644696, 0.989980213, 0.988165472, 
                                   0.986200747, 0.984086337, 0.981822563, 0.979409768, 0.976848318, 0.974138602, 0.971281032, 
                                   0.968276041, 0.965124085, 0.961825643, 0.958381215, 0.954791325, 0.951056516, 0.947177357, 
                                   0.943154434, 0.938988361, 0.934679767, 0.930229309, 0.92563766, 0.920905518, 0.916033601, 
                                   0.911022649, 0.905873422, 0.900586702, 0.895163291, 0.889604013, 0.88390971, 0.878081248, 
                                   0.872119511, 0.866025404, 0.859799851, 0.853443799, 0.846958211, 0.840344072, 0.833602385, 
                                   0.826734175, 0.819740483, 0.812622371, 0.805380919, 0.798017227, 0.790532412, 0.78292761, 
                                   0.775203976, 0.767362681, 0.759404917, 0.75133189, 0.743144825, 0.734844967, 0.726433574, 
                                   0.717911923, 0.709281308, 0.700543038, 0.691698439, 0.682748855, 0.673695644, 0.664540179, 
                                   0.65528385, 0.645928062, 0.636474236, 0.626923806, 0.617278221, 0.607538946, 0.597707459, 
                                   0.587785252, 0.577773831, 0.567674716, 0.557489439, 0.547219547, 0.536866598, 0.526432163, 
                                   0.515917826, 0.505325184, 0.494655843, 0.483911424, 0.473093557, 0.462203884, 0.451244057, 
                                   0.440215741, 0.429120609, 0.417960345, 0.406736643, 0.395451207, 0.384105749, 0.372701992, 
                                   0.361241666, 0.349726511, 0.338158275, 0.326538713, 0.314869589, 0.303152674, 0.291389747, 
                                   0.279582593, 0.267733003, 0.255842778, 0.24391372, 0.231947641, 0.219946358, 0.207911691, 
                                   0.195845467, 0.183749518, 0.171625679, 0.159475791, 0.147301698, 0.135105247, 0.122888291, 
                                   0.110652682, 0.098400278, 0.086132939, 0.073852527, 0.061560906, 0.049259941, 0.036951499, 
                                   0.024637449, 0.01231966, 0};

// spectral specific vars  
const static arm_cfft_instance_f32 *spec_FFT;
float32_t DMAMEM NR_Hk_old[NR_FFT_L / 2];
float32_t DMAMEM NR_Nest[NR_FFT_L / 2][2]; //
float32_t DMAMEM NR_SNR_post[NR_FFT_L / 2];
float32_t DMAMEM NR_SNR_prio[NR_FFT_L / 2];
float32_t DMAMEM NR_long_tone_gain[NR_FFT_L / 2];

// NR stuff - shared by Kim and spectral at least.
const static arm_cfft_instance_f32 *NR_FFT;
const static arm_cfft_instance_f32 *NR_iFFT;
float32_t DMAMEM NR_last_sample_buffer_L [NR_FFT_L / 2];
float32_t DMAMEM NR_M[NR_FFT_L / 2]; // minimum of the 20 last values of E
//now define const uint8_t NR_N_frames = 15; // default 24 //40 //12 //20 //18//12 //20
float32_t DMAMEM NR_E[NR_FFT_L / 2][NR_N_frames]; // averaged (over the last four values) X values for the last 20 FFT frames
float32_t DMAMEM NR_X[NR_FFT_L / 2][3]; // magnitudes (fabs) of the last four values of FFT results for 128 frequency bins
float32_t DMAMEM NR_G[NR_FFT_L / 2]; // preliminary gain factors (before time smoothing) and after that contains the frequency smoothed gain factors
float32_t DMAMEM NR_FFT_buffer [NR_FFT_L] __attribute__ ((aligned (4)));  //was 512
float32_t NR_alpha = 0.95; // default 0.99 --> range 0.98 - 0.9999; 0.95 acts much too hard: reverb effects
float32_t DMAMEM NR_last_iFFT_result [NR_FFT_L / 2];
float32_t DMAMEM NR_Gts[NR_FFT_L / 2][2]; // time smoothed gain factors (current and last) for each of the 128 bins
uint8_t   NR_first_time = 1;
float lastNRSetting =0;

void spectral_noise_reduction_init()
{
  // Init code imported from main init routine.
  /****************************************************************************************
     init complex FFTs
  ****************************************************************************************/
  //spec_FFT = &arm_cfft_sR_f32_len256;
  NR_FFT = &arm_cfft_sR_f32_len256;
  NR_iFFT = &arm_cfft_sR_f32_len256;
    
  for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++)
  {
    NR_last_sample_buffer_L[bindx] = 0.1;
    NR_Hk_old[bindx] = 0.1; // old gain
    NR_Nest[bindx][0] = 0.01;
    NR_Nest[bindx][1] = 0.015;
    NR_Gts[bindx][1] = 0.1;
    NR_M[bindx] = 500.0;
    NR_E[bindx][0] = 0.1;
    NR_X[bindx][1] = 0.5;
    NR_SNR_post[bindx] = 2.0;
    NR_SNR_prio[bindx] = 1.0;
    NR_first_time = 2;
    NR_long_tone_gain[bindx] = 1.0;
  }

  //These bits taken from the Convolution 'DMA var' init code - they should probably really be in a
  //global 'NR' init section - a bunch is shared between Kim and spectral.
  for (unsigned i = 0; i < NR_FFT_L / 2; i++)
  {
      NR_last_iFFT_result[i] = 0.0;
      NR_last_sample_buffer_L [i] = 0.0;
      NR_M[i] = 0.0;
      NR_G[i] = 0.0;
      NR_SNR_prio[i] = 0.0;
      NR_SNR_post[i] = 0.0;
      NR_Hk_old[i] = 0.0;
  }

  for (unsigned j = 0; j < 3; j++)
  {
      for(unsigned i=0; i<NR_FFT_L / 2; i++)
      {
          NR_X[i][j] = 0.0;
      }
  }

  for (unsigned j = 0; j < 2; j++)
  {
      for(unsigned i=0; i<NR_FFT_L / 2; i++)
      {
          NR_Nest[i][j] = 0.0;
          NR_Gts[i][j] = 0.0;
      }
  }

  for (unsigned j = 0; j < NR_N_frames; j++)
  {
      for(unsigned i=0; i<NR_FFT_L / 2; i++)
      {
          NR_E[i][j] = 0.0;
      }
  }
}

void spectral_noise_reduction ()
{ 
  static uint8_t NR_init_counter = 0;
  uint8_t VAD_low = 0;
  uint8_t VAD_high = 127;
  float32_t lf_freq; // = (offset - width/2) / (12000 / NR_FFT_L); // bin BW is 46.9Hz [12000Hz / 256 bins] @96kHz
  float32_t uf_freq; //= (offset + width/2) / (12000 / NR_FFT_L);
  static float32_t ax; //=0.8; // ax=exp(-tinc/tax); % noise output smoothing factor
  static float32_t ap; //=0.9; // ap=exp(-tinc/tap); % noise output smoothing factor
  static float32_t xih1; // = 31.6;
  ax = expf(-tinc / tax);
  ap = expf(-tinc / tap);
  xih1 = powf(10, (float32_t)asnr / 10.0);
  static float32_t xih1r = 1.0 / (1.0 + xih1) - 1.0;
  static float32_t pfac = (1.0 / pspri - 1.0) * (1.0 + xih1);
  float32_t snr_prio_min = powf(10, - (float32_t)20 / 20.0);
  static float32_t pslp[NR_FFT_L / 2];
  static float32_t xt[NR_FFT_L / 2];
  static float32_t xtr;
  static float32_t pre_power;
  static float32_t post_power;
  static float32_t power_ratio;
  static int16_t NN;
  float32_t ph1y[NR_FFT_L / 2];
  static int NR_first_time_2 = 1;
  
  // Audio Band some lo/hi cutoff freqs - vaguely SSB
  lf_freq = 100;
  uf_freq = 3600; 

  //Our sample rate is fixed...
  lf_freq /= (SAMPLE_RATE / DF) / NR_FFT_L;
  uf_freq /= (SAMPLE_RATE / DF) / NR_FFT_L;

  // INITIALIZATION ONCE 1
  if (NR_first_time_2 == 1)
  { // Initialize all the variables
    for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++)
    {
      NR_last_sample_buffer_L[bindx] = 0.0;
      NR_G[bindx] = 1.0;
      NR_Hk_old[bindx] = 1.0; // old gain or xu in development mode
      NR_Nest[bindx][0] = 0.0;
      NR_Nest[bindx][1] = 1.0;
      pslp[bindx] = 0.5;
    }
    NR_first_time_2 = 2; // we need to do some more a bit later down
  }

  for (int k = 0; k < 2; k++)
  {
    // NR_FFT_buffer is 256 floats big
    // interleaved r, i, r, i . . .
    // fill first half of FFT_buffer with last events audio samples
    for (int i = 0; i < NR_FFT_L / 2; i++)
    {
      NR_FFT_buffer[i * 2] = NR_last_sample_buffer_L[i]; // real
      NR_FFT_buffer[i * 2 + 1] = 0.0; // imaginary
    }
    // copy recent samples to last_sample_buffer for next time!
    for (int i = 0; i < NR_FFT_L  / 2; i++)
    {
      NR_last_sample_buffer_L [i] = float_buffer_L[i + k * (NR_FFT_L / 2)];
    }
    // now fill recent audio samples into second half of FFT_buffer
    for (int i = 0; i < NR_FFT_L / 2; i++)
    {
      NR_FFT_buffer[NR_FFT_L + i * 2] = float_buffer_L[i + k * (NR_FFT_L / 2)]; // real
      NR_FFT_buffer[NR_FFT_L + i * 2 + 1] = 0.0;
    }
    /////////////////////////////////
    // WINDOWING
    // perform windowing on samples in the NR_FFT_buffer
    
    for (int idx = 0; idx < NR_FFT_L; idx++)
    { // sqrt Hann window
      NR_FFT_buffer[idx * 2] *= sqrtHann[idx];
    }

    // NR_FFT
    // calculation is performed in-place the FFT_buffer [re, im, re, im, re, im . . .]
    arm_cfft_f32(NR_FFT, NR_FFT_buffer, 0, 1);

    //##########################################################################################################################################
    //##########################################################################################################################################
    //##########################################################################################################################################

    for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++)
    {
      // this is squared magnitude for the current frame
      NR_X[bindx][0] = (NR_FFT_buffer[bindx * 2] * NR_FFT_buffer[bindx * 2] + NR_FFT_buffer[bindx * 2 + 1] * NR_FFT_buffer[bindx * 2 + 1]);
    }

    if (NR_first_time_2 == 2)
    { // Initialize all the variables
      for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++)
      {
        NR_Nest[bindx][0] = NR_Nest[bindx][0] + 0.05 * NR_X[bindx][0]; //  we do it 20 times to average over 20 frames 
                                                                       // for app. 100ms only on NR_on/bandswitch/modeswitch,...
        xt[bindx] = psini * NR_Nest[bindx][0];
      }
      NR_init_counter++;
      if (NR_init_counter > 19)//average over 20 frames for app. 100ms
      {
        NR_init_counter = 0;
        NR_first_time_2 = 3;  // now we did all the necessary initialization to actually start the noise reduction
      }
    }

    if (NR_first_time_2 == 3)
    {
      // Execute noise estimate MMSE based.
      for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++) // 1. Step of NR - calculate the SNR's
      {
        ph1y[bindx] = 1.0 / (1.0 + pfac * expf(xih1r * NR_X[bindx][0] / xt[bindx]));
        pslp[bindx] = ap * pslp[bindx] + (1.0 - ap) * ph1y[bindx];

        if (pslp[bindx] > psthr)
        {
          ph1y[bindx] = 1.0 - pnsaf;
        }
        else
        {
          ph1y[bindx] = fmin(ph1y[bindx] , 1.0);
        }
        xtr = (1.0 - ph1y[bindx]) * NR_X[bindx][0] + ph1y[bindx] * xt[bindx];
        xt[bindx] = ax * xt[bindx] + (1.0 - ax) * xtr;
      }

      for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++) // 1. Step of NR - calculate the SNR's
      {
        NR_SNR_post[bindx] = fmax(fmin(NR_X[bindx][0] / xt[bindx], 1000.0), snr_prio_min); 
        // limited to +30 /-15 dB, might be still too much of reduction, let's try it?
        NR_SNR_prio[bindx] = fmax(NR_alpha * NR_Hk_old[bindx] + (1.0 - NR_alpha) * fmax(NR_SNR_post[bindx] - 1.0, 0.0), 0.0);
      }

      VAD_low = (int)lf_freq;
      VAD_high = (int)uf_freq;
      if (VAD_low == VAD_high)
      {
        VAD_high++;
      }
      if (VAD_low < 1)
      {
        VAD_low = 1;
      }
      else if (VAD_low > NR_FFT_L / 2 - 2)
      {
        VAD_low = NR_FFT_L / 2 - 2;
      }
      if (VAD_high < 1)
      {
        VAD_high = 1;
      }
      else if (VAD_high > NR_FFT_L / 2)
      {
        VAD_high = NR_FFT_L / 2;
      }
      // (eq. 12 of Schmitt et al. 2002, eq. 9 of Romanin et al. 2009)
      // 4    calculate v = SNRprio(n, bin[i]) / (SNRprio(n, bin[i]) + 1) * SNRpost(n, bin[i]) 
      //      and calculate the HK's

      for (int bindx = VAD_low; bindx < VAD_high; bindx++) // maybe we should limit this to the signal containing bins (filtering!!)
      {
        float32_t v = NR_SNR_prio[bindx] * NR_SNR_post[bindx] / (1.0 + NR_SNR_prio[bindx]);

        NR_G[bindx] = 1.0 / NR_SNR_post[bindx] * sqrtf((0.7212 * v + v * v));

        NR_Hk_old[bindx] = NR_SNR_post[bindx] * NR_G[bindx] * NR_G[bindx]; //
      }

      // MUSICAL NOISE TREATMENT HERE, DL2FW
      // musical noise "artefact" reduction by dynamic averaging - depending on SNR ratio
      pre_power = 0.0;
      post_power = 0.0;
      for (int bindx = VAD_low; bindx < VAD_high; bindx++)
      {
        pre_power += NR_X[bindx][0];
        post_power += NR_G[bindx] * NR_G[bindx]  * NR_X[bindx][0];
      }

      power_ratio = post_power / pre_power;
      if (power_ratio > power_threshold)
      {
        power_ratio = 1.0;
        NN = 1;
      }
      else
      {
        NN = 1 + 2 * (int)(0.5 + NR_width * (1.0 - power_ratio / power_threshold));
      }

      for (int bindx = VAD_low + NN / 2; bindx < VAD_high - NN / 2; bindx++)
      {
        NR_Nest[bindx][0] = 0.0;
        for (int m = bindx - NN / 2; m <= bindx + NN / 2; m++)
        {
          NR_Nest[bindx][0] += NR_G[m];
        }
        NR_Nest[bindx][0] /= (float32_t)NN;
      }

      // and now the edges - only going NN steps forward and taking the average
      // lower edge
      for (int bindx = VAD_low; bindx < VAD_low + NN / 2; bindx++)
      {
        NR_Nest[bindx][0] = 0.0;
        for (int m = bindx; m < (bindx + NN); m++)
        {
          NR_Nest[bindx][0] += NR_G[m];
        }
        NR_Nest[bindx][0] /= (float32_t)NN;
      }

      // upper edge - only going NN steps backward and taking the average
      for (int bindx = VAD_high - NN; bindx < VAD_high; bindx++)
      {
        NR_Nest[bindx][0] = 0.0;
        for (int m = bindx; m > (bindx - NN); m--)
        {
          NR_Nest[bindx][0] += NR_G[m];
        }
        NR_Nest[bindx][0] /= (float32_t)NN;
      }

      // end of edge treatment
      for (int bindx = VAD_low + NN / 2; bindx < VAD_high - NN / 2; bindx++)
      {
        NR_G[bindx] = NR_Nest[bindx][0];
      }
      // end of musical noise reduction

    } //end of "if ts.nr_first_time == 3"


    //##########################################################################################################################################
    //##########################################################################################################################################
    //##########################################################################################################################################

    // FINAL SPECTRAL WEIGHTING: Multiply current FFT results with NR_FFT_buffer for 128 bins with the 128 bin-specific gain factors G
    //              for(int bindx = 0; bindx < NR_FFT_L / 2; bindx++) // try 128:
    for (int bindx = 0; bindx < NR_FFT_L / 2; bindx++) // try 128:
    {
      NR_FFT_buffer[bindx * 2] = NR_FFT_buffer [bindx * 2] * NR_G[bindx] * NR_long_tone_gain[bindx]; // real part
      NR_FFT_buffer[bindx * 2 + 1] = NR_FFT_buffer [bindx * 2 + 1] * NR_G[bindx] * NR_long_tone_gain[bindx]; // imag part
      NR_FFT_buffer[NR_FFT_L * 2 - bindx * 2 - 2] = NR_FFT_buffer[NR_FFT_L * 2 - bindx * 2 - 2] * NR_G[bindx] * NR_long_tone_gain[bindx]; // real part conjugate symmetric
      NR_FFT_buffer[NR_FFT_L * 2 - bindx * 2 - 1] = NR_FFT_buffer[NR_FFT_L * 2 - bindx * 2 - 1] * NR_G[bindx] * NR_long_tone_gain[bindx]; // imag part conjugate symmetric
    }

    /*****************************************************************
       NOISE REDUCTION CODE ENDS HERE
     *****************************************************************/

    // NR_iFFT
    // perform iFFT (in-place)
    arm_cfft_f32(NR_iFFT, NR_FFT_buffer, 1, 1);

    // perform windowing on samples in the NR_FFT_buffer
    for (int idx = 0; idx < NR_FFT_L; idx++)
    { // sqrt Hann window
      NR_FFT_buffer[idx * 2] *= sqrtHann[idx];
    }

    // do the overlap & add
    for (int i = 0; i < NR_FFT_L / 2; i++)
    { // take real part of first half of current iFFT result and add to 2nd half of last iFFT_result
      //              NR_output_audio_buffer[i + k * (NR_FFT_L / 2)] = NR_FFT_buffer[i * 2] + NR_last_iFFT_result[i];
      float_buffer_L[i + k * (NR_FFT_L / 2)] = NR_FFT_buffer[i * 2] + NR_last_iFFT_result[i];
      float_buffer_R[i + k * (NR_FFT_L / 2)] = float_buffer_L[i + k * (NR_FFT_L / 2)];
    }
    for (int i = 0; i < NR_FFT_L / 2; i++)
    {
      NR_last_iFFT_result[i] = NR_FFT_buffer[NR_FFT_L + i * 2];
    }
    // end of "for" loop which repeats the FFT_iFFT_chain two times !!!
  }
} // end of Romanin algorithm
