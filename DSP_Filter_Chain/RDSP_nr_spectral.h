/**
  ******************************************************************************
  * @file    
  * @author  Giuseppe Callipo - IK8YFW - ik8yfw@libero.it
  * @version V2.0.0
  * @date    29-01-2022
  * @brief    file
  *
  ******************************************************************************/
#ifndef SPECTRAL_NOISE_REDUCTION_H_INCLUDED
#define SPECTRAL_NOISE_REDUCTION_H_INCLUDED

#include <Audio.h>
#include <arm_math.h>
#include <arm_const_structs.h>
#include "RDSP_dsp_chain_exchange.h"
 
/************************************************************************/
//                FUNCTIONS 
/************************************************************************/

// Takes input from float_buffer_L.
// Leaves output in both R and L
extern void spectral_noise_reduction_init();
extern void spectral_noise_reduction (void);

#endif
