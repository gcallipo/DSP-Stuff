// SPDX-License-Identifier: GNU General Public License v3.0 or later
#ifndef GLOBAL_NOISE_REDUCTION_H_INCLUDED
#define GLOBAL_NOISE_REDUCTION_H_INCLUDED

#include <Audio.h>
#include <arm_math.h>

//this is the raw hardware rate, before decimation
#define SAMPLE_RATE ((float32_t)AUDIO_SAMPLE_RATE_EXACT)

// Decimate down to 11kHz - see if that helps the NR systems perform, as they will then not be trying
// to operate on the 5-20kHz data, which we never listen to anyway!
#define DF 4
#define FFT_length 256
#define BUFFER_SIZE 128
#define N_BLOCKS (FFT_length / BUFFER_SIZE * DF)

extern float32_t DMAMEM float_buffer_L[];
extern float32_t DMAMEM float_buffer_R[];
 
#endif
