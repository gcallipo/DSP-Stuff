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

#ifndef RDSP_FIR_H_INCLUDED
#define RDSP_FIR_H_INCLUDED

#include <Audio.h>
#include <arm_math.h>
#include <arm_const_structs.h>

extern void calc_FIR_coeffs (float * coeffs_I, int numCoeffs, float32_t fc, float32_t Astop, int type, float dfc, float Fsamprate);


#endif /* RDSP_FIR_H_INCLUDED */

/**************************************END OF FILE****/
