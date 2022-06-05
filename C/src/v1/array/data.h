#pragma once
#include "../common.h"

/** Note
 * Unlike v2, v1 does not use the welch method, thus
 * input_max and cols should be the same size
 */


#define POW_OF_2 6

#define ROWS 64
#define COLS (1<<POW_OF_2)


F32 Table[ROWS][COLS] = {0};


#define INPUTS 3
#define INPUT_MAX COLS
F32 Input[INPUTS][INPUT_MAX] = {0};