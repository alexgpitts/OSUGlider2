#pragma once

#include "../common.h"

// this controls the size of the welch window
// COLS must be power of 2, hence 'POW_OF_2'
#define POW_OF_2 5

#define ROWS 24
#define COLS (1<<POW_OF_2)


F32 Table[ROWS][COLS] = {0};

// Input buffer used for raw acceleration data in x,y,z directions
// ECE team defined their own input buffer for 'process' to access
// thus, when using C implementation on microcontroller, this buffer
// should be commented out. This buffer is mainly for desktop testing
// CSV parser assumes access to this buffer
// Input[0] = x
// Input[1] = y
// Input[2] = z
#define INPUTS 3
#define INPUT_MAX (1<<12)
F32 Input[INPUTS][INPUT_MAX] = {0};