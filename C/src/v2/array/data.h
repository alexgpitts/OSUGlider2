#pragma once

#include "../common.h"

#define POW_OF_2 11

#define ROWS 64
#define COLS (1<<POW_OF_2)
// #define COLS (10)


// F32 Table[ROWS][COLS] = {0};
F32 Table[ROWS][COLS] = {0};

// C64 xTable[ROWS][COLS] = {0};

#define INPUTS 3
#define INPUT_MAX (1<<12)
F32 Input[INPUTS][INPUT_MAX] = {0};