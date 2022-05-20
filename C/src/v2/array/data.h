#pragma once

#include "../common.h"

// this controls the size of the welch window
#define POW_OF_2 5

#define ROWS 24
#define COLS (1<<POW_OF_2)


F32 Table[ROWS][COLS] = {0};


#define INPUTS 3
#define INPUT_MAX (1<<12)
F32 Input[INPUTS][INPUT_MAX] = {0};