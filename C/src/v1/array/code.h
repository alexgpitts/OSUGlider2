#pragma once

#include "data.h"
#include "FFT.h"



/** array/code.h assumes access to global memory table 'Table' which is
 * defined in array/data.h
 * 
 * the functionality of the operators and the general workflow they imply
 * is explained in /C/readme.md
 */




typedef enum { ADD, SUB, DIV, MUL } Op;

Index Mov(Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = Table[s_r][c];
	}

	return t_r;
}

Index Iota(Index r) {
	for (Index c = 0; c < COLS; c++) {
		Table[r][c] = (F32) c;
	}
	return r;
}


Index Scale(F32 scalar, Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = scalar*Table[s_r][c];
	}
	return t_r;
}

Coord ScaleCell(F32 scalar, Coord s, Coord t) {
	Table[t.row][t.col] = scalar*Table[s.row][s.col];
	return t;
}

Index Inc(F32 scalar, Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = scalar+Table[s_r][c];
	}
	return t_r;
}

Coord IncCell(F32 scalar, Coord s, Coord t) {
	Table[t.row][t.col] = scalar+Table[s.row][s.col];
	return t;
}



Index Reciprocal(Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = 1.0/Table[s_r][c];
	}

	return t_r;
}




Index Add(Index s_r_a, Index s_r_b, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = Table[s_r_a][c]+Table[s_r_b][c];
	}
	return t_r;
}

Index Sub(Index s_r_a, Index s_r_b, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = Table[s_r_a][c]-Table[s_r_b][c];
	}
	return t_r;
}

Index Mul(Index s_r_a, Index s_r_b, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = Table[s_r_a][c]*Table[s_r_b][c];
	}
	return t_r;
}

Index Div(Index s_r_a, Index s_r_b, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = Table[s_r_a][c]/Table[s_r_b][c];
	}
	return t_r;
}

Index ComplexMul(Index s_r_a, Index s_r_b, Index t_r) {
	// C128 C = {0};

	// C.real = (A.real * B.real) - (A.imag * B.imag);
	// C.imag = (A.real * B.imag) + (A.imag * B.real);

	// return C;

	C64 A = {0};
	C64 B = {0};
	C64 C = {0};
	for (Index c = 0; c < COLS/2; c++) {
		A = (((C64*) Table[s_r_a])[c]);
		B = (((C64*) Table[s_r_b])[c]);
		C.real = (A.real * B.real) - (A.imag * B.imag);
		C.imag = (A.real * B.imag) + (A.imag * B.real);

		((C64*) Table[t_r])[c] = C;
	}
	return t_r;
}

Index ComplexScale(C64 scalar, Index s_r, Index t_r) {
	
	C64 A = {0};
	C64 C = {0};
	for (Index c = 0; c < COLS/2; c++) {
		A = (((C64*) Table[s_r])[c]);
		
		C.real = (A.real * scalar.real) - (A.imag * scalar.imag);
		C.imag = (A.real * scalar.imag) + (A.imag * scalar.real);

		((C64*) Table[t_r])[c] = C;
	}
	return t_r;
}

Index SetImag(float imag, Index s_r, Index t_r) {
	for (Index c = 0; c < COLS/2; c++) {
		((C64*) Table[t_r])[c] = (C64) {
			.imag = 0,
			.real = ((C64*) Table[s_r])[c].real
		};
	}

	return t_r;
}

Index SetReal(float real, Index s_r, Index t_r) {
	for (Index c = 0; c < COLS/2; c++) {
		((C64*) Table[s_r])[c].real = real;
	}

	return t_r;
}


Index Shift(I32 shift, Index s_r, Index t_r) {
	if (shift > 0)	{
		for (I32 c = COLS-shift-1; c >= 0; c--) {
			Table[t_r][c+shift] = Table[s_r][c];
		}
		for (UZ c = 0; c < shift; c++) {
			Table[t_r][c] = 0;
		}
	}
	else if (shift < 0) {
		for (I32 c = 0; c < COLS+shift; c++) {
			Table[t_r][c] = Table[s_r][c-shift];
		}
		for (I32 c = shift; c < 0; c++) {
			Table[t_r][COLS+c] = 0;
		}
	}
	else Mov(s_r, t_r);

	return t_r;
}


Index Sqrt(Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = sqrtf(Table[s_r][c]);
	}
	return t_r;
}

Coord SqrtCell(Coord s, Coord t) {
	Table[t.row][t.col] = sqrtf(Table[s.row][s.col]);
	// printf("%g %g\n", Table[t.row][t.col], sqrtf(Table[s.row][s.col]));
		
	return t;
}


float GetCell(Coord s) {
	return Table[s.row][s.col];
}

Coord DivCell(Coord s_a, Coord s_b, Coord t) {
	Table[t.row][t.col] = Table[s_a.row][s_a.col]/Table[s_b.row][s_b.col];
	return t;
}

Coord MovCell(Coord s, Coord t) {
	Table[t.row][t.col] = Table[s.row][s.col];
	return t;
}

Coord SetCell(float val, Coord t) {
	Table[t.row][t.col] = val;
	return t;
}


Coord ArcTanCell(Coord s_a, Coord s_b, Coord t){
	Table[t.row][t.col] = atan2f(
		Table[s_a.row][s_a.col],
		Table[s_b.row][s_b.col]
	);
	// printf("%g %g %g\n",
	// 	Table[s_a.row][s_a.col],
	// 	Table[s_b.row][s_b.col], atan2f(
	// 	Table[s_a.row][s_a.col],
	// 	Table[s_b.row][s_b.col]
	// ));
	return t;
}


Coord MaxCoord (Index s_r) {
	Index max = 0;

	for (UZ c = 1; c < COLS; c++) {
		max = Table[s_r][c] > Table[s_r][max] ? c : max;
	}

	return POS(max, s_r);
}

Coord MaxCell (Index s_r, Coord t) {
	return MovCell(
		MaxCoord(s_r),
		t
	);
}

Coord MaxCoordReal(Index s_r) {
	Index max = 0;

	for (UZ c = 2; c < COLS; c+=2) {
		max = Table[s_r][c] > Table[s_r][max] ? c : max;
	}
	return POS(max, s_r);
}

Coord RadToDegreeCell(Coord s, Coord t){
	Table[t.row][t.col] = 180*Table[s.row][s.col]/M_PI;
	return t;
}

Coord ModCell_immediate(Coord s, float mod, Coord t){
	Table[t.row][t.col] = fmodf(Table[s.row][s.col], mod);
	return t;
}

Index Scan(Op op, Index s_r, Index t_r) {
	Table[t_r][0] = Table[s_r][0];

	switch (op)	{
	case ADD:
		for (Index c = 1; c < COLS; c++) {
			Table[t_r][c] = Table[t_r][c-1]+Table[s_r][c];
		}
		break;

	case SUB:
		for (Index c = 1; c < COLS; c++) {
			Table[t_r][c] = Table[t_r][c-1]-Table[s_r][c];
		}
		break;

	case MUL:
		for (Index c = 1; c < COLS; c++) {
			Table[t_r][c] = Table[t_r][c-1]*Table[s_r][c];
		}
		break;

	case DIV:
		for (Index c = 1; c < COLS; c++) {
			Table[t_r][c] = Table[t_r][c-1]/Table[s_r][c];
		}
		break;
	
	default:
		printf("Error/scan, Invalid operation: %d\n", op);
		exit(1);
		break;
	}

	return t_r;
}


Coord FOLD(Op op, Index s_r, Coord t) {
	// float acc = Table[s_r][0];
	Table[t.row][t.col] = Table[s_r][0];
	switch (op)	{
	case ADD:
		for (Index c = 1; c < COLS; c++) {
			Table[t.row][t.col] += Table[s_r][c];
		}
		break;

	case SUB:
		for (Index c = 1; c < COLS; c++) {
			Table[t.row][t.col] -= Table[s_r][c];
		}
		break;

	case MUL:
		for (Index c = 1; c < COLS; c++) {
			Table[t.row][t.col] *= Table[s_r][c];
		}
		break;

	case DIV:
		for (Index c = 1; c < COLS; c++) {
			Table[t.row][t.col] /= Table[s_r][c];
		}
		break;


	
	default:
		printf("Error/scan, Invalid operation: %d\n", op);
		exit(1);
		break;
	}

	return t;
}





Index ComplexConj(Index s_r, Index t_r) {

	for (UZ c = 1; c < COLS; c+= 2) {
		Table[t_r][c-1] = Table[s_r][c-1];
		Table[t_r][c] = Table[s_r][c] * (-1);
	}

	return t_r;
}


Index FreqSpace(float freq, Index t_r) {
	for (UZ c = 0; c < COLS; c++) {
		Table[t_r][c] = c*(freq/COLS);
	}
	return t_r;
}





Index FFT(Index s_r, Index t_r) {
	
	static C64 TMP[COLS] = {0};
	static C64 data[COLS] = {0};

	memset(TMP, 0, x64_size * COLS);
	memset(data, 0, x64_size * COLS);

	for (UZ i = 0; i < COLS; i++) {
		data[i].real = Table[s_r][i];
	}
	lil_FFT(data, TMP, COLS);				// future: may be replaced with OEM driver
	for (UZ i = 0; i < COLS/2; i++) {
		((C64*) Table[t_r])[i] = data[i];
	}

	return t_r;
}



Index Load(F32* input, Index offset, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = input[offset+c];
	}
	return t_r;
}


Index Cos(Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = cosf(Table[s_r][c]);
	}

	return t_r;
}

