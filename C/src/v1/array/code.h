#pragma once

#include "data.h"
#include "meow_fft.h"


// typedef Index (*OP2)(Index a, Index b, Index z);

typedef enum { ADD, SUB, DIV, MUL } Op;

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

Index Inc(F32 scalar, Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = scalar+Table[s_r][c];
	}
	return t_r;
}

Coord INC(F32 scalar, Coord s, Coord t) {
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




// float Mean(Index s_r) {
// 	return Fold(ADD, s_r) / ((float)COLS);
// }


Index Mov(Index s_r, Index t_r) {
	for (Index c = 0; c < COLS; c++) {
		Table[t_r][c] = Table[s_r][c];
	}

	return t_r;
}


// float SampleStandardDeviation();



Index meow_rFFT(Index s_r, Index t_r) {
	

	// Index Upper = Iota(ROW(0));
	// Index Lower = Inc(1, ROW(1), ROW(1));

	// Scale(1.0/0.0, Mean(Upper,Lower, ROW(2)), ROW(2));





	unsigned          N    = COLS;
	float*            in   = Table[s_r]; //malloc(sizeof(float) * N);
	// Meow_FFT_Complex* out  = malloc(sizeof(Meow_FFT_Complex) * N);
	Meow_FFT_Complex* temp = malloc(sizeof(Meow_FFT_Complex) * N);

	// prepare data for "in" array.
	// ...

	size_t workset_bytes = meow_fft_generate_workset_real(N, NULL);

	// printf("worksetbytes %lu\n", workset_bytes);
	// Get size for a N point fft working on non-complex (real) data.

	Meow_FFT_Workset_Real* fft_real = (Meow_FFT_Workset_Real*) malloc(workset_bytes);

	meow_fft_generate_workset_real(N, fft_real);

	// meow_fft_real(fft_real, in, out);
	meow_fft_real(fft_real, in, (Meow_FFT_Complex*) Table[t_r]);
	// out[0].r == out[0  ].r
	// out[0].j == out[N/2].r

	// for (UZ i = 0; i < COLS/2; i ++) {
	// 	printf("%f %f\n", out[i].r, out[i].j);
	// }

	// printf("%d\n", fft_real->half);

	// meow_fft_real_i(fft_real, in, temp, out);
	// result is not scaled, need to divide all values by N

	free(fft_real);
	// free(out);
	free(temp);
	// free(in);

	return t_r;
}




Index ComplexConj(Index s_r, Index t_r) {

	for (UZ c = 1; c < COLS; c+= 2) {
		Table[t_r][c-1] = Table[s_r][c-1];
		Table[t_r][c] = Table[s_r][c] * (-1);
	}

	return t_r;
}













// Index meow_rFFT(Index s_r, Index t_r) {
	

// 	// Index Upper = Iota(ROW(0));
// 	// Index Lower = Inc(1, ROW(1), ROW(1));

// 	// Scale(1.0/0.0, Mean(Upper,Lower, ROW(2)), ROW(2));





// 	unsigned          N    = COLS;
// 	float*            in   = Table[s_r]; //malloc(sizeof(float) * N);
// 	Meow_FFT_Complex* out  = malloc(sizeof(Meow_FFT_Complex) * N);
// 	Meow_FFT_Complex* temp = malloc(sizeof(Meow_FFT_Complex) * N);

// 	// prepare data for "in" array.
// 	// ...

// 	size_t workset_bytes = meow_fft_generate_workset_real(N, NULL);
// 	// Get size for a N point fft working on non-complex (real) data.

// 	Meow_FFT_Workset_Real* fft_real = (Meow_FFT_Workset_Real*) malloc(workset_bytes);

// 	meow_fft_generate_workset_real(N, fft_real);

// 	meow_fft_real(fft_real, in, out);
// 	// out[0].r == out[0  ].r
// 	// out[0].j == out[N/2].r

// 	for (UZ i = 0; i < COLS/2; i ++) {
// 		printf("%f %f\n", out[i].r, out[i].j);
// 	}

// 	// printf("%d\n", fft_real->half);

// 	// meow_fft_real_i(fft_real, in, temp, out);
// 	// result is not scaled, need to divide all values by N

// 	free(fft_real);
// 	free(out);
// 	free(temp);
// 	// free(in);

// 	return t_r;
// }