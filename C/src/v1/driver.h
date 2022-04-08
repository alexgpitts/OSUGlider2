#pragma once

#include "common.h"
// #include "array/meow_fft.h"
#include "array/code.h"
#include "read_csv.h"
#include "processdata.h"
#include "array/FFT.h"

void process() {

	Index timespace = ROW(0);
	Scale(0.781250024, Iota(timespace), timespace);

	Index t = ROW(0);
	Index x = ROW(1);
	Index y = ROW(2);
	Index z = ROW(3);

	Rolling_mean(2, ROW(1), ROW(5));
	Rolling_mean(2, ROW(2), ROW(6));
	Rolling_mean(2, ROW(3), ROW(7));


	// test();

	// C128 H[COLS] = {0};
	// C128 TMP[COLS] = {0};
	// for (UZ i = 0; i < COLS; i++) {
	// 	H[i].real = Table[4][i];
	// }

	// // Mov(ROW(4), ROW(8));
	// FFT(H, TMP, COLS);
	// for (UZ i = 0; i < COLS; i++) {
	// 	Table[8][i] = H[i].real;
	// 	Table[9][i] = H[i].imag;
	// }

	Index x_fft = meow_rFFT(ROW(5), ROW(9));
	Index y_fft = meow_rFFT(ROW(6), ROW(10));
	Index z_fft = meow_rFFT(ROW(7), ROW(11));

	
	
	CalcPSD(20, x_fft, x_fft, ROW(13));
	CalcPSD(20, y_fft, y_fft, ROW(14));
	CalcPSD(20, z_fft, z_fft, ROW(15));

	CalcPSD(20, x_fft, y_fft, ROW(16));
	CalcPSD(20, x_fft, z_fft, ROW(17));
	CalcPSD(20, y_fft, z_fft, ROW(18));

	// Index Upper = Iota(ROW(0));
	// Index Lower = Inc(1, ROW(1), ROW(1));

	// Scale(1.0/0.0, Mean(Upper,Lower, ROW(2)), ROW(2));




	// unsigned          N    = 1024;
	// float*            in   = malloc(sizeof(float) * N);
	// Meow_FFT_Complex* out  = malloc(sizeof(Meow_FFT_Complex) * N);
	// Meow_FFT_Complex* temp = malloc(sizeof(Meow_FFT_Complex) * N);

	// // prepare data for "in" array.
	// // ...

	// size_t workset_bytes = meow_fft_generate_workset_real(N, NULL);
	// // Get size for a N point fft working on non-complex (real) data.

	// Meow_FFT_Workset_Real* fft_real =
	// (Meow_FFT_Workset_Real*) malloc(workset_bytes);

	// meow_fft_generate_workset_real(N, fft_real);

	// meow_fft_real(fft_real, in, out);
	// // out[0].r == out[0  ].r
	// // out[0].j == out[N/2].r

	// // meow_fft_real_i(fft_real, in, temp, out);
	// // result is not scaled, need to divide all values by N

	// free(fft_real);
	// free(out);
	// free(temp);
	// free(in);
}