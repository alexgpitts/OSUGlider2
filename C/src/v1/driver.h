#pragma once

#include "common.h"
// #include "array/meow_fft.h"
#include "array/code.h"
#include "read_csv.h"
#include "processdata.h"
#include "array/FFT.h"

void process() {

	Index timespace = ROW(0);
	// Scale(0.781250024, Iota(timespace), timespace);

	F32 freq = 1.2799999713897705;
	Scale(freq, Iota(timespace), timespace);

	Index t = ROW(0);
	Index x = ROW(1);
	Index y = ROW(2);
	Index z = ROW(3);

	Rolling_mean(2, ROW(1), ROW(5));
	Rolling_mean(2, ROW(2), ROW(6));
	Rolling_mean(2, ROW(3), ROW(7));

	C128 H[COLS] = {0};
	C128 TMP[COLS] = {0};
	for (UZ i = 0; i < COLS; i++) {
		H[i].real = Table[5][i];
	}

	// Mov(ROW(4), ROW(8));
	FFT(H, TMP, COLS);
	// for (UZ i = 0; i < COLS; i++) {
	// 	note: index0's real is at 0, but the imag is at n/2+1 (8)
	// 	Table[22][i] = H[i].real;
	// 	Table[23][i] = H[i].imag;
	// }
	// for (UZ i = 0; i < COLS/2; i++) {
	// 	Table[22][i*2] = H[i].real;
	// 	Table[22][i*2+1] = H[i].imag;
	// }
	for (UZ i = 0; i < COLS/2; i++) {
		((C64*) Table[22])[i] = (C64) {
			.real = (float) H[i].real,
			.imag = (float) H[i].imag,
		};
	}

	Index x_fft = meow_rFFT(ROW(5), ROW(9));
	Index y_fft = meow_rFFT(ROW(6), ROW(10));
	Index z_fft = meow_rFFT(ROW(7), ROW(11));

	
	
	CalcPSD(freq, x_fft, x_fft, ROW(13));
	// CalcPSD(20, y_fft, y_fft, ROW(14));
	// CalcPSD(20, z_fft, z_fft, ROW(15));

	// CalcPSD(20, x_fft, y_fft, ROW(16));
	// CalcPSD(20, x_fft, z_fft, ROW(17));
	// CalcPSD(20, y_fft, z_fft, ROW(18));

}