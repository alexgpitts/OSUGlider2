#pragma once

#include "common.h"
#include "array/code.h"
#include "read_csv.h"
#include "processdata.h"
#include "array/FFT.h"



void process(
	UZ input_max,
	F32* x_input,
	F32* y_input,
	F32* z_input,
	F32 freq
) {

	// make sure INPUT_MAX is greater than or equal to COLS
	for (UZ i = 0; i <= (input_max-COLS); i += (COLS>>1)) {	// more to read

		Load(x_input, i, ROW(0));
		Load(y_input, i, ROW(1));
		Load(z_input, i, ROW(2));

		Rolling_mean(2, ROW(0), ROW(0));
		Rolling_mean(2, ROW(1), ROW(1));
		Rolling_mean(2, ROW(2), ROW(2));

		Mul(Bias(ROW(3)), ROW(0), ROW(0));
		Mul(Bias(ROW(4)), ROW(1), ROW(1));
		Mul(Bias(ROW(5)), ROW(2), ROW(2));

		Index x_fft = FFT(ROW(0), ROW(0));
		Index y_fft = FFT(ROW(1), ROW(1));
		Index z_fft = FFT(ROW(2), ROW(2));

		CalcPSD(freq, x_fft, x_fft, ROW(3));
		CalcPSD(freq, y_fft, y_fft, ROW(4));
		CalcPSD(freq, z_fft, z_fft, ROW(5));
		
		CalcPSD(freq, x_fft, y_fft, ROW(6));
		CalcPSD(freq, x_fft, z_fft, ROW(7));
		CalcPSD(freq, y_fft, z_fft, ROW(8));

		// // add to accumulator
		Add(ROW(3), ROW(9), ROW(9));
		Add(ROW(4), ROW(10), ROW(10));
		Add(ROW(5), ROW(11), ROW(11));

		Add(ROW(6), ROW(12), ROW(12));
		Add(ROW(7), ROW(13), ROW(13));
		Add(ROW(8), ROW(14), ROW(14));

		// skip first row
		// normalize
		if (i) {
			Scale(0.5, ROW(9), ROW(9));
			Scale(0.5, ROW(10), ROW(10));
			Scale(0.5, ROW(11), ROW(11));

			Scale(0.5, ROW(12), ROW(12));
			Scale(0.5, ROW(13), ROW(13));
			Scale(0.5, ROW(14), ROW(14));
		}
			
	}


	WaveCoefficients(freq);

	return;




}
