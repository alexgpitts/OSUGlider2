#pragma once

#include "common.h"
// #include "array/meow_fft.h"
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

	// timespace is never actually used
	Index timespace = ROW(0);
	Scale(freq, Iota(timespace), timespace);

	Load(x_input, 0, ROW(1));
	Load(y_input, 0, ROW(2));
	Load(z_input, 0, ROW(3));

	Index t = ROW(0);
	Index x = ROW(1);
	Index y = ROW(2);
	Index z = ROW(3);

	Rolling_mean(2, ROW(1), ROW(5));
	Rolling_mean(2, ROW(2), ROW(6));
	Rolling_mean(2, ROW(3), ROW(7));


	Index x_fft = FFT(ROW(5), ROW(9));
	Index y_fft = FFT(ROW(6), ROW(10));
	Index z_fft = FFT(ROW(7), ROW(11));

	
	
	CalcPSD(freq, x_fft, x_fft, ROW(13));
	CalcPSD(freq, y_fft, y_fft, ROW(14));
	CalcPSD(freq, z_fft, z_fft, ROW(15));
	
	CalcPSD(freq, x_fft, y_fft, ROW(16));
	CalcPSD(freq, x_fft, z_fft, ROW(17));
	CalcPSD(freq, y_fft, z_fft, ROW(18));


	WaveCoeffients(freq);
}
