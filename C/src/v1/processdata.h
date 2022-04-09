#pragma once


#include "common.h"
#include "array/data.h"
#include "array/code.h"

// works with a windows sizes that are even
Index Rolling_mean(UZ window_size, Index s_r, Index t_r) {

	if (window_size % 2) {
		printf("Rolling_mean requires an even window size\n");
		exit(1);
	}

	Scan(ADD, s_r, t_r);

	for (Index c = COLS-window_size-1; c >= 1; c --) {
		Table[t_r][c+window_size] =
			Table[t_r][c+window_size] - Table[t_r][c];
	}

	Scale(1.0/((float) window_size), t_r, t_r);
	// center
	Index hl = (window_size-1)/2;
	Index hr = (window_size)/2;
	for (Index c = hr; c < COLS; c++) {
		Table[t_r][c-hr] = Table[t_r][c];
	}
	for (Index c = 0; c < hl; c++) {
		Table[t_r][c] = 0;
	}
	for (Index c = COLS-hr; c < COLS; c++) {
		Table[t_r][c] = 0;
	}
	return t_r;
}





Index CalcPSD(
	float freq,
	Index s_r_x,
	Index s_r_y, 
	Index t_r
) {
	//  "calculates the PSD on an output of a FFT"
   //  nfft = xFFT.size
   //  qOdd = nfft % 2
   //  n = (nfft - qOdd) * 2  # Number of data points input to FFT
   //  w = Bias(n, window)  # Get the window used
   //  wSum = (w * w).sum()		// 1*1 ??
   //  psd = (xFFT.conjugate() * yFFT) / (fs * wSum)
   //  if not qOdd:       # Even number of FFT bins
   //      psd[1:] *= 2   # Real FFT -> double for non-zero freq
   //  else:              # last point unpaired in Nyquist freq
   //      psd[1:-1] *= 2  # Real FFT -> double for non-zero freq
   //  return psd

	UZ n = COLS;
	UZ wSum = COLS;


	// ComplexMul(ComplexConj(s_r_x, t_r), s_r_y, t_r);
	Scale((1. / (freq*wSum)),
		ComplexMul(ComplexConj(s_r_x, t_r), s_r_y, t_r), t_r);

	C64 first = ((C64*) Table[t_r])[0];	// save the first val
	Scale(2, t_r, t_r);
	((C64*) Table[t_r])[0] = first;

	return t_r;
}


