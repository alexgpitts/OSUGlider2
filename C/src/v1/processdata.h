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

	// for (Index c = COLS-window_size-1; c >= 1; c --) {
	// 	Table[t_r][c+window_size] =
	// 		Table[t_r][c+window_size] - Table[t_r][c];
	// }
	// crucial that c is signed here
	for (I32 c = COLS-window_size-1; c >= 0; c --) {
		if (c == 0) {
			printf("%g %g\n", Table[t_r][c+window_size], Table[t_r][c]);
		}
		Table[t_r][c+window_size] =
			Table[t_r][c+window_size] - Table[t_r][c];
	}

	Scale(1.0/((float) window_size), t_r, t_r);
	// center
	Index hl = (window_size)/2;
	Index hr = (window_size-1)/2;
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



/*
    a0 = PSD["zz"][1:] / \
        np.square(np.square(2 * np.pi * PSD["freq_space"][1:]))
    m0 = (a0 * PSD["freq_space"][1]).sum()
    m1 = (a0*PSD["freq_space"][1:]*PSD["freq_space"][1]).sum()
    mm1 = (a0/PSD["freq_space"][1:]*PSD["freq_space"][1]).sum()
    te = mm1/m0  # mean energy period
    m2 = (a0*np.square(PSD["freq_space"][1:]) * PSD["freq_space"][1]).sum()
    tp = 1/PSD["freq_space"][1:][a0.argmax()]

    # directional
    denom = np.sqrt(PSD["zz"] * (PSD["xx"] + PSD["yy"]))
    a1 = PSD["zx"].imag / denom
    b1 = -PSD["zy"].imag / denom
    denom = PSD["xx"] + PSD["yy"]
    dp = np.arctan2(b1[a0.argmax()], a1[a0.argmax()])  # radians

    output = {
        "Hs": 4 * np.sqrt(m0),
        "Ta": m0/m1,  # average period
        "Tp": tp,  # peak wave period
        "wave_energy_ratio": te/tp,
        "Tz": np.sqrt(m0/m2),
        "PeakPSD": a0.max(),
        "Te": te,  # mean energy period
        "Dp": np.arctan2(b1[a0.argmax()], a1[a0.argmax()]),
        "Dp_mag": np.degrees(dp+Data["Meta"]["declination"]) % 360,
        "Dp_true": np.degrees(dp) % 360,
        "A1": a1,
        "B1": b1,
        "A2": (PSD["xx"] - PSD["yy"]) / denom,
        "B2": -2 * PSD["xy"].real / denom,
    }
*/
Index WaveCoeffients(float freq) {
	// FreqSpace populates a row with F32 increments, but I need a
	// freq space for C64, so by halving the frequency, things line up properly
	Index freq_space = FreqSpace(freq*0.5, ROW(20));
	// Inc(freq/COLS, freq_space, freq_space);

	Index psd_xx = ROW(13);
	Index psd_yy = ROW(14);
	Index psd_zz = ROW(15);

	Index psd_xy = ROW(16);
	Index psd_xz = ROW(17);
	Index psd_yz = ROW(18);

	Index a0 = ROW(21);
	Scale(2*M_PI, freq_space, a0);
	Div(psd_zz, Mul(a0, a0, Mul(a0, a0, a0)), a0);
	Table[a0][0] = 0;

	Coord m0 = POS(0, 22);
	Index m0_row = ROW(22);
	FOLD(ADD, Scale(freq/COLS, a0, m0_row), m0);

	Coord m1 = POS(0, 23);
	Index m1_row = ROW(23);
	FOLD(ADD, Scale(freq/COLS, Mul(a0, freq_space, m1_row), m1_row), m1);

	Coord mm1 = POS(0, 24);
	Index mm1_row = ROW(24);
	Scale(freq/COLS, Div(a0, freq_space, mm1_row), mm1_row);
	Table[mm1_row][0] = 0;	// clear nan
	FOLD(ADD, mm1_row, mm1);

	Coord te = DivCell(mm1,m0, POS(0, 25));

	Coord m2 = POS(0, 26);
	Index m2_row = ROW(26);
	FOLD(ADD,
		Scale(freq/COLS,
			Mul(a0,
				Mul(freq_space, freq_space, m2_row),
			 m2_row
			),
			m2_row),
		m2)
	;

	Coord a0_max_coord = MaxCoordReal(a0);
	Coord tp = POS(0, 27);
	SetCell(1.0/(Table[freq_space][a0_max_coord.col]), tp);

	Index denom = ROW(28);
	Sqrt(Mul(psd_zz, Add(psd_xx, psd_yy, denom), denom), denom);

	Index a1 = ROW(29);
	Scale(-1.0, SetImag(0,
		Div(Shift(-1, Mov(psd_xz, a1), a1), denom, a1),
		a1
	), a1);

	Index b1 = ROW(30);
	// Scale(-1.0, (SetImag(0,
	// 	Div(Shift(-1, Mov(psd_yz, b1), b1), denom, b1),
	// 	b1
	// )), b1);
	SetImag(0,
		Div(Shift(-1, Mov(psd_yz, b1), b1), denom, b1),
		b1
	);


	Index denom_2 = Add(psd_xx, psd_yy, ROW(31));

	// printf("a0_max_col %u\n", a0_max_coord.col-2);
	// printf("a0_max_col %u\n", a0_max_coord.col);

	// doesn't match python
	Coord dp = ArcTanCell(
		POS(a0_max_coord.col, b1),
		POS(a0_max_coord.col, a1),
		POS(0, 32)
	);

	Coord Hs = POS(0, 33);
	ScaleCell(4.0, SqrtCell(m0, Hs), Hs);

	Coord Ta = POS(0, 34);
	DivCell(m0, m1, Ta);

	Coord wave_energy_ratio = POS(0, 35);
	DivCell(te, tp, wave_energy_ratio);

	Coord Tz = POS(0, 36);
	SqrtCell(DivCell(m0, m2, Tz), Tz);

	Coord PeakPSD = POS(0, 37);
	MaxCell(a0, PeakPSD);


	// doesn't match python
	Coord Dp_true = POS(0, 38);
	ModCell_immediate(RadToDegreeCell(dp, Dp_true), 360.0, Dp_true);



	Index A2 = ROW(39);
	SetImag(0,
		Div(Sub(psd_xx, psd_yy, A2), denom_2, A2),
	A2);


	Index B2 = ROW(40);
	SetImag(0, Div(Scale(-2.0, psd_xy, B2), denom_2, B2), B2);
	
	

	return 0;
}