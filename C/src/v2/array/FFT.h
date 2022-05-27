#pragma once
#include "../common.h"

/**
 * custom FFT implementation
 * no dynamic memory allocation
 * lil_fft is the interface function used by operator 'FFT' in array/code.h
 */

/*function to calculate the log2(.) of int numbers*/
// int log2(int N) {
//   int k = N, i = 0;
//   while(k) {
//     k >>= 1;
//     i++;
//   }
//   return i - 1;
// }


// GCC specific log2
U32 gcc_log2(U32 n){
	return (sizeof(U32)*8 - 1) - __builtin_clz(n);
}





//checking if the number of element is a power of 2
int check(int n) {
	return n > 0 && (n & (n - 1)) == 0;
}


C64 polar(F32 magnitude, F32 phase) {
	return (C64) {
		.real = magnitude * cos(phase),
		.imag = magnitude * sin(phase)
	};
}

//calculating revers number
int reverse(int N, int n) {
	int j, p = 0;
	for(j = 1; j <= gcc_log2(N); j++) {
		if(n & (1 << (gcc_log2(N) - j)))
			p |= 1 << (j - 1);
		}
	return p;
}

//using the reverse order in the array
void ordina(C64* f1, int N) {
	static C64 f2[COLS];

	for(int i = 0; i < N; i++) {
		f2[i] = f1[reverse(N, i)];
	}
	for(int j = 0; j < N; j++) {
		f1[j] = f2[j];
	}
}

// wouldn't hurt to verify that I did this right
C64 complex_mult(C64 A, C64 B) {
	C64 C = {0};

	C.real = (A.real * B.real) - (A.imag * B.imag);
	C.imag = (A.real * B.imag) + (A.imag * B.real);

	return C;
}

C64 complex_add(C64 A, C64 B) {
	return (C64) {
		.real = A.real + B.real,
		.imag = A.imag + B.imag
	};
}

C64 complex_sub(C64 A, C64 B) {
	return (C64) {
		.real = A.real - B.real,
		.imag = A.imag - B.imag
	};
}


C64 complex_pow(C64 A, F32 B) {
	C64 C = A;
	for (UZ i = 1; i < B; i++) {
		C = complex_mult(A, C);
	}

	return C;
}




void transform(C64* f, C64*W, int N) {
	ordina(f, N);    //first: reverse order
	W[1] = polar(1., -2. * M_PI / N);
	W[0] = (C64) { .real = 1, .imag = 0};
	for(int i = 2; i < N / 2; i++) {
		W[i] = complex_pow(W[1], i);
	}
	int n = 1;
	int a = N / 2;
	for(int j = 0; j < gcc_log2(N); j++) {
		for(int i = 0; i < N; i++) {
			if(!(i & n)) {
				C64 temp = f[i];
				// C128 Temp = W[(i * a) % (n * a)] * f[i + n];
				C64 Temp = complex_mult(W[(i * a) % (n * a)], f[i + n]);
				f[i] = complex_add(temp, Temp);
				f[i + n] = complex_sub(temp, Temp);
			}
		}
		n *= 2;
		a = a / 2;
	}

}


// interface function used by operator FFT defined in array/code.h
void lil_FFT(
	C64* in,
	C64* tmp,
	int N
) {
	transform(in, tmp, N);
}

void rFFT(
	C64* in,
	C64* tmp,
	F32* out, 
	int N
) {
	transform(in, tmp, N);

	for (int i = 0 ; i < N; i++) {
		out[i] = in[i].real;
	}
}

