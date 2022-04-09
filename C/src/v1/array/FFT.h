// #include <iostream>
// #include <complex>
#include "../common.h"
#define MAX 200

// using namespace std;

#define M_PI 3.1415926535897932384

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


C128 polar(F64 magnitude, F64 phase) {
	return (C128) {
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
void ordina(C128* f1, int N) {
	C128 f2[MAX];

	for(int i = 0; i < N; i++) {
		f2[i] = f1[reverse(N, i)];
	}
	for(int j = 0; j < N; j++) {
		f1[j] = f2[j];
	}
}

// wouldn't hurt to verify that I did this right
C128 complex_mult(C128 A, C128 B) {
	C128 C = {0};

	C.real = (A.real * B.real) - (A.imag * B.imag);
	C.imag = (A.real * B.imag) + (A.imag * B.real);

	return C;
}

C128 complex_add(C128 A, C128 B) {
	return (C128) {
		.real = A.real + B.real,
		.imag = A.imag + B.imag
	};
}

C128 complex_sub(C128 A, C128 B) {
	return (C128) {
		.real = A.real - B.real,
		.imag = A.imag - B.imag
	};
}


C128 complex_pow(C128 A, F64 B) {
	C128 C = A;
	for (UZ i = 1; i < B; i++) {
		C = complex_mult(A, C);
	}

	return C;
}


// void transform(C128* f, int N) {
// 	ordina(f, N);    //first: reverse order
// 	C128 *W;
// 	W = (C128 *)malloc(N / 2 * sizeof(C128));
// 	W[1] = polar(1., -2. * M_PI / N);
// 	W[0] = (C128) { .real = 1, .imag = 0};
// 	for(int i = 2; i < N / 2; i++) {
// 		W[i] = complex_pow(W[1], i);
// 	}
// 	int n = 1;
// 	int a = N / 2;
// 	for(int j = 0; j < gcc_log2(N); j++) {
// 		for(int i = 0; i < N; i++) {
// 			if(!(i & n)) {
// 				C128 temp = f[i];
// 				// C128 Temp = W[(i * a) % (n * a)] * f[i + n];
// 				C128 Temp = complex_mult(W[(i * a) % (n * a)], f[i + n]);
// 				f[i] = complex_add(temp, Temp);
// 				f[i + n] = complex_sub(temp, Temp);
// 			}
// 		}
// 		n *= 2;
// 		a = a / 2;
// 	}
// 	free(W);
// }


void transform(C128* f, C128*W, int N) {
	ordina(f, N);    //first: reverse order
	W[1] = polar(1., -2. * M_PI / N);
	W[0] = (C128) { .real = 1, .imag = 0};
	for(int i = 2; i < N / 2; i++) {
		W[i] = complex_pow(W[1], i);
	}
	int n = 1;
	int a = N / 2;
	for(int j = 0; j < gcc_log2(N); j++) {
		for(int i = 0; i < N; i++) {
			if(!(i & n)) {
				C128 temp = f[i];
				// C128 Temp = W[(i * a) % (n * a)] * f[i + n];
				C128 Temp = complex_mult(W[(i * a) % (n * a)], f[i + n]);
				f[i] = complex_add(temp, Temp);
				f[i + n] = complex_sub(temp, Temp);
			}
		}
		n *= 2;
		a = a / 2;
	}

	// fix missing complex val... I don't know why this works... soooo... :D
	// printf(">> %d %f %f\n", N, f[0].real, f[N/2].real);
	f[0].imag = f[N/2].real;
}

// void FFT(C128* f, int N, double d) {
// 	transform(f, N);
// 	for(int i = 0; i < N; i++) {
// 		// f[i] *= d; //multiplying by step
// 		f[i] = (C128){
// 			.real = f[i].real * d,
// 			.imag = f[i].imag * d
// 		}; //multiplying by step
// 	}
// }

void FFT(
	C128* in,
	C128* tmp,
	int N
) {
	transform(in, tmp, N);
}

void rFFT(
	C128* in,
	C128* tmp,
	F64* out, 
	int N
) {
	transform(in, tmp, N);

	for (int i = 0 ; i < N; i++) {
		out[i] = in[i].real;
	}
}

