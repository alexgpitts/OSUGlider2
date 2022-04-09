
#include "src/v1/array/data.h"
#include "src/v1/array/code.h"
#include "src/v1/array/meta.h"
#include "src/v1/driver.h"


// // .e r4 !! 1 + (*) // %%  (+) / 1 +

int main(int argc, char const *argv[]) {

	Coord e =
	INC(1,
		FOLD(ADD,
			Reciprocal(
				Scan(MUL,
					Inc(1
					,	Iota(ROW(20))
					,	ROW(20))
				,	ROW(20))
			,	ROW(20))
		,	POS(0,20))
	,	POS(0,20))
	;
	
	// read_csv("./acc.csv");
	read_csv("./067.20201225_1200.20201225_1600_xyz.csv");

	process();

	print_tab(20, ROWS);
	return 0;
}

// #define MEOW_FFT_IMPLEMENTATION
// #include "src/v1/array/meow_fft.h"

// #include <malloc.h>

// int main(int argc, char** argv) {
// 	// (void) argv;
// 	// (void) argc;

// 	unsigned          N    = 1024;
// 	float*            in   = malloc(sizeof(float) * N);
// 	Meow_FFT_Complex* out  = malloc(sizeof(Meow_FFT_Complex) * N);
// 	Meow_FFT_Complex* temp = malloc(sizeof(Meow_FFT_Complex) * N);

// 	// prepare data for "in" array.
// 	// ...

// 	size_t workset_bytes = meow_fft_generate_workset_real(N, NULL);
// 	// Get size for a N point fft working on non-complex (real) data.

// 	Meow_FFT_Workset_Real* fft_real =
// 	(Meow_FFT_Workset_Real*) malloc(workset_bytes);

// 	meow_fft_generate_workset_real(N, fft_real);

// 	meow_fft_real(fft_real, in, out);
// 	// out[0].r == out[0  ].r
// 	// out[0].j == out[N/2].r

// 	// meow_fft_real_i(fft_real, in, temp, out);
// 	// result is not scaled, need to divide all values by N

// 	free(fft_real);
// 	free(out);
// 	free(temp);
// 	free(in);
// }