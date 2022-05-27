

#include "src/v2/array/data.h"
#include "src/v2/array/code.h"
#include "src/v2/array/meta.h"
#include "src/v2/driver.h"



int main(int argc, char const *argv[]) {

	// used for testing on desktop, replaced by ece implementation
	read_csv("./067.20201225_1200.20201225_1600_xyz.csv");


	// entry point to processing pipeline
	process(
		INPUT_MAX,				// length of input array
		Input[0],				// x acceleration input array
		Input[1],				// y acceleration input array
		Input[2],				// z acceleration input array
		1.2799999713897705	// frequency
	);

	// rows, columns. determines the size of the data window viewable in console, not the size of the actual memory table (this is defined in v2/array/data.h)
	print_table(20, 41);	
	return 0;
}
