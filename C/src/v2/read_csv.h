

#pragma once

#include "common.h"
#include "array/data.h"

#define MAX_LINE 1<<10
char buffer[MAX_LINE] = {0};


#define SKIP_COLS 1

#define TOKEN_MAX (INPUTS + SKIP_COLS)
char tokens[TOKEN_MAX][MAX_LINE] = {0};


void read_csv(char*filename)	{

	UZ skip_lines = 1;

	FILE*file = fopen(filename, "r");
	if (!file) {
      perror("Error opening file");
      exit(1);
	}

	Index col = 0;
	while (fgets(buffer, MAX_LINE, file)) {
		// puts(buffer);

		UZ T = 0;
		UZ C = 0;
		for (UZ i = 0; (i < MAX_LINE) && (T < TOKEN_MAX); i++) {
			char c = buffer[i];

			switch (c) {
			case '\n':
				T = 0;
				C = 0;
				if (skip_lines > 0) {
					skip_lines--;
					goto SKIP;
				}
				else {
					goto NEXT_LINE;
				}
			break;

			case ' ': case '\t':
			case '\r':
			break;

			case ',':
				T += 1;
				C = 0;
			break;

			default:
				tokens[T][C++] = c;
			}
		}

		NEXT_LINE:


		// skip timestamp
		for (UZ i = SKIP_COLS; i < TOKEN_MAX; i++) {
			Input[i-SKIP_COLS][col] = (float) atof(tokens[i]);
			// printf("%s %f\n", tokens[i], Table[i][col]);
		}

		col += 1;

		SKIP:
		if (col >= INPUT_MAX) {
			break;
		}
		memset(buffer, 0, MAX_LINE);
		for (UZ i = 0; i < TOKEN_MAX; i++) {
			memset(tokens[i], 0, MAX_LINE);
		}
	}


	fclose(file);



}



