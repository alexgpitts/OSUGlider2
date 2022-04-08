

#pragma once

#include "common.h"
#include "array/data.h"

#define MAX_LINE 1<<10
char buffer[MAX_LINE] = {0};

#define TOKEN_MAX 10
char tokens[TOKEN_MAX][MAX_LINE] = {0};

void read_csv
(	char*filename
// ,	XYZ xyz
)	{

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
			goto NEXT_LINE;

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

		for (UZ i = 1; i < TOKEN_MAX; i++) {

			Table[i][col] = (float) atof(tokens[i]);
			// printf("%.15lF %s\n", atof(tokens[i]), tokens[i]);
		}

		col += 1;

		if (col >= COLS) {
			break;
		}
		memset(buffer, 0, MAX_LINE);
	}


	fclose(file);



}



