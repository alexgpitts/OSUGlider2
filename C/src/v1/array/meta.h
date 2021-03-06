#pragma once

#include "../common.h"
#include "data.h"
/**
 * meta.h defines functions for printing the contents of Table in a table-like format
 * print_table is the interface function used in main
 * does not do anything fancy to determine the width of your terminal window, so if
 * it looks a mess, make your terminal window wider or WIDTH smaller
 * output in MS Powershell is frunked up. Best use a unix terminal
 */


// Number of columns printed before moving to 'newline'
#define WIDTH 10

void print_val(Coord p);
void print_row(Index r);

void print_table(Index, Index);
void print_col_headers(Index);
void print_line();
void print_line2();



void print_val(Coord p) {
	printf("│ %-12g ", Table[p.row][p.col]);
}

void print_row(Index r) {
	for (Index c = 0; c < COLS; c++) {
		print_val(POS(c, r));
		printf(" ");
	}
	printf("│\n");
}


// prints the contents of 'Table' to stdin
void print_table(Index max_cols, Index max_rows) {
	for (Index k = 0; k < COLS && k < max_cols; k+=WIDTH) {
		print_line2();
		print_col_headers(k);
		print_line();
		for (Index r = 0; r < ROWS && r < max_rows; r++) {
			printf("r%-2d ", r);
			for (Index c = k; c < k+WIDTH; c++) {
				if (c < COLS) {
					print_val(POS(c, r));
				}
				else {
					printf("│ x            ");
				}
			}
			printf("\n");
		}
		printf("\n");
	}
}

void print_col_headers(Index k) {
	printf("    ");
	for (Index c = k; c < k+WIDTH; c++) {
		printf("│ %-12d ", c);
	}
		printf("│\n");

}

void print_line() {
	printf("────");
	for (Index i = 0; i < WIDTH; i++) {
		printf("┼──────────────");
	}
	printf("┘\n");
}

void print_line2() {
	printf("    ┌");
	for (Index i = 0; i < WIDTH-1; i++) {
		printf("──────────────┬");
	}
	printf("──────────────┐");
	printf("\n");
}