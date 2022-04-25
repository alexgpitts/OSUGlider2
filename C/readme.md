# Ocean wave measurements (C Implementation - On Glider)


## Prerequisites 
Assumes GCC compiler\
Make runs with C11

## Basics

Most functionality assume access to a global table of floating point values called `Table`.\
`Table` is defined in **/src/v1/array/data.h**\
`Table` named as such because it is a rectangular array or equally sized float arrays.\
from data.h:
```
F32 Table[ROWS][COLS] = {0};
```
`ROWS` controls the number of rows available to the program, and `COLS` controls the length of each row. \
Inputs and Outputs to functions are generally either an `Index` or a `Coord`. An `Index` is simple unsigned integer used to specify a row of `Table`. Thus, `Index` is defined as:
```
typedef u32 Index;
```
 `Coord` is used to specify individual values in Table. Thus `Coord` is defined as:
 ```
typedef struct { Index row; Index col } Coord;
 ```

In much of the code you will see two helper macros `ROW` and `POS`. These are used to help annotate the intention of the code.
\
`ROW(x) => Index x`\
`POS(x,y) => Coord { row = y, col = x }`
\
\
Most functions have the interface:
```
Index function(Index source_row_1, ... Index source_row_N, Index target_row)
```
Functions take as input the index of zero, one, or more rows in the table.\
Functions return as outout the index of the row wherein the result was written.
\
\
For example, given the operation:
```
A = B + C - D
```
This expression must be written such that the target of each operation is explicit.
```
B + C → A - D → A
```
Which, when written using the Array.h functions:
```
Index A = ROW(0);
Index B = ROW(1);
Index C = ROW(2);
Index D = ROW(3);

Sub(
	Add(B, C, A),		// Add returns Index A
	D,
	A
);
```
*Every function in array is written such that the input source arrays are not mutated as a result of the operation.*

## File Tour
`Main` calls three functions: `read_csv`, `process`, and `print_table`

`read_csv` exects as input a well formatted CSV files with following header:
```
t, x, y, z
```
`read_csv` expects this format, because this is the format of the output of the Python program (described below).\
The header is discarded, then each row is parsed into 4 data points. The time (t) value is discarded.
Each the x, y, and z columns are then translated from ascii to 32-bit floats.
\
`read_csv` writes the
- `x` column to `ROW(1)`
- `y` column to `ROW(2)`
- `z` column to `ROW(3)`
 
\
\
Before describing the function `process` which concerns the bulk of the programs functionality, I will describe `print_table` which is used for viewing the state of the program.

`print_table` is defined in **/src/v1/array/meta.h**
```
print_table(max_columns, max_rows)
```
`max_columns` and `max_rows` specify the size of the view window into the global `Table`.

\
\
`process` runs after read_csv and assumes that rows 1, 2, and 3 are already filled with acceleration data (x,y,z).
```
TABLE
Row  0 ← Time Series

Row  1 ← X acceleration data
Row  2 ← Y acceleration data
Row  3 ← Z acceleration data

Row  4 ← Nothing

Row  5 ← Rolling mean of X acceleration data
Row  6 ← Rolling mean of Y acceleration data
Row  7 ← Rolling mean of Z acceleration data

Row  8 ← Nothing

Row  9 ← FFT(X), FFT of Row 5
Row 10 ← FFT(Y), FFT of Row 6
Row 11 ← FFT(Z), FFT of Row 7

Row 12 ← Nothing

Row 13 ← PSD(X, X), PSD of Row 9  and Row  9
Row 14 ← FFT(Y, Y), PSD of Row 10 and Row 10
Row 15 ← FFT(Z, Z), PSD of Row 11 and Row 11

Row 16 ← FFT(X, Y), PSD of Row  9 and Row 10
Row 17 ← FFT(X, Z), PSD of Row  9 and Row 11
Row 18 ← FFT(Y, Z), PSD of Row 10 and Row 11

Row 19 ← Nothing

Row 20 ← Frequency Space

Row 21 ← a0

Pos 0,22 ← m0
Pos 0,23 ← m1
Pos 0,24 ← mm1
Pos 0,25 ← te
Pos 0,26 ← m2
Pos 0,27 ← tp

Row 28 ← denom (denominator)
Row 29 ← a1
Row 30 ← b1

Row 31 ← denom2 (denominator 2)

Pos 0,32 ← dp
Pos 0,33 ← Hs
Pos 0,34 ← Ta
Pos 0,35 ← wave energy ratio
Pos 0,36 ← Tz
Pos 0,37 ← PeakPSD
Pos 0,38 ← dp true

Row 39 ← A2
Row 40 ← B2

```

## Tutorial:
Run python program to generate CSV file from CDIP data
```
python .\Driver.py ".\ncFiles\067.20201225_1200.20201225_1600_output.nc"
```

