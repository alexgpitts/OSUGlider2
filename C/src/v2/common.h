#pragma once


#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#define M_PI 3.1415926535897932384

#define NUMARGS(...)  (sizeof((int[]){__VA_ARGS__})/sizeof(int))
#define SUM(...)  (sum(NUMARGS(__VA_ARGS__), __VA_ARGS__))

#define KB (1<<10)
#define MB (1<<20)

typedef uint8_t   U8;
typedef uint16_t  U16;
typedef uint32_t  U32;
typedef uint64_t  U64;

typedef int8_t    I8;
typedef int16_t   I16;
typedef int32_t   I32;
typedef int64_t   I64;

// typedef __int128 i128;
// typedef unsigned __int128 u128;

typedef volatile U32 R32;
typedef volatile U32 Reg32;

typedef volatile U64 R64;
typedef volatile U64 Reg64;

typedef float F32;
typedef double F64;

typedef size_t UZ;
typedef U16 Index;
typedef struct { Index col; Index row; } Coord;


typedef struct { F32 real; F32 imag; } C64;
typedef struct { F64 real; F64 imag; } C128;




#define x8_size  (1)
#define x16_size (2)
#define x32_size (4)
#define x64_size (8)

#define Index_size sizeof(Index)
#define Coord_size sizeof(Coord)




typedef struct {	Index x; Index y; Index z; } XYZ;


Index inline MIN(Index a, Index b) {
	return a < b ? a : b;
}
Index inline MAX(Index a, Index b) {
	return a > b ? a : b;
}

Index UPPER_ROUND(UZ num, UZ base) {
	return (num + base - 1) / base * base;
}



#define ROW(x) (x)
// #define COL(x) (x)
#define POS(x,y) ((Coord) { .col = x, .row = y })

#define mXYZ(a,b,c) ((XYZ) { .x = a, .y = b, .z = c })