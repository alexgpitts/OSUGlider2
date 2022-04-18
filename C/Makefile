

# -fshort-enums is critial, otherwise enums will be integer sized

all:
	gcc -std=c11 -Wno-unused-variable -Wall ./main.c -g -o ./build/main -lm
	./build/main

# all:
# 	gcc -Wno-unused-variable -fsanitize=undefined -fshort-enums -DWORD_ALIGN -O3 -Wall ./main.c -g -o ./build/main
# 	./build/main


valgrind:
	gcc -Wno-unused-variable -fsanitize=undefined -fshort-enums -DWORD_ALIGN -O3 -Wall ./main.c -g -o ./build/main
	valgrind ./build/main



packed:
	gcc -Wno-unused-variable -O3 -Wall ./main.c -g -o ./build/array
	./build/array



clang:
	clang -O0 -Wall ./main.c -g -o ./build/array
	./build/array


gdb:
	gcc -O0 -Wall ./main.c -g -o ./build/array
	gdb ./build/array

