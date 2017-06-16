build: tema2.c utils.c utils.h
	gcc -O -ftree-vectorize -funroll-loops -Wall -o tema2 utils.c tema2.c
run: build
	./tema2
verbose: tema2.c utils.c utils.h
	gcc -D VERBOSE -Wall -o tema2 utils.c tema2.c

clean:
	-rm -rf tema2 a.out
