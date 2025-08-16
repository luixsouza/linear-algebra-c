CC=gcc
CFLAGS=-O2 -std=c11 -Wall -Wextra -Iinclude
LDFLAGS=
OBJ=src/main.o src/linalg.o src/parser.o src/io_utils.o

all: linalg_app

linalg_app: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

src/main.o: src/main.c include/linalg.h include/parser.h include/io_utils.h
src/linalg.o: src/linalg.c include/linalg.h
src/parser.o: src/parser.c include/parser.h include/linalg.h
src/io_utils.o: src/io_utils.c include/io_utils.h

clean:
	rm -f $(OBJ) linalg_app
