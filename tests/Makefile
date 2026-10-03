CC ?= gcc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Werror -pedantic
all: test_engine
test_engine: test_engine.c ../src/engine.c ../src/engine.h
	$(CC) $(CFLAGS) test_engine.c ../src/engine.c -lm -o $@
run: test_engine
	./test_engine
clean:
	rm -f test_engine
