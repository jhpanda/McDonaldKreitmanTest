
CFLAGS = -O3 -w -funroll-all-loops
LDFLAGS = -lm -fopenmp
CC = gcc

all : ../bin/fisher_exact_test ../bin/cmktest ../bin/cmktest_list

../bin/fisher_exact_test : fisher_exact_test.c fisher_exact.c
	$(CC) ${CFLAGS} ${LDFLAGS} -o $@ $^
../bin/cmktest : cmktest.c fasta.c polymorph.c fisher_exact.c
	$(CC) ${CFLAGS} ${LDFLAGS} -o $@ $^
../bin/cmktest_list : cmktest_list.c fasta.c polymorph.c fisher_exact.c
	$(CC) ${CFLAGS} ${LDFLAGS} -o $@ $^

clean:
	rm -f *.o
