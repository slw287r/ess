ARCH := $(shell arch)
ifeq ($(ARCH),x86_64)
CFLAGS=-Wall -O3 -static
LIBS=-lhts -lpthread -ldeflate -lisal -lcairo -lfontconfig -lpixman-1 -lfreetype -lzip -lbz2 -llzma -lpng -lm -ldl
else
CFLAGS=-Wall -O3 -march=native
LIBS=-lhts -lisal -lm -lcairo
endif

CC=cc

SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

ess: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

test: ess
	$^ -i ../dat/toy.bam -r ../dat/ref.fna.gz

clean:
	rm -f ess *.o
