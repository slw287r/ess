ARCH := $(shell arch)
ifeq ($(ARCH),x86_64)
CFLAGS=-Wall -O3 -static
LIBS=-lhts -lpthread -ldeflate -lisal -lcairo -lfontconfig -lpixman-1 -lfreetype -lexpat -lzip -lbz2 -llzma -lpng -lm -ldl -lz
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
	$^ -i ../dat/toy.bam -r ../dat/ref.fna.gz -p ../dat/toy.png -s 'E. coli'
	$^ -i ../dat/toy.bam -r ../dat/ref.fna.gz -p ../dat/toy.svg -s 'E. coli'

clean:
	rm -f ess *.o ../dat/toy.png ../dat/toy.svg
