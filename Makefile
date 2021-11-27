CC=gcc

IDIR=include
CFLAGS=-I$(IDIR) -lm -lz -O3 -march=native -lpthread -llzma -lbz2

_DEPS = data_t.h destroy.h loglik.h reading.h sampler.h usage.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = chromopainter.o destroy.o loglik.o main.o reading.o sampler.o usage.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

SRC=src

BIN=bin

ODIR=obj

$(ODIR)/%.o: $(SRC)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BIN)/chromopainter: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -f $(BIN)/*
