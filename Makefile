.PHONY: all clean

OBJ = src/matrix.o src/bam_filt.o src/count_motif.o src/dna.o src/read_count.o
SRC = src/matrix.c src/bam_filt.c src/count_motif.c src/dna.c src/read_count.c

CFLAG=-Wall -O3

INCLUDE=-Ihtslib -Isrc

all: bin/matlock

htslib/libhts.a:
	cd htslib && $(MAKE)

bin/matlock: bin htslib/libhts.a src/matlock.c $(OBJ)
	gcc $(CFLAG) src/matlock.c $(OBJ) $(INCLUDE) htslib/libhts.a  -lbz2 -llzma -lz -lm -lpthread -lgsl -lgslcblas -o bin/matlock

$(OBJ) : src/%.o : src/%.c htslib/libhts.a
	gcc $(CFLAG) $(INCLUDE) -c $< -o $@ 


clean:
	rm -rf src/*.o; rm -rf bin/matlock
