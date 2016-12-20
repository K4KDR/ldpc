CC=g++
#CFLAGS=-Wformat=0 -ggdb -O0
CFLAGS=-Wformat=0 -O2

LIBS = -lntl
PREFIX=/usr/local

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< 

all: ldpc_compute_generator

ldpc_compute_generator: ldpc_compute_generator.o
	$(CC) $(CFLAGS) $(LIBS) ldpc_compute_generator.cpp -o $@

clean:
	rm -f *.o
	rm -f ldpc_compute_generator

install:
	cp ldpc_compute_generator $(PREFIX)/bin/
	
uninstall:
	rm -rf $(PREFIX)/bin/ldpc_compute_generator
