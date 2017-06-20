CC=g++
#CFLAGS=-fPIC -Wall -Wformat=0 -O0 -ggdb -I./include/
#CFLAGS=-fPIC -Wall -Wformat=0 -O2 -ggdb -I./include/
CFLAGS=-fPIC -Wall -Wformat=0 -O5 -I./include/

PREFIX=/usr/local/

ldpc.o: ldpc.cpp include/ldpc/ldpc.h
	$(CC) $(CFLAGS) -c -o ldpc.o ldpc.cpp

encoder.o: encoder.cpp include/ldpc/encoder.h include/ldpc/ldpc.h
	$(CC) $(CFLAGS) -c -o encoder.o encoder.cpp

decoder.o: decoder.cpp include/ldpc/decoder.h include/ldpc/ldpc.h
	$(CC) $(CFLAGS) -c -o decoder.o decoder.cpp

all: libldpc

libldpc: encoder.o decoder.o ldpc.o
	$(CC) -shared $(CFLAGS) decoder.o encoder.o ldpc.o -o libldpc.so

clean:
	rm -f *.o
	rm -f libldpc.so

install:
	cp -R include/ldpc $(PREFIX)include/
	cp libldpc.so $(PREFIX)lib/
	
uninstall:
	rm -rf $(PREFIX)/include/ldpc
	rm -rf $(PREFIX)/lib/libldpc.so
