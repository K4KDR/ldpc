CC=g++
#CFLAGS=-fPIC -Wall -Wformat=0 -O0 -ggdb -I./include/
#CFLAGS=-fPIC -Wall -Wformat=0 -O2 -ggdb -I./include/
CFLAGS=-fPIC -Wall -Wformat=0 -O5 -I./include/

PREFIX=/usr/local/

%.o: %.cpp include/ldpc/%.h include/ldpc/ldpc.h $
	$(CC) $(CFLAGS) -c -o $@ $< 

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
