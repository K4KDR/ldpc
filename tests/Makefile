CC=g++
CFLAGS=-fPIC -Wall -Wformat=0 -O0 -ggdb -std=c++11
#CFLAGS=-fPIC -Wall -Wformat=0 -O2 -std=c++11

%.o: %.cpp $
	$(CC) $(CFLAGS) -c -o $@ $< 

all: test_encoder test_decoder test_chain

test_encoder: test_encoder.o
	$(CC) $(CFLAGS) -lldpc test_encoder.cpp -o $@

test_decoder: test_decoder.o
	$(CC) $(CFLAGS) -lldpc test_decoder.cpp -o $@

test_chain: test_chain.o
	$(CC) $(CFLAGS) -lldpc test_chain.cpp -o $@

clean:
	rm -f *.o
	rm -f test_encoder
	rm -f test_decoder
	rm -f test_chain
