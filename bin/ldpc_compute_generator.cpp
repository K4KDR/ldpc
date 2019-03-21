#include <NTL/GF2.h>
#include <NTL/mat_GF2.h>
#include <vector>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>

void error(const char* err_str) {
    fprintf(stderr, "ERROR: ");
    fprintf(stderr, "%s", err_str);
    fprintf(stderr, "\n");
    exit( EXIT_FAILURE );
}

template <typename T> std::vector<T> parse_digits_from_line(FILE *f, const char *line_descr, size_t num, bool ignore_zeros) {
    std::vector<T> ret;
    
    char *line = NULL;
    char *line_alloc;
    char *end;
    size_t len = 0;
    ssize_t read;
    
    read = getline(&line, &len, f);
    if(read == -1) {
        fprintf(stderr, "EOF reached while reading %s.\n", line_descr);
        exit( EXIT_FAILURE );
    }
    line_alloc = line;
    
    //fprintf(stdout, "Parsing line: '%s'\n", line);
    for (long long i = strtoll(line, &end, 10); line != end; i = strtol(line, &end, 10)) {
        //fprintf(stdout, "'%.*s' -> %lld", (int)(end-line), line, i);
        line = end;
        if (errno == ERANGE){
            fprintf(stderr, "Read number is out of range.\n");
            exit( EXIT_FAILURE );
        } else if(i<0) {
            fprintf(stderr, "Number is negative.\n");
            exit( EXIT_FAILURE );
        } else {
            // Value in range
            
            if(!ignore_zeros || i != 0) {
                ret.push_back( (T) i );
            }
        }
    }
    
    free(line_alloc);
    
    if(ret.size() != num) {
        fprintf(stderr, "%lu numbers read in %s line, but %lu were expected.\n", ret.size(), line_descr, num);
        exit( EXIT_FAILURE );
    }
    
    return ret;
}

typedef unsigned long val_t;

void compute_generator(NTL::Mat<NTL::GF2> &G, const char* alist_file) {
    NTL::Mat<NTL::GF2> Q, P;
    
    ////
    //// Read alist file and output matrix as wide [Q P] matrix
    ////
    // M > N
    // K := M-N
    // Q: NxK matrix
    // P: NxN matrix
    
    FILE* f = fopen(alist_file, "r");
    
    if(!f) {
        fprintf(stderr, "Cannot open file %s\n", alist_file);
        exit( EXIT_FAILURE );
    }
    
    std::vector<val_t> buf;
    std::vector<val_t> num_n, num_m;
    
    // Read N M
    buf = parse_digits_from_line<val_t>(f, "dimensions", 2, false);
    
    const val_t alist_N = buf[0];
    const val_t alist_M = buf[1];
    
    bool do_trans = (alist_N > alist_M);
    
    const val_t N = do_trans ? alist_M : alist_N;
    const val_t M = do_trans ? alist_N : alist_M;
    
    if(N > M) {
        error("Transpose logic failed.");
    }
    
    const val_t K = M-N;
    Q.SetDims(N, K);
    P.SetDims(N, N);
    
    
    // Read biggest_num_n biggest_num_m (ignored)
    buf = parse_digits_from_line<val_t>(f, "maximum elements", 2, false);
    
    // Read num_n
    num_n = parse_digits_from_line<val_t>(f, "nlist count", alist_N , false);
    
    // Read num_m
    num_m = parse_digits_from_line<val_t>(f, "mlist count", alist_M , false);
    
    // If we have to transpose, skip nlist
    if(do_trans) {
        for(size_t i=0; i<alist_N ;i++) {
            buf = parse_digits_from_line<val_t>(f, "nlist", num_n[i], true);
        }
        num_n = num_m;
    }
    
    // Now read either nlist or mlist as new nlist (num_n was updated before)
    for(size_t i=0; i<N ;i++) {
        buf = parse_digits_from_line<val_t>(f, "n/m-list", num_n[i], true);
        
        for(size_t j=0; j<num_n[i]; j++) {
            if(buf[j] <= K) {
                // Entry belongs to Q matrix
                Q[i][buf[j]-1] = 1;
            } else {
                // Entry belongs to P matrix
                P[i][buf[j]-K-1] = 1;
            }
        }
    }
    
    // If we have not transposed, read mlist now
    if(!do_trans) {
        for(size_t i=0; i<alist_M ;i++) {
            buf = parse_digits_from_line<val_t>(f, "mlist", num_m[i], true);
        }
    }
    
    // Read until EOF, ignore spaces and newlines
    int c;
    while(true) {
        c = getc(f);
        
        if(c == -1) {
            // EOF
            break;
        } else if( c == ' ' || c == '\n' || c == '\r' ) {
            continue;
        } else {
            fprintf(stderr, "alist file contains illegal character '%c' after matrix is read in.\n", c);
            exit( EXIT_FAILURE );
        }
    }
    
    fclose(f);
    // Reading of alist finished
    
    // print alist_ matrix
    /*
    printf("[Q | P] = \n");
    for(size_t i=0; i<N; i++) {
        for(size_t j=0; j<K; j++) {
            printf("%1d ", NTL::IsOne(Q[i][j]) ? 1 : 0);
        }
        printf("| ");
        for(size_t j=0; j<N; j++) {
            printf("%1d ", NTL::IsOne(P[i][j]) ? 1 : 0);
        }
        printf("\n");
    }
    //*/
    
    NTL::Mat<NTL::GF2> Pi;
    NTL::inv(Pi, P);
    
    // Free P
    P.kill();
    
    G = Pi * Q;
}

void write_matrix_txt(FILE *out, NTL::Mat<NTL::GF2> &G) {
    const val_t N = G.NumRows();
    const val_t K = G.NumCols();
    
    for(size_t i=0; i<N; i++) {
        for(size_t j=0; j<K; j++) {
            fprintf(out, "%1d ", NTL::IsOne(G[i][j]) ? 1 : 0);
        }
        fprintf(out, "\n");
    }
}

void write_matrix_bin(FILE *out, NTL::Mat<NTL::GF2> &G) {
    
    // G is a NxK matrix
    const val_t N = G.NumRows();
    const val_t K = G.NumCols();
    
    if(K%8 != 0) {
        error("Number of information bits is not a multiple of 8 bit. Cannot write binary compressed form.");
    }
    const val_t num_bytes = K/8;
    uint8_t buf[8];
    uint8_t tmp;
    
    // Write N as 8 byte
    buf[0] = (uint8_t) (0x00000000000000FF & (N>>7*8));
    buf[1] = (uint8_t) (0x00000000000000FF & (N>>6*8));
    buf[2] = (uint8_t) (0x00000000000000FF & (N>>5*8));
    buf[3] = (uint8_t) (0x00000000000000FF & (N>>4*8));
    buf[4] = (uint8_t) (0x00000000000000FF & (N>>3*8));
    buf[5] = (uint8_t) (0x00000000000000FF & (N>>2*8));
    buf[6] = (uint8_t) (0x00000000000000FF & (N>>1*8));
    buf[7] = (uint8_t) (0x00000000000000FF & (N>>0*8));
    fwrite(buf, sizeof(uint8_t), 8, out);
    printf("N = %lu = %2X%2X%2X%2X%2X%2X%2X%2X\n", N, buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
    
    // Write K as 8 byte
    buf[0] = (uint8_t) (0x00000000000000FF & (K>>7*8));
    buf[1] = (uint8_t) (0x00000000000000FF & (K>>6*8));
    buf[2] = (uint8_t) (0x00000000000000FF & (K>>5*8));
    buf[3] = (uint8_t) (0x00000000000000FF & (K>>4*8));
    buf[4] = (uint8_t) (0x00000000000000FF & (K>>3*8));
    buf[5] = (uint8_t) (0x00000000000000FF & (K>>2*8));
    buf[6] = (uint8_t) (0x00000000000000FF & (K>>1*8));
    buf[7] = (uint8_t) (0x00000000000000FF & (K>>0*8));
    fwrite(buf, sizeof(uint8_t), 8, out);
    printf("K = %lu = %2X%2X%2X%2X%2X%2X%2X%2X\n", K, buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
    
    // Write byte masks
    for(size_t i=0; i<N; i++) {
        
        // Write mask for parity check i
        for(size_t j=0; j<num_bytes; j++) {
            
            // Write j-th byte mask of parity check i
            tmp = 0x00;
            buf[0] = 0x00;
            
            for(uint8_t k=0; k<8; k++) {
                if(NTL::IsOne(G[i][j*8+k])) {
                    // if k-th bit (0:MSB, 7:LSB) is one, set it in mask
                    
                    // Set bit to one
                    tmp = 0x01 << (7-k);
                    
                    // Toggle bit in buffer
                    buf[0] ^= tmp;
                }
            }
            fwrite(buf, sizeof(uint8_t), 1, out);
        }
    }
}


int main(int argc, char* argv[]) {
    
    // Check number of provided arguments and print help if number is not correct.
    if(argc < 2 || argc > 4) {
        fprintf(stdout, "Usage: %s alist_file <output_txt> <output_gen>\n\n", argv[0]);
        fprintf(stdout, "Read parity check matrix from `alist_file`. If matrix is tall (i.e. has more rows than columns) it is transposed to always yield a NxM matrix with M >= N and K=M-N >= 0. This matrix is split into [Q P] with Q: NxK and P: NxN. The generator matrix is then found by inverting P and computing G=inv(P)*Q. The generator has dimensions NxK.\n\n");
        fprintf(stdout, "The computed generator is written in ASCII format (N lines of K '0' or '1's, separated by spaces). If `output_txt` is provided, the matrix is written into this file, otherwise it is printed on stdout. The output can be read in by matlab's load command by using the '-ascii' option.\n\n");
        fprintf(stdout, "If `output_gen` is provided, a binary form of the generator matrix is produced and written into this file. The binary form consists of 2*8 bytes containing N and K followed by N*ceil(K/8) bytes. The first ceil(K/8) bytes belong to the first row of G, the next ceil(K/8) bytes to the second row, etc. The first byte of each row contains the first 8 columns of G, the second the next 8 columns, etc. The MSB belongs to the column with the lowest index covered by the byte. E.g: the MSB for the first byte belongs to the first column of G.\n");
        return 1;
    }
    
    // Compute generator matrix G
    NTL::Mat<NTL::GF2> G;
    printf("Computing generator matrix\n");
    compute_generator(G, argv[1]);
    printf(" ... complete\n\n");
    
    // Produce ASCII output
    FILE *out_txt;
    if(argc >= 3) {
        printf("Writing ASCII generator to  %s\n", argv[2]);
        out_txt = fopen(argv[2], "w");
    } else {
        printf("Computed ASCII generator:\n");
        out_txt = stdout;
    }
    write_matrix_txt(out_txt, G);
    if(argc >= 3) {
        fclose(out_txt);
    }
    printf(" ... complete\n\n");
    
    // Produce binary output
    if(argc >= 4) {
        printf("Writing binary generator to %s\n", argv[3]);
        
        FILE *out_bin = fopen(argv[3], "wb");
        write_matrix_bin(out_bin, G);
        fclose(out_bin);
        printf(" ... complete\n\n");
    }
}
