#include <ldpc/encoder.h>
#include <ldpc/decoder.h>

#include <stdlib.h>
#include <math.h>
#include <random>

#define LDPC_CODES_PATH "/home/v1tzl1/work/MOVE/LDPC/code/codes/"

/** Convert received BPSK symbol (-1.0=1, 1.0=0) to LLR ( log10(p(0)/p(1)) ) */
ldpc::softbit_t sym2llr(ldpc::softbit_t x, const float sigma) {
    // p_1(x) = 0.5 + 0.25*( erf( (x-1)/scale ) + erf( (x+1)/scale ) )      with scale = sqrt(2)*sigma
    
    const ldpc::softbit_t p1 = 0.5f-0.5f*tanh(x/sigma/sigma);
    const ldpc::softbit_t llr = log10( (1.0f-p1)/p1 );
    
    //printf("Symb: %12f => p(1)==%12f => LLR=%12f\n", x, p1, llr);
    return llr;
}

/** count nuber of active bits in byte */
uint8_t count_bits(uint8_t b) {
    uint8_t count = 0;
    for(uint8_t i=0; i<8; i++) {
        count += (b>>i & 0x01) ? 1 : 0;
    }
    return count;
}

void test(const char *code_name, ldpc::systematic::systematic_t systype, ldpc::puncturing::conf_t pconf, float sigma) {
    
    //
    //// Generate paths to generator and parity check matrix
    //
    const uint64_t STRLEN_GEN = strlen(LDPC_CODES_PATH)+strlen(code_name)+strlen("g.gen");
    const uint64_t STRLEN_PAR = strlen(LDPC_CODES_PATH)+strlen(code_name)+strlen(".a");
    
    char filename_gen[STRLEN_GEN] = LDPC_CODES_PATH;
    strcat(filename_gen, "g");
    strcat(filename_gen, code_name);
    strcat(filename_gen, ".gen");
    
    char filename_par[STRLEN_PAR] = LDPC_CODES_PATH;
    strcat(filename_par, code_name);
    strcat(filename_par, ".a");
    
    printf("  using generator matrix from    %s\n", filename_gen);
    printf("  using parity check matrix from %s\n", filename_par);
    
    //
    //// Create encoder and decoder
    //
    ldpc::encoder enc(filename_gen, systype, &pconf);
    ldpc::decoder dec(filename_par, systype, &pconf);
    ldpc::decoder::metadata_t meta;
    
    //
    //// Define dimmensions
    //
    const uint64_t K_bytes = enc.get_num_input();
    const uint64_t M_punct_bytes = enc.get_num_output();
    const uint64_t K = dec.get_num_output();
    const uint64_t M_punct = dec.get_num_input();
    
    if((uint64_t)ceil(K/8.0) != K_bytes) {
        fprintf(stderr, "Number of information bits is inconsistent between encoder (%u Bytes = %u bits) and decoder (%u bits)\n", K_bytes, K_bytes*8, K);
        exit( EXIT_FAILURE );
    }
    
    if((uint64_t)ceil(M_punct/8.0) != M_punct_bytes) {
        fprintf(stderr, "Number of total bits is inconsistent between encoder (%u Bytes = %u bits) and decoder (%u bits)\n", M_punct_bytes, M_punct_bytes*8, M_punct);
        exit( EXIT_FAILURE );
    }
    
    //
    //// Allocate memory
    //
    uint8_t data_send[K_bytes];
    uint8_t data_encoded[M_punct_bytes];
    ldpc::softbit_t sym_send[M_punct];
    bool mask_ber[M_punct];
    ldpc::softbit_t sym_recv[M_punct];
    ldpc::softbit_t sbits_recv[M_punct];
    ldpc::softbit_t sbits_dec[K];
    uint8_t data_dec[K_bytes];
    
    uint64_t ber_counter = 0;
    uint64_t sym_ber_counter_total = 0;
    uint64_t sym_ber_counter_bits = 0;
    
    size_t i, j;
    uint8_t tmp_byte, tmp_bit;
    std::default_random_engine gen;
    std::uniform_int_distribution<uint8_t> rng_uint8(0,255);
    std::normal_distribution<ldpc::softbit_t> rng_softbit_normal(0.0f,sigma);
    //
    //// Generate send data
    //
    for(i=0; i<K_bytes; i++) {
        data_send[i] = rng_uint8(gen);
    }
    
    //
    //// Encode data
    //
    enc.encode(data_encoded, data_send);
    
    //
    //// Modulate data
    //
    for(i=0; i<M_punct_bytes; i++) {
        for(j=0; j<8; j++) {
            sym_send[i*8+j] = (data_encoded[i] & (0x01<<(7-j))) ? -1.0f : 1.0f;
        }
    }
    
    //
    //// Add noise
    //
    for(i=0; i<M_punct; i++) {
        sym_recv[i] = sym_send[i] + rng_softbit_normal(gen);
        mask_ber[i] = (sym_recv[i]*sym_send[i] < 0.0f);
    }
    
    //
    //// Detector
    //
    for(i=0; i<M_punct; i++) {
        sbits_recv[i] = sym2llr(sym_recv[i], sigma);
    }
    
    //
    //// Decoder
    //
    //dec.decode(sbits_dec, sbits_recv, &meta, "/tmp/decoder_tree.txt");
    dec.decode(sbits_dec, sbits_recv, &meta);
    
    /*
    for(i=0; i<100; i++) {
        printf("% 4u: %12f => %12f\n", i, sbits_recv[i], sbits_dec[i]);
    }
    exit( EXIT_SUCCESS );
    //*/
    
    //
    //// Convert and compare data
    //
    for(i=0; i<K_bytes; i++) {
        tmp_byte = 0x00;
        for(j=0; j<8; j++) {
            // Recover bit
            tmp_bit = (sbits_dec[i*8+j] < 0.0f) ? 0x01 : 0x00;
            
            // Set bit at right position
            tmp_byte ^= tmp_bit << (7-j);
        }
        data_dec[i] = tmp_byte;
        ber_counter += count_bits(data_send[i]^data_dec[i]);
    }
    
    //
    //// Status output
    //
    
    if(meta.num_corrected > 0) {
        printf("Send => Encoded => Decoded\n");
        for(i=0; i<K_bytes; i++) {
            printf("Bit % 4u - % 4u: %02X => %02X => %02X\n", i*8, i*8+7, data_send[i], data_encoded[i], data_dec[i]);
        }
    }
    
    printf("Decoding %s after %u iterations. %u bits corrected.\n", meta.success ? "SUCCESSFUL" : "FAILED", meta.num_iterations, meta.num_corrected);
    
    for(i=0; i<K; i++) {
        sym_ber_counter_bits += mask_ber[i] ? 1 : 0;
    }
    for(i=0; i<M_punct; i++) {
        sym_ber_counter_total += mask_ber[i] ? 1 : 0;
    }
    printf("Simulation introduced %u BERs in total.\n%u bits were changed in information part, send and decoded information differ in %u bits.\n", sym_ber_counter_total, sym_ber_counter_bits, ber_counter);
    
    
}

int main(void) {
    
    printf("Testing rate 1/2 k=1024 block code without noise.\n");
    test("AR4JA_r12_k1024n", ldpc::systematic::FRONT, ldpc::puncturing::conf_t(ldpc::puncturing::BACK, 512, NULL), 0.0f);
    
    /*
    printf("\n\n");
    printf("Testing rate 1/2 k=1024 block code with sigma 0.1.\n");
    test("AR4JA_r12_k1024n", ldpc::systematic::FRONT, ldpc::puncturing::conf_t(ldpc::puncturing::BACK, 512, NULL), 0.1f);
    ///*/
    
    printf("Finished.\n");
}
