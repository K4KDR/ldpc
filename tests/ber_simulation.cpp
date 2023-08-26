#include <ldpc/encoder.h>
#include <ldpc/decoder.h>

#include <cstdlib>
#include <cmath>
#include <random>

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

int main(int argc, char** argv) {
    if(argc != 6) {
        printf("usage: %s generator_bin parity_alist min_Eb_N0 max_Eb_N0 setp_Eb_N0\n", argv[0]);
        exit( EXIT_FAILURE );
    }
    const char *generator_matrix_file = argv[1];
    const char *parity_matrix_file = argv[2];
    float Eb_N0_start = static_cast<float>(std::atof(argv[3]));
    float Eb_N0_stop = static_cast<float>(std::atof(argv[4]));
    float Eb_N0_step = static_cast<float>(std::atof(argv[5]));
    
    ldpc::systematic::systematic_t systype = ldpc::systematic::FRONT;
    ldpc::puncturing::conf_t pconf = ldpc::puncturing::conf_t(ldpc::puncturing::BACK, 512, NULL);
    
    const size_t min_errors = 10lu;
    const size_t update_rate = 1lu;
    
    //
    //// Create encoder and decoder
    //
    ldpc::encoder enc(generator_matrix_file, systype, &pconf);
    ldpc::decoder dec(parity_matrix_file, systype, &pconf);
    ldpc::decoder::metadata_t meta;
    
    //
    //// Define dimmensions
    //
    const uint64_t K_bytes = enc.get_num_input();
    const uint64_t M_punct_bytes = enc.get_num_output();
    const uint64_t K = dec.get_num_output();
    const uint64_t M_punct = dec.get_num_input();
    
    if((uint64_t)ceil(K/8.0) != K_bytes) {
        fprintf(stderr, "Number of information bits is inconsistent between encoder (%lu Bytes = %lu bits) and decoder (%lu bits)\n", K_bytes, K_bytes*8u, K);
        exit( EXIT_FAILURE );
    }
    
    if((uint64_t)ceil(M_punct/8.0) != M_punct_bytes) {
        fprintf(stderr, "Number of total bits is inconsistent between encoder (%lu Bytes = %lu bits) and decoder (%lu bits)\n", M_punct_bytes, M_punct_bytes*8, M_punct);
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
    
    uint64_t sym_ber_counter_total = 0;
    uint64_t sym_ber_counter_bits = 0;
    
    size_t i, j, k;
    uint8_t tmp_byte, tmp_bit;
    ldpc::softbit_t tmp_softbit;
    std::default_random_engine gen;
    std::uniform_int_distribution<uint8_t> rng_uint8(0,255);
    
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
    //// Run simulation
    //
    const size_t num_steps = static_cast<std::size_t>(std::floor((Eb_N0_stop-Eb_N0_start)/Eb_N0_step));
    float EbN0vec[num_steps];
    double BERvec[num_steps];
    for (i = 0lu; i<num_steps; i++) {
        float EbN0 = Eb_N0_start + static_cast<float>(i)*Eb_N0_step;
        float sigma = std::sqrt((std::pow(10.0f, -EbN0/10.0f)) / 2.0f);
        std::normal_distribution<ldpc::softbit_t> rng_softbit_normal(0.0f,sigma);
        
        size_t num_errors = 0lu;
        size_t num_errors_uncoded = 0lu;
        
        size_t num_bits_total = 0lu;
        size_t num_bits_total_uncoded = 0lu;
        size_t num_iteration = 0lu;
        double ber = -1.0f;
        double ber_uncoded = -1.0f;
        
        printf("%3lu/%3lu: Eb/N0 = %+7.3f dB", i, num_steps, EbN0);
        while (true) {
            
            //
            //// Add noise
            //
            //printf("noise vector: ");
            for(j=0; j<M_punct; j++) {
                tmp_softbit = rng_softbit_normal(gen);
                //printf("%f ", tmp_softbit);
                sym_recv[j] = sym_send[j] + tmp_softbit;
                mask_ber[j] = (sym_recv[j]*sym_send[j] < 0.0f);
                num_errors_uncoded += mask_ber[j];
            }
            
            //
            //// Detector
            //
            for(j=0; j<M_punct; j++) {
                sbits_recv[j] = sym2llr(sym_recv[j], sigma);
            }
            
            //
            //// Decoder
            //
            dec.decode(sbits_dec, sbits_recv, &meta);
            
            //
            //// Convert and compare data
            //
            for(j=0; j<K_bytes; j++) {
                tmp_byte = 0x00;
                for(k=0; k<8; k++) {
                    // Recover bit
                    tmp_bit = (sbits_dec[j*8+k] < 0.0f) ? 0x01 : 0x00;
                    
                    // Set bit at right position
                    tmp_byte ^= tmp_bit << (7-k);
                }
                data_dec[j] = tmp_byte;
                num_errors += count_bits(data_send[j]^data_dec[j]);
            }
            num_bits_total += K_bytes*8u;
            num_bits_total_uncoded += M_punct_bytes*8u;
            
            num_iteration++;
            
            const bool step_complete = num_errors > min_errors;
            
            ber = static_cast<double>(num_errors)/static_cast<double>(num_bits_total);
            ber_uncoded = static_cast<double>(num_errors_uncoded)/static_cast<double>(num_bits_total_uncoded);
            if ((num_iteration-1u) % update_rate == 0lu || step_complete) {
                printf("\r%3lu/%3lu: Eb/N0 = %+7.3f dB, BER=%5.3e, %3lu errors, %8lu bits, %6lu iterations, BER uncoded=%5.3e, %8lu/%8lu", i, num_steps, EbN0, ber, num_errors, num_bits_total, num_iteration, ber_uncoded, num_errors_uncoded, num_bits_total_uncoded);
            }
            
            
            if (step_complete) {
                break;
            }
        }
        
        EbN0vec[i] = EbN0;
        BERvec[i] = ber;
        printf("\n");
    }
    
    printf("Finished.\n");
}
