#include <ldpc/encoder.h>
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace ldpc;

encoder::encoder(const char* generator_file, systematic::systematic_t systype, puncturing::conf_t *punctconf) {
    
    FILE *fgen = fopen(generator_file, "rb");
    if(!fgen) {
        fprintf(stderr, "Cannot open generator file %s\n", generator_file);
        exit( EXIT_FAILURE );
    }
    
    // Read N
    uint64_t tmp = 0;
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (7*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (6*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (5*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (4*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (3*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (2*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (1*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (0*8);
    this->N = tmp;
    
    // Read K
    tmp = 0;
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (7*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (6*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (5*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (4*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (3*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (2*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (1*8);
    tmp |= ((uint64_t) read_byte(fgen, "generator dimensions")) << (0*8);
    this->K = tmp;
    
    // Check puncturing
    if(punctconf->type == puncturing::NONE && punctconf->num_punct != 0) {
        fprintf(stderr, "Puncturing was set to none, but non-zero number of puncturing positions was given.\n");
        exit( EXIT_FAILURE );
    }
    if(this->N <= punctconf->num_punct) {
        fprintf(stderr, "After puncturing %lu positions from the %lu parity checks, no parity check would remain.\n", punctconf->num_punct, this->N);
        exit( EXIT_FAILURE );
    }
    
    this->N_punct = this->N - punctconf->num_punct;
    this->N_punct_bytes = (uint64_t) ceil(this->N_punct/8.0);
    this->K_bytes = (uint64_t) ceil(this->K/8.0);
    
    uint8_t buf;
    this->parity_checks = new uint8_t*[this->N_punct];
    size_t i_local = 0;
    bool punct;
    for(size_t i=0; i<this->N; i++) {
        
        punct = punctconf->is_punctured(i+this->K,this->N+this->K);
        
        if(!punct) {
            this->parity_checks[i_local] = new uint8_t[this->K_bytes];
        }
        
        for(size_t j=0; j<this->K_bytes; j++) {
            buf = read_byte(fgen, "Parity generator byte");
            
            if(!punct) {
                this->parity_checks[i_local][j] = buf;
            }
        }
        
        if(!punct) {
            i_local++;
        }
    }
    if(i_local!=this->N_punct) {
        fprintf(stderr, "Allocated %u parity checks, but code has %u.\n", i_local, this->N_punct);
        exit( EXIT_FAILURE );
    }
    
    this->systype = systype;
}

encoder::~encoder() {
    for(size_t i=0; i<this->N_punct; i++) {
        delete[] this->parity_checks[i];
    }
    delete[] this->parity_checks;
}

uint64_t encoder::get_num_input(void) const {
    return this->K_bytes;
}

uint64_t encoder::get_num_output(void) const {
    if(this->systype == systematic::NONE) {
        return this->N_punct_bytes;
    } else {
        return this->N_punct_bytes + this->K_bytes;
    }
}

void encoder::encode(uint8_t *out, const uint8_t *input) {
    
    
    size_t par_ofst = (this->systype == systematic::FRONT) ? this->K_bytes : 0;
    
    uint8_t parity_buf;
    uint8_t output_buf = 0;
    bool flag_buf;
    for(size_t i=0; i<this->N_punct_bytes; i++) {
        output_buf = 0;
        for(uint8_t j=0; j<8 && i*8+j<this->N_punct; j++) {
            parity_buf = 0;
            for(size_t k=0; k<this->K_bytes; k++) {
                parity_buf ^= input[k] & this->parity_checks[i*8+j][k];
            }
            flag_buf = byte_parity(parity_buf);
            
            modify_bit(&output_buf, j, flag_buf);
        }
        out[i+par_ofst] = output_buf;
    }
    
    if(this->systype == systematic::FRONT) {
        std::memcpy(out, input, this->K_bytes);
    } else if(this->systype == systematic::BACK) {
        std::memcpy(&out[this->N_punct_bytes], input, this->K_bytes);
    }
}

uint8_t encoder::read_byte(FILE *fp, const char *descr) {
    uint8_t buf;
    
    if(fread(&buf, sizeof(uint8_t), 1, fp) != 1) {
        fprintf(stderr, "Unable to read byte for %s.\n");
        exit( EXIT_FAILURE );
    }
    
    return buf;
}

bool encoder::byte_parity(uint8_t byte) {
    uint8_t par = 0;
    
    for(uint8_t i=0; i<8; i++) {
        par ^= ( 0x01 & (byte>>i));
    }
    
    return (par > 0);
}

void encoder::modify_bit(uint8_t *byte, uint8_t pos, bool value) {
    uint8_t buf = (0x01 << (7-pos)); // one at pos
    
    if(value) {
        *byte |= buf; // Set to one, by OR
    } else {
        *byte &= (0xFF ^ buf); // Set to zero by AND (0xFF^buf gives all ones except at position pos
    }
}
