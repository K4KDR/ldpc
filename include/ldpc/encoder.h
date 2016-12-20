#ifndef __LIBLDPC_ENCODER_H__DEFINED__
#define __LIBLDPC_ENCODER_H__DEFINED__

#include <stdint.h>
#include <stdio.h>
#include <ldpc/ldpc.h>

namespace ldpc {
    class encoder {
    private:
        uint64_t N;
        uint64_t K;
        uint64_t N_punct;
        uint64_t N_punct_bytes;
        uint64_t K_bytes;
        
        uint8_t **parity_checks;
        systematic::systematic_t systype;
        
    public:
        encoder(const char* generator_file, systematic::systematic_t systype, puncturing::conf_t *punctconf);
        ~encoder();
        
        uint64_t get_num_input(void) const;
        uint64_t get_num_output(void) const;
        
        void encode(uint8_t *out, const uint8_t *input);
        
    private:
        static uint8_t read_byte(FILE *fp, const char *descr);
        static bool byte_parity(uint8_t byte);
        static void modify_bit(uint8_t *byte, uint8_t pos, bool value);
    };
}

#endif /* __LIBLDPC_ENCODER_H__DEFINED__ */
