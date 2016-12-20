#ifndef __LIBLDPC_DECODER_H__DEFINED__
#define __LIBLDPC_DECODER_H__DEFINED__

#include <stdint.h>
#include <stdio.h>
#include <ldpc/ldpc.h>

#define DECODER_MAX_ITERATIONS 50

namespace ldpc {
    
    class decoder {
    private:
        class check_node;
        class bit_node;

        /** Number of parity checks without puncturing */
        uint64_t N;
        
        /** Number of total bits without puncturing */
        uint64_t M;
        
        /** Number of information bits */
        uint64_t K;
        
        uint64_t *nlist_num;
        uint64_t *mlist_num;
        uint64_t **nlist;
        uint64_t **mlist;
        
        systematic::systematic_t systype;
        puncturing::conf_t *punctconf;

        check_node **check_nodes;
        bit_node **bit_nodes;
        
    public:
        decoder(const char* alist_file, systematic::systematic_t systype, puncturing::conf_t *punctconf);
        ~decoder();
        
        struct metadata_t {
            /** Number of decoding iterations */
            uint64_t num_iterations;
            
            /** Whether or not decoding was successful */
            bool success;
            
            /** Number of bits that have been corrected
             * 
             * Bit is counted as corrected, if it was part of the input (i.e. not punctured), the
             * bit had an information assigned (i.e. LLR != 0.0) and the decoded softbit has a
             * different sign than the original input (i.e. decoder decided the original bit was wrong).
             */
            uint64_t num_corrected;
        };
        
        uint64_t get_num_input(void) const;
        uint64_t get_num_output(void) const;
        
        bool decode(softbit_t *out, const softbit_t *input, metadata_t *metadata=NULL, const char *debugout=NULL); // decode M bits from M inputs
        
    private:
        void parse_alist(const char* alist_file);
        void parse_numbers_from_file(uint64_t *ret, FILE *f, const char *line_descr, uint64_t num, bool ignore_zeros);
        uint64_t get_syndrome_count(void);
    };
    
    class decoder::check_node {
    public:
        const uint64_t NUM_BITS;
        softbit_t *bit_values_tanh;
        uint64_t tmp_indx;
        bool syndrome_def;
        bool syndrome; /** Syndrome value, true = violation, false = fulfilment */
        
        check_node(const uint64_t NUM_BITS);
        ~check_node(void);
        
        bool isFullfilled(void);
        softbit_t computeValForMessage(uint64_t indx);
        void new_round(void); // Reset node to new round of value inputs
        void set_bit_value(softbit_t val); // Set a new bit value at tmp_indx
    };
    
    class decoder::bit_node {
    public:
        const uint64_t NUM_CHECKS;
        softbit_t *check_values;
        softbit_t channel_value;
        uint64_t tmp_indx;
        
        bit_node(const uint64_t NUM_CHECKS);
        ~bit_node(void);
        
        void reset(softbit_t channel_val); // Reset decoder to decode new message
        void new_round(void); // Reset node to new round of value inputs
        void set_check_value(softbit_t val); // Set a new check value at tmp_indx
        
        softbit_t computeValForCheck(uint64_t indx, bool final); // Compute probability (LLR) that bit should be one, based on all check node informations, but the one indicated by indx. In the final iteration all check messages are considered and indx is ignored
        
    private:
        struct llrsum_t {
            /** Counter for infinite terms.
             * 
             * If negative there are more -inf terms in the sum than +inf. If zero, there are either
             * equal ammounts of positive and negative infinity terms in the sum, or none at all.
             */
            int inf_count;
            
            /** Sum of finite sum elements */
            softbit_t fin_sum;
        };
        
        /** Set LLR sum to zero */
        void llrsum_reset(llrsum_t *l) const;
        
        /** Add softbit to sum */
        void llrsum_add(llrsum_t *l, softbit_t val) const;
        
        /** Return sum as softbit */
        softbit_t llrsum_get(llrsum_t *l) const;
        
        
        
    };
    
}

#endif /* __LIBLDPC_DECODER_H__DEFINED__ */
