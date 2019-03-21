#ifndef __LIBLDPC_DECODER_H__DEFINED__
#define __LIBLDPC_DECODER_H__DEFINED__

#include <stdint.h>
#include <stdio.h>
#include <ldpc/ldpc.h>
#include <string>

/** Maximum number of decoding iterations before decoding failure is declared */
#define DECODER_MAX_ITERATIONS 1000

/** Maximum number of iterations until the AWRM has to reach a new minimum or the codeword is considered undecodable
 * 
 * Set to DECODER_MAX_ITERATIONS or larger to disable early return due to AWRM heuristic
 */
#define DECODER_MAX_AWRM_ITERATIONS 25

/** Minimum LLR magnitude to consider bit defined
 * 
 * LLRs with lower abs value are considered undefined and will always generate a positive syndrome. Consequently the output LLRs of the decoder will all have a magnitude greater than this value if decoding is succesfull.
 */
#define DECODER_MIN_LLR_MAG 0.000001f

namespace ldpc {
    
    class decoder {
    private:
        class check_node;
        class bit_node;
        class guess_tree;

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
        
        enum fail_t { NONE=0x00, MAX_ITERATIONS=(0x01<<0), AWRM_STOP=(0x01<<1), NO_SOFTBITS_CHANGE=(0x01<<2) };
        
        struct metadata_t {
            /** Number of decoding iterations */
            uint64_t num_iterations;
            
            /** Whether or not decoding was successful */
            bool success;
            
            /** Flags to indicate failure reason(s) from the fail_t enum */
            uint8_t failure_flags;
            
            /** Number of violated syndrome elements */
            uint64_t syndrome_count;
            
            /** Number of bits that have been corrected
             * 
             * Bit is counted as corrected, if it was part of the input (i.e. not punctured), the
             * bit had an information assigned (i.e. LLR != 0.0) and the decoded softbit has a
             * different sign than the original input (i.e. decoder decided the original bit was wrong).
             */
            uint64_t num_corrected;
            
            /** Number of bit guesses performed by the decoder */
            uint64_t num_guesses;
            
            /** Maximum number of iterations */
            uint64_t num_iterations_max;
            
            /** Number of total bits
             * 
             * used to determin BER from \ref num_corrected
             */
            uint64_t num_bits_total;
        };
        
        uint64_t get_num_input(void) const;
        uint64_t get_num_output(void) const;
        
        bool decode(softbit_t *out, const softbit_t *input, metadata_t *metadata=NULL, const char *debugout=NULL); // decode K bits from M inputs
        bool decode_guess(softbit_t *out, const softbit_t *input, metadata_t *metadata=NULL, const char *debugprefix=NULL); // decode K bits from M inputs
        
    private:
        void parse_alist(const char* alist_file);
        void parse_numbers_from_file(uint64_t *ret, FILE *f, const char *line_descr, uint64_t num, bool ignore_zeros);
        bool get_syndrome(const uint64_t check_indx, bool *defined=NULL) const;
        uint64_t get_syndrome_count(void) const;
        double get_awrm(void) const;
        softbit_t llrdiff(const softbit_t a, softbit_t b) const;
        void debug_check(const uint64_t check_indx);
        
    };
    
    class decoder::guess_tree {
    public:
        guess_tree *parent;
        uint64_t level;
        
        guess_tree **children;
        uint64_t num_children;
        uint64_t num_children_alloc;
        
        
        uint64_t guess_pos;
        bool pref_pos;
        enum result_t { PENDING=0, WORSE=1, BETTER=2, SUCCESS=3 };
        
        result_t results[2];
        
        uint64_t traverse_counter;
        
        guess_tree(guess_tree *parent, uint64_t guess_pos, bool pref_pos, uint64_t max_children=0);
        ~guess_tree(void);
        
        void add_child(uint64_t guess_pos, bool pref_pos, uint64_t max_children=0);
        guess_tree* traverse(void);
        void reset_traverse(void);
        std::string get_str(void);
    };
    
    class decoder::check_node {
    public:
        const uint64_t NUM_BITS;
        softbit_t *bit_values_tanh;
        uint64_t tmp_indx;
        
        check_node(const uint64_t NUM_BITS);
        ~check_node(void);
        
        softbit_t computeValForMessage(uint64_t indx) const;
        void new_round(void); // Reset node to new round of value inputs
        void set_bit_value(softbit_t val); // Set a new bit value at tmp_indx
    };
    
    class decoder::bit_node {
    public:
        const uint64_t NUM_CHECKS;
        softbit_t *check_values;
        softbit_t channel_value;
        softbit_t final_value;
        uint64_t tmp_indx;
        
        bit_node(const uint64_t NUM_CHECKS);
        ~bit_node(void);
        
        void reset(softbit_t channel_val); // Reset decoder to decode new message
        void new_round(void); // Reset node to new round of value inputs
        void set_check_value(softbit_t val); // Set a new check value at tmp_indx
        softbit_t get_buffered_final_value(void) const; // Return buffered final estimate (include all check-nodes. Throws an error is not computed yet.
        
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
