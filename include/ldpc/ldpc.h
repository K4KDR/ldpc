#ifndef __LIBLDPC_LDPC_H__DEFINED__
#define __LIBLDPC_LDPC_H__DEFINED__

#include <cstring>
#include <math.h>
#include <stdint.h>

namespace ldpc {
    
    typedef float softbit_t;
    
    namespace systematic {
        enum systematic_t { NONE=0, FRONT=1, BACK=2 }; /** Whether the a copy of the original message should be added at the front, back, ot not at all */
    }
    
    namespace puncturing {
        enum puncturing_t { NONE=0, FRONT=1, BACK=2, CUSTOM=3 };
        
        class conf_t {
        public:
            const puncturing_t type;
            const uint64_t num_punct;
            uint64_t *punct_pos;
        
            // Default constructor, NO puncturing
            conf_t(void);
            
            conf_t(const puncturing_t type, const uint64_t num_punct, uint64_t *punct_pos);
            
            conf_t(const conf_t &cpy);
            
            ~conf_t(void);
            
            bool is_punctured(uint64_t indx, const uint64_t M) const;
        };
    }
    
    softbit_t prob2llr(const softbit_t prob_one);
    
    softbit_t llr2prob(const softbit_t llr);
    
    softbit_t addllrs(const softbit_t val1, softbit_t val2);
    
    softbit_t my_abs(const softbit_t v);
}

#endif /* __LIBLDPC_LDPC_H__DEFINED__ */
