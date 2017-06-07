#include <ldpc/ldpc.h>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

ldpc::puncturing::conf_t::conf_t(void) : type(NONE), num_punct(0) {}
            
ldpc::puncturing::conf_t::conf_t(const puncturing_t type, const uint64_t num_punct, uint64_t *punct_pos)
    : type(type), num_punct(num_punct) {
        
    if(this->type == CUSTOM) {
        this->punct_pos = new uint64_t[num_punct];
        std::memcpy(this->punct_pos, punct_pos, num_punct*sizeof(uint64_t));
    }
}
        
ldpc::puncturing::conf_t::conf_t(const conf_t &cpy)
    : type(cpy.type), num_punct(cpy.num_punct) {
        
    if(this->type == CUSTOM) {
        this->punct_pos = new uint64_t[this->num_punct];
        std::memcpy(this->punct_pos, cpy.punct_pos, this->num_punct*sizeof(uint64_t));
    }
}
        
ldpc::puncturing::conf_t::~conf_t(void) {
    if(this->type == CUSTOM) {
        delete[] this->punct_pos;
    }
}

bool ldpc::puncturing::conf_t::is_punctured(uint64_t indx, const uint64_t M) const {
    bool ret = false;;
    
    if(this->type == puncturing::NONE) {
        ret = false;
    } else if(this->type == puncturing::FRONT) {
        ret = (indx < this->num_punct);
    } else if(this->type == puncturing::BACK) {
        ret = (indx >= M-this->num_punct);
    } else if(this->type == puncturing::CUSTOM) {
        for(size_t i=0; i<this->num_punct; i++) {
            if(indx == this->punct_pos[i]) {
                ret = true;
                break;
            }
        }
    } else {
        fprintf(stderr, "State machine error.\n");
        exit( EXIT_FAILURE );
    }
    
    //printf("bit %4u of %4u %s punctured.\n", indx, M, ret ? "IS" : "is NOT");
    return ret;
}

ldpc::softbit_t ldpc::prob2llr(const softbit_t prob_one) {
    return log10((1.0f-prob_one)/prob_one);
}
    
ldpc::softbit_t ldpc::llr2prob(const softbit_t llr) {
    return 1.0f/(pow(10.0f, llr) + 1.0f);
}

ldpc::softbit_t ldpc::addllrs(const softbit_t val1, softbit_t val2) {
    if(isinf(val1) && isinf(val2)) {
        // both values are infinite, if the have the same sign, set value to infinite, otherwise set to zero
        if( (val1<0.0f && val2<0.0f) || (val1>0.0f && val2>0.0f) ) {
            return val1;
        } else {
            return 0.0f;
        }
    } else if(isinf(val1)) {
        // val1 is not finite (but val2 is, so val1 dominates)
        return val1;
    } else if(isinf(val2)) {
        // val2 is not finite (but val1 is, so val2 dominates)
        return val2;
    } else {
        // No value is infinite, just add them
        return val1+val2;
    }
}

//inline ldpc::softbit_t ldpc::my_abs(const softbit_t v) {
//    return (v>=0.0f) ? v : -v;
    /*if(v >= 0.0f) {
        return v;
    } else {
        return -v;
    }*/
//}
