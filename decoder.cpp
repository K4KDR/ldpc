#include <ldpc/decoder.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>


void ldpc::decoder::parse_alist(const char* alist_file) {
    FILE* f = fopen(alist_file, "r");
    
    if(!f) {
        fprintf(stderr, "Cannot open file %s\n", alist_file);
        exit( EXIT_FAILURE );
    }
    
    uint64_t buf[2];
    
    // Read N M
    parse_numbers_from_file(buf, f, "dimensions", 2, false);
    this->N = buf[0];
    this->M = buf[1];
    
    // Read biggest_num_n biggest_num_m (ignored)
    parse_numbers_from_file(buf, f, "maximum elements", 2, false);
    
    // Read num_n
    this->nlist_num = new uint64_t[this->N];
    parse_numbers_from_file(this->nlist_num, f, "nlist count", this->N , false);
    
    // Read num_m
    this->mlist_num = new uint64_t[this->M];
    parse_numbers_from_file(this->mlist_num, f, "mlist count", this->M , false);
    
    // Read nlist
    this->nlist = new uint64_t*[this->N];
    for(size_t i=0; i<N ;i++) {
        this->nlist[i] = new uint64_t[this->nlist_num[i]];
        
        parse_numbers_from_file(this->nlist[i], f, "n-list", this->nlist_num[i], true);
    }

    // Read mlist
    this->mlist = new uint64_t*[this->M];
    for(size_t i=0; i<M ;i++) {
        this->mlist[i] = new uint64_t[this->mlist_num[i]];
        
        parse_numbers_from_file(this->mlist[i], f, "m-list", this->mlist_num[i], true);
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
    
    //// Alist read, transpose if necessary
    
    if(this->N > this->M) {
        // Swap N and M
        uint64_t tmp = this->N;
        this->N = this->M;
        this->M =tmp;
        
        // Swap nlist_num and mlist_num
        uint64_t *tmpp = this->nlist_num;
        this->nlist_num = this->mlist_num;
        this->mlist_num = tmpp;
        
        // Swap nlist and mlist
        uint64_t **tmppp = this->nlist;
        this->nlist = this->mlist;
        this->mlist = tmppp;
    }
    
    // Compute K
    this->K = M-N;
    
}


void ldpc::decoder::parse_numbers_from_file(uint64_t *ret, FILE *f, const char *line_descr, uint64_t num, bool ignore_zeros) {
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
    
    len = 0;
    errno = 0;
    //fprintf(stdout, "Parsing line: '%s'\n", line);
    for (unsigned long long i = strtoull(line, &end, 10); line != end; i = strtoull(line, &end, 10)) {
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
                ret[len++] = (uint64_t) i;
            }
        }
        errno = 0;
    }
    
    free(line_alloc);
    
    if(len != num) {
        fprintf(stderr, "%lu numbers read in %s line, but %lu were expected.\n", len, line_descr, num);
        exit( EXIT_FAILURE );
    }
    
    return;
}

uint64_t ldpc::decoder::get_num_input(void) const {
    return this->M - this->punctconf->num_punct;
}

uint64_t ldpc::decoder::get_num_output(void) const {
    return (this->systype == systematic::NONE) ? this->M : this->K;
}
        
ldpc::decoder::decoder(const char* alist_file, systematic::systematic_t systype, puncturing::conf_t *punctconf) {
    
    // Read in N, M, K, nlist, mlist, nlist_num, mlist_num
    this->parse_alist(alist_file);
    
    // Store systematics configuration
    this->systype = systype;
    
    // Store puncturing configuration
    this->punctconf = punctconf;
    
    // Initialize parity check nodes
    this->check_nodes = new check_node*[this->N];
    for(size_t i=0; i<this->N; i++) {
        this->check_nodes[i] = new check_node(this->nlist_num[i]);
    }
    
    // Initialize bit nodes
    this->bit_nodes = new bit_node*[this->M];
    for(size_t i=0; i<this->M; i++) {
        this->bit_nodes[i] = new bit_node(this->mlist_num[i]);
    }
    
}

ldpc::decoder::~decoder() {
    for(size_t i=0; i<this->N ;i++) {
        delete[] this->nlist[i];
    }
    delete[] this->nlist;
    delete[] this->nlist_num;
    
    for(size_t i=0; i<this->M ;i++) {
        delete[] this->mlist[i];
    }
    delete[] this->mlist;
    delete[] this->mlist_num;
    
}

uint64_t ldpc::decoder::get_syndrome_count(void) {
    uint64_t count=0;
    
    for(size_t i=0; i<this->N; i++) {
        count += (this->check_nodes[i]->isFullfilled()) ? 0 : 1;
    }
    
    return count;
}

bool ldpc::decoder::decode(softbit_t *out, const softbit_t *input, metadata_t *meta, const char *debugout) {
    uint64_t i, j;
    
    j=0;
    for(i=0; i<this->M; i++) {
        if(this->punctconf->is_punctured(i, this->M)) {
            this->bit_nodes[i]->reset(0.0f);
        } else {
            this->bit_nodes[i]->reset(input[j++]);
        }
    }
    
    FILE *debugf = NULL;
    if(debugout) {
        debugf = fopen(debugout,"w");
        
        if(!debugf) {
            fprintf(stderr, "Cannot open debug file %s\n", debugout);
            exit( EXIT_FAILURE );
        }
    }
    
    uint64_t syndrome_count;
    uint64_t bit_indx;
    uint64_t check_indx;
    uint64_t iteration_counter = 0;
    
    uint64_t syndrome_num_fulfilled=0;
    uint64_t syndrome_num_failed=0;
    uint64_t syndrome_num_undef=0;
    uint64_t syndrome_num_undef_punct=0;
    
    do {
        // Propagate new round
        for(i=0;i<this->N;i++) {
            this->check_nodes[i]->new_round();
        }
        for(i=0;i<this->M;i++) {
            this->bit_nodes[i]->new_round();
        }
        
        // propagate values from bit nodes to check nodes
        for(bit_indx=0; bit_indx<this->M; bit_indx++) {
            for(i=0; i<this->mlist_num[bit_indx]; i++) {
                check_indx = this->mlist[bit_indx][i]-1;
                this->check_nodes[check_indx]->set_bit_value(this->bit_nodes[bit_indx]->computeValForCheck(i, false));
            }
        }
        
        // propagate check values back to bit nodes
        for(check_indx=0; check_indx<this->N; check_indx++) {
            for(i=0; i<this->nlist_num[check_indx]; i++) {
                bit_indx = this->nlist[check_indx][i]-1;
                this->bit_nodes[bit_indx]->set_check_value(this->check_nodes[check_indx]->computeValForMessage(i));
            }
            
            if(iteration_counter==0) {
                bool contains_punct=false;
                for(i=0; i<this->nlist_num[check_indx]; i++) {
                    bit_indx = this->nlist[check_indx][i]-1;
                    contains_punct = this->punctconf->is_punctured(bit_indx, this->M) ? true : contains_punct;
                }
                if(this->check_nodes[check_indx]->syndrome_def) {
                    if(this->check_nodes[check_indx]->syndrome) {
                        syndrome_num_failed++;
                    } else {
                        syndrome_num_fulfilled++;
                    }
                    if(contains_punct) {
                        printf("Syndrome %u is defined although it contains punctured bits.\n", check_indx);
                    }
                } else {
                    syndrome_num_undef++;
                    syndrome_num_undef_punct += (contains_punct) ? 1 : 0;
                }
            }
        }
        
        if(iteration_counter==0) {
            printf("== Check results: ==\n");
            printf("  %4u checks SUCCEDED\n", syndrome_num_fulfilled);
            printf("  %4u checks FAILED\n", syndrome_num_failed);
            printf("  %4u checks are undefined due to puncturing\n", syndrome_num_undef_punct);
            printf("  %4u checks are undefined due to receiving undefined bits\n", syndrome_num_undef-syndrome_num_undef_punct);
            printf("====================\n");
        }
        
        // Compute number of unfulfilled syndromes
        syndrome_count = this->get_syndrome_count();
        
        //printf("  Decoding round gave %u syndrome errors.\n", syndrome_count);
        iteration_counter++;
        
        
        // Print probability of ones to debug file
        if(debugf) {
            for(bit_indx=0; bit_indx<this->M; bit_indx++) {
                fprintf(debugf, "%f ", ldpc::llr2prob(this->bit_nodes[bit_indx]->computeValForCheck(0, true)));
            }
            fprintf(debugf, "\n");
            printf("  decoding round %4u, %4u syndrome errors.\n", iteration_counter, syndrome_count);
        }
        
    } while(syndrome_count > 0 && iteration_counter < DECODER_MAX_ITERATIONS);
    
    if(debugf) {
        fclose(debugf);
    }
    
    uint64_t index_out_first;
    uint64_t index_out_last;
    switch(this->systype) {
        case systematic::NONE:
            // Output all M bits
            index_out_first = 0;
            index_out_last = this->M;
            break;
        case systematic::FRONT:
            // Output first K bits
            index_out_first = 0;
            index_out_last = this->K;
            break;
        case systematic::BACK:
            // Output last K bits
            index_out_first = this->M-this->K;
            index_out_last = this->M;
            break;
        default:
            fprintf(stderr, "State machine error.\n");
            exit( EXIT_FAILURE );
    }
    
    // Final iteration
    uint64_t ber_counter = 0;
    softbit_t tmp_bit;
    j=0;
    for(i=0; i<this->M; i++) {
        tmp_bit = this->bit_nodes[i]->computeValForCheck(0, true);
        
        if(i>=index_out_first && i<index_out_last) {
            out[j++] = tmp_bit;
        }
        
        ber_counter += (!this->punctconf->is_punctured(i,this->M) && my_abs(input[i])>1e-10 && input[i]*tmp_bit<0.0f) ? 1 : 0;
    }
    
    bool success = (syndrome_count==0 && iteration_counter<DECODER_MAX_ITERATIONS);
    if(meta) {
        meta->num_iterations = iteration_counter;
        meta->success = success;
        meta->num_corrected = ber_counter;
    }
    
    return success;
}

////
//////  Parity check node
////
ldpc::decoder::check_node::check_node(const uint64_t NUM_BITS) : NUM_BITS(NUM_BITS) {
    this->bit_values_tanh = new softbit_t[NUM_BITS];
    this->tmp_indx = 0;
}

ldpc::decoder::check_node::~check_node(void) {
    delete[] this->bit_values_tanh;
}

bool ldpc::decoder::check_node::isFullfilled(void) {
    // Syndrome must be fulfilled and zero (false)
    return this->syndrome_def && !(this->syndrome);
}

void ldpc::decoder::check_node::new_round(void) {
    if(this->tmp_indx != 0 && this->tmp_indx != this->NUM_BITS) {
        fprintf(stderr, "Resetting parity check node, although not all values are read out (%u/%u)",this->tmp_indx, this->NUM_BITS );
        exit( EXIT_FAILURE );
    }
    
    this->tmp_indx = 0;
    this->syndrome = false;
    this->syndrome_def = true;
}

void ldpc::decoder::check_node::set_bit_value(softbit_t val) {
    const bool debug=false;
    
    // compute tanh(LLR/2) and store it
    this->bit_values_tanh[this->tmp_indx] = tanh(val/2.0f);
    
    if(debug) {
        printf("Set value for parity check tanh(%12f / 2) = %12f\n", val, tanh(val/2.0f));
        printf(":: set value %f, syndrome %s%1u%s => ", val, this->syndrome_def ? "" : "(", this->syndrome ? 1 : 0, this->syndrome_def ? "" : ")");
    }
    
    // if bit is one (i.e. LLR < 0), flip syndrome
    this->syndrome ^= (val < 0.0f) ? true : false;
    
    // if bit is undefined (i.e. |LLR| < threshold), mark syndrome ans undefined
    this->syndrome_def = (my_abs(val) > 1e-5) ? this->syndrome_def : false;
    
    if(debug) {
        printf("%s%1u%s\n", this->syndrome_def ? "" : "(", this->syndrome ? 1 : 0, this->syndrome_def ? "" : ")");
    }
    
    // increase counter
    this->tmp_indx++;
}

ldpc::softbit_t ldpc::decoder::check_node::computeValForMessage(uint64_t indx) {
    const bool debug=false;
    
    softbit_t tmp_prod = 1.0f;
    
    if(debug)
        printf("  compute value for bit node %u\n", indx);
    
    for(size_t i=0; i<this->NUM_BITS; i++) {
        if(i==indx) {
            
            if(debug)
                printf("  ([%3u]) = tanh: %12f\n", i, this->bit_values_tanh[i]);
            
            continue;
        }
        
        if(debug)
            printf("   [%3u]  = tanh: %12f\n", i, this->bit_values_tanh[i]);
        
        tmp_prod *= this->bit_values_tanh[i];
    }
    
    if(debug)
        printf("  ======> %12f\n\n", log10( (1.0f+tmp_prod) / (1.0f-tmp_prod) ));
    
    return log10( (1.0f+tmp_prod) / (1.0f-tmp_prod) );
}


////
//////  Bit node
////
ldpc::decoder::bit_node::bit_node(const uint64_t NUM_CHECKS) : NUM_CHECKS(NUM_CHECKS){
    this->check_values = new softbit_t[NUM_CHECKS];
}

ldpc::decoder::bit_node::~bit_node(void){
    delete[] this->check_values;
}

void ldpc::decoder::bit_node::reset(softbit_t channel_val) {
    this->channel_value = channel_val;
    
    for(size_t i=0; i<this->NUM_CHECKS; i++) {
        this->check_values[i] = 0.0f;
    }
}

void ldpc::decoder::bit_node::new_round(void) {
    if(this->tmp_indx != 0 && this->tmp_indx != this->NUM_CHECKS) {
        fprintf(stderr, "Resetting bit node, although not all values are read out");
        exit( EXIT_FAILURE );
    }
    
    this->tmp_indx = 0;
}

void ldpc::decoder::bit_node::set_check_value(softbit_t val) {
    // store new check value
    this->check_values[this->tmp_indx] = val;
    
    // increase counter
    this->tmp_indx++;
}

ldpc::softbit_t ldpc::decoder::bit_node::computeValForCheck(uint64_t indx, bool final) {
    const bool debug=false;
    
    if(debug) {
        printf("  compute value for check node %u (%sfinal)\n", indx, final ? "" : "not ");
        printf("    values:\t channel : %12f\n", this->channel_value);
    
        for(size_t i=0; i<this->NUM_CHECKS; i++) {
            if(!final && i==indx)
                printf("          \t([% 5u]): %12f\n", i, this->check_values[i]);
            else
                printf("          \t [% 5u] : %12f\n", i, this->check_values[i]);
        }
    }
    
    llrsum_t sum;
    llrsum_reset(&sum);
    
    llrsum_add(&sum, this->channel_value);
    
    //softbit_t ret = this->channel_value;    
    for(size_t i=0; i<this->NUM_CHECKS; i++) {
        //printf("           \t [% 5u]: %f\n", i, this->check_values[i]);
        if(!final && i==indx) {
            continue;
        }
        //printf(":: %f +", ret);
        //ret = addllrs(ret, this->check_values[i]);
        llrsum_add(&sum, this->check_values[i]);
        //printf(" %f = %f\n", this->check_values[i], ret);
    }

    if(debug) {
        printf("    \t\t\t   -------------\n");
        printf("    \t\t\t   %12f\n\n", llrsum_get(&sum));
    }
    
    //return ret;
    return llrsum_get(&sum);
}

void ldpc::decoder::bit_node::llrsum_reset(llrsum_t *l) const {
    l->inf_count = 0;
    l->fin_sum = 0.0f;
}
        
void ldpc::decoder::bit_node::llrsum_add(llrsum_t *l, softbit_t val) const {
    if(isinf(val)) {
        l->inf_count += (val > 0) ? 1 : -1;
    } else {
        l->fin_sum += val;
    }
}
        
ldpc::softbit_t ldpc::decoder::bit_node::llrsum_get(llrsum_t *l) const {
    if(l->inf_count == 0) {
        // Finite
        return l->fin_sum;
    } else if(l->inf_count > 0) {
        // + Inf
        return 1.0f/0.0f;
    } else {
        // - Inf
        return -1.0f/0.0f;
    }
}

