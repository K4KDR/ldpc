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

bool ldpc::decoder::get_syndrome(const uint64_t check_indx, bool *defined) const {
    softbit_t tmp_bit;
    bool s_i = false;
    
#if LDPC_DO_SANITY_CHECKS
    if(check_indx >= this->N) {
        fprintf(stderr, "ERROR: Check index too large in ldpc::decoder::get_syndrome(%lu)\n", check_indx);
        exit( EXIT_FAILURE );
    }
#endif

    for(uint64_t j=0; j<this->nlist_num[check_indx]; j++) {
        
        tmp_bit = this->bit_nodes[this->nlist[check_indx][j]-1]->get_buffered_final_value();
        
        if(my_abs(tmp_bit) < DECODER_MIN_LLR_MAG) {
            // bit undefined, set syndrome to false
            if(defined) {
                *defined = false;
            }
            return true;
        } else if(tmp_bit < 0.0f) {
            s_i = !s_i;
        }
    }
    
    if(defined) {
        *defined = true;
    }
    return s_i;
}

uint64_t ldpc::decoder::get_syndrome_count(void) const {
    uint64_t count=0;
    for(size_t i=0; i<this->N; i++) {
        // do not trust check nodes syndrome, because they do not contain the final bit falue estmates
        //count += (this->check_nodes[i]->isFullfilled()) ? 0 : 1;
        
        count += (this->get_syndrome(i)) ? 1 : 0;
    }
    
    return count;
}

double ldpc::decoder::get_awrm(void) const {
    double ret=0.0;
    double e_i;
    double s_j;
    double abs_yi;
    double w_ij_tmp;
    double w_ij = 1.0/0.0;
    uint64_t j_indx, j;
    uint64_t k_indx, k;
    
    for(uint64_t i=0; i<this->M; i++) {
        //printf("Computing AWRM for bit %lu\n", i);
        abs_yi = my_abs(tanh(this->bit_nodes[i]->channel_value));
        //printf("  |y_i| = %lf\n", abs_yi);
        
        e_i = 0.0;
        for(j_indx=0; j_indx<this->mlist_num[i]; j_indx++) {
            j=this->mlist[i][j_indx]-1;
            s_j = (this->get_syndrome(j)) ? 0.0 : 1.0;
            
            w_ij = 1.0/0.0; // +inf
            for(k_indx=0; k_indx<this->nlist_num[j]; k_indx++) {
                k = this->nlist[j][k_indx]-1;
                if(k==i) {
                    continue;
                }
                w_ij_tmp = my_abs(tanh(this->bit_nodes[k]->channel_value));
                w_ij = (w_ij <= w_ij_tmp) ? w_ij : w_ij_tmp;
            }
            //printf("  j=%4lu: s_j=%3lf w_ij=%12lf, ()=%12lf\n", j, s_j, w_ij, (2.0*s_j-1.0)*w_ij);
            e_i += (2.0*s_j-1.0)*w_ij;
        }
        //printf("sum(.) = %12lf, e_{%4lu} = %12lf\n", e_i, i, e_i-abs_yi);
        e_i -= abs_yi;
        
        ret += e_i;
    }
    
    return ret/this->M;
}

ldpc::softbit_t ldpc::decoder::llrdiff(const ldpc::softbit_t a, const ldpc::softbit_t b) const {
    if(isinf(a) && isinf(b)) {
        // Both infinite
        
        if(a*b>0.0f) {
            // both the same, difference is zero
            return 0.0f;
        } else {
            // If a is -inf and b=+Inf => a-b=-Inf=a, else a-b=+Inf=b
            return (a<0) ? a : b;
        }
    } else if(!isinf(a) && !isinf(b)) {
        // Both values are finite, calculate difference
        return a-b;
    } else if(isinf(a)) {
        // a is +/- inf, result is still a
        return a;        
    } else if(isinf(b)) {
        // b is +/- inf, result is -b;
        return -b;
    } else {
        fprintf(stderr, "State machine error in llrdiff\n");
        exit( EXIT_FAILURE );
    }
}

void ldpc::decoder::debug_check(const uint64_t check_indx) {
    bool tmp_syn_def;
    bool tmp_syn = this->get_syndrome(check_indx, &tmp_syn_def);
    printf("Debug check node %lu\n", check_indx);
    for(uint64_t i=0; i<this->nlist_num[check_indx]; i++) {
        uint64_t j=this->nlist[check_indx][i]-1;
        
	bool flag = false;
        uint64_t check_number=0;
        for(size_t k=0; k<this->mlist_num[j]; k++) {
            if(this->mlist[j][k]-1 == check_indx) {
                check_number = k;
		flag = true;
                break;
            }
        }
        if(flag) {
	    printf("  no check number found.\n");
	} else {
            printf("  connected to bit %lu: %12f final: %12f\n", j, this->bit_nodes[j]->computeValForCheck(check_number, false), this->bit_nodes[j]->get_buffered_final_value());
	}
    }
    printf("  syndrome: %s%1u%s\n", tmp_syn_def ? " " : "(", tmp_syn ? 1u : 0u, tmp_syn_def ? " " : ")");
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
    
    uint64_t awrm_counter = 0;
    double awrm_tmp = nan("");
#if DECODER_MAX_AWRM_ITERATIONS<DECODER_MAX_ITERATIONS
    double awrm_min = 0.0;
#endif
    
    softbit_t bits_last_it[this->M];
    std::memcpy(bits_last_it, input, this->M*sizeof(softbit_t));
    softbit_t tmp_softbit, delta_bits_sum;
    
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
            // Compute final estimate
            this->bit_nodes[bit_indx]->computeValForCheck(0, true);
        }
        
        // Evaluate syndromes at first iteration
        if(iteration_counter==0) {
            bool contains_punct=false;
            bool tmp_syndrome;
            bool tmp_syndrome_def;
            
            for(check_indx=0; check_indx<this->N; check_indx++) {
                for(i=0; i<this->nlist_num[check_indx]; i++) {
                    bit_indx = this->nlist[check_indx][i]-1;
                    contains_punct = this->punctconf->is_punctured(bit_indx, this->M) ? true : contains_punct;
                }
                
                if(debugout) {
                    tmp_syndrome = this->get_syndrome(check_indx, &tmp_syndrome_def);
                    if(tmp_syndrome_def) {
                        if(tmp_syndrome) {
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
	    
            if(debugout) {
                    printf("== Check results: ==\n");
                    printf("  %4u checks SUCCEDED\n", syndrome_num_fulfilled);
                    printf("  %4u checks FAILED\n", syndrome_num_failed);
                    printf("  %4u checks are undefined due to puncturing\n", syndrome_num_undef_punct);
                    printf("  %4u checks are undefined due to receiving undefined bits\n", syndrome_num_undef-syndrome_num_undef_punct);
                    printf("====================\n");
            }
        }

        // propagate check values back to bit nodes
        for(check_indx=0; check_indx<this->N; check_indx++) {
            for(i=0; i<this->nlist_num[check_indx]; i++) {
                bit_indx = this->nlist[check_indx][i]-1;
                this->bit_nodes[bit_indx]->set_check_value(this->check_nodes[check_indx]->computeValForMessage(i));
            }
        }
        
        // Compute number of unfulfilled syndromes
        syndrome_count = this->get_syndrome_count();
        
        // AWRM stopping criterion
#if DECODER_MAX_AWRM_ITERATIONS<DECODER_MAX_ITERATIONS
        awrm_tmp = this->get_awrm();
        if(awrm_tmp < awrm_min) {
            awrm_min = awrm_tmp;
            awrm_counter = 0;
        } else {
            awrm_counter++;
        }
#endif

        //printf("  Decoding round gave %u syndrome errors.\n", syndrome_count);
        iteration_counter++;
        
        // See if bits have settled
        //printf("Compute LLR delta sum\n");
        delta_bits_sum = 0.0f;
        for(bit_indx=0; bit_indx<this->M; bit_indx++) {
            tmp_softbit = this->bit_nodes[bit_indx]->get_buffered_final_value();
            //printf(" %4lu: old %12f => new %12f => |diff| %12f\n", bit_indx, bits_last_it[bit_indx], tmp_softbit, my_abs(llrdiff(bits_last_it[bit_indx],tmp_softbit)));
            //printf("%12.4f + %12.4f = %12.4f\n", delta_bits_sum, my_abs(llrdiff(bits_last_it[bit_indx],tmp_softbit)), delta_bits_sum+my_abs(llrdiff(bits_last_it[bit_indx],tmp_softbit)) );
            delta_bits_sum += my_abs(llrdiff(bits_last_it[bit_indx],tmp_softbit));
            bits_last_it[bit_indx] = tmp_softbit;
        }
        //printf("%12.4f / %12.4f = %12.4f\n", delta_bits_sum, (softbit_t)this->M, delta_bits_sum / (softbit_t)this->M );
        delta_bits_sum /= (softbit_t)this->M;
        
        
        // Print probability of ones to debug file
        if(debugf) {
            for(bit_indx=0; bit_indx<this->M; bit_indx++) {
                fprintf(debugf, "%f ", ldpc::llr2prob(this->bit_nodes[bit_indx]->get_buffered_final_value()));
            }
            fprintf(debugf, "\n");
            printf("  decoding round %4u/%4u, %4u syndrome errors, AWRM = %12lf (%4u/%4u), delta LLRs=%12le.\n", iteration_counter, DECODER_MAX_ITERATIONS, syndrome_count, awrm_tmp, awrm_counter, DECODER_MAX_AWRM_ITERATIONS, delta_bits_sum);
            
            /*
            if(syndrome_count > 0) {
                if(syndrome_count < 5) {
                    for(size_t i=0; i<this->N; i++) {
                        if(this->get_syndrome(i)) {
                            this->debug_check(i);
                        }
                    }
                } else {
                    printf("  unfulfilled syndromes are: ");
                    for(size_t i=0; i<this->N; i++) {
                        if(this->get_syndrome(i)) {
                            printf("%lu ", i);
                        }
                    }
                    printf("\n");
                }
            }
            //*/
            
            /*
            const uint64_t debug_check = 227;
            bool tmp_syn_def;
            bool tmp_syn = this->get_syndrome(debug_check, &tmp_syn_def);
            printf("Tracking check node %lu\n", debug_check);
            for(uint64_t i=0; i<this->nlist_num[debug_check]; i++) {
                uint64_t j=this->nlist[debug_check][i]-1;
                
                uint64_t check_number;
                for(size_t k=0; k<this->mlist_num[j]; k++) {
                    if(this->mlist[j][k]-1 == debug_check) {
                        check_number = k;
                        break;
                    }
                }
                
                printf("  connected to bit %lu: %12f final: %12f\n", j, this->bit_nodes[j]->computeValForCheck(check_number, false), this->bit_nodes[j]->get_buffered_final_value());
            }
            printf("  syndrome: %s%1u%s\n", tmp_syn_def ? " " : "(", tmp_syn ? 1u : 0u, tmp_syn_def ? " " : ")");
            //*/
            
        }
        
    } while(syndrome_count > 0 && iteration_counter < DECODER_MAX_ITERATIONS && awrm_counter < DECODER_MAX_AWRM_ITERATIONS && (isnan(delta_bits_sum) || delta_bits_sum>0.0f));
    
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
        tmp_bit = this->bit_nodes[i]->get_buffered_final_value();
        
        if(i>=index_out_first && i<index_out_last) {
            out[j++] = tmp_bit;
        }
        
        ber_counter += (!this->punctconf->is_punctured(i,this->M) && input[i]*tmp_bit<0.0f) ? 1 : 0;
    }
    
    uint8_t fail_flags = NONE;
    fail_flags |= (iteration_counter>=DECODER_MAX_ITERATIONS) ? MAX_ITERATIONS     : NONE;
    fail_flags |= (awrm_counter>=DECODER_MAX_AWRM_ITERATIONS) ? AWRM_STOP          : NONE;
    fail_flags |= (delta_bits_sum==0.0f)                      ? NO_SOFTBITS_CHANGE : NONE;
    
    bool success = (syndrome_count==0 && fail_flags==NONE);
    
    if(meta) {
        meta->num_iterations = iteration_counter;
        meta->num_iterations_max = DECODER_MAX_ITERATIONS;
        meta->success = success;
        meta->failure_flags = fail_flags;
        meta->syndrome_count = syndrome_count;
        meta->num_corrected = ber_counter;
        meta->num_bits_total = this->M;
        meta->num_guesses = 0;
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

void ldpc::decoder::check_node::new_round(void) {
#if LDPC_DO_SANITY_CHECKS
    if(this->tmp_indx != 0 && this->tmp_indx != this->NUM_BITS) {
        fprintf(stderr, "Resetting parity check node, although not all values are read out (%u/%u)",this->tmp_indx, this->NUM_BITS );
        exit( EXIT_FAILURE );
    }
#endif

    this->tmp_indx = 0;
}

void ldpc::decoder::check_node::set_bit_value(softbit_t val) {
    // compute tanh(LLR/2) and store it
    this->bit_values_tanh[this->tmp_indx] = tanh(val/2.0f);
    
    // increase counter
    this->tmp_indx++;
}

ldpc::softbit_t ldpc::decoder::check_node::computeValForMessage(uint64_t indx) const {
    softbit_t tmp_prod = 1.0f;
    
    for(size_t i=0; i<this->NUM_BITS; i++) {
        if(i==indx) {
            continue;
        }
        
        tmp_prod *= this->bit_values_tanh[i];
    }
    
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
    
    this->final_value = nan(""); // mark value as unset
}

void ldpc::decoder::bit_node::new_round(void) {
#if LDPC_DO_SANITY_CHECKS
    if(this->tmp_indx != 0 && this->tmp_indx != this->NUM_CHECKS) {
        fprintf(stderr, "Resetting bit node, although not all values are read out");
        exit( EXIT_FAILURE );
    }
#endif

    this->final_value = nan(""); // mark value as unset
    this->tmp_indx = 0;
}

void ldpc::decoder::bit_node::set_check_value(softbit_t val) {
    // store new check value
    this->check_values[this->tmp_indx] = val;
    
    // increase counter
    this->tmp_indx++;
}

ldpc::softbit_t ldpc::decoder::bit_node::get_buffered_final_value(void) const {
#if LDPC_DO_SANITY_CHECKS
    if(isnan(this->final_value)) {
        fprintf(stderr, "ERROR: Access to bits final estimate, before it is computed.\n");
        exit( EXIT_FAILURE );
    }
#endif

    return this->final_value;
}

ldpc::softbit_t ldpc::decoder::bit_node::computeValForCheck(uint64_t indx, bool final) {
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

    // Store final value in class
    if(final) {
        this->final_value = llrsum_get(&sum);
    }
    
    //return ret;
    return llrsum_get(&sum);
}

void ldpc::decoder::bit_node::llrsum_reset(llrsum_t *l) const {
    l->inf_count = 0;
    l->fin_sum = 0.0f;
}
        
void ldpc::decoder::bit_node::llrsum_add(llrsum_t *l, softbit_t val) const {
    l->inf_count += isinf(val) ? ((val > 0) ? 1 : -1) : 0;
    l->fin_sum += isinf(val) ? 0 : val;
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


//
//// Guess tree
//

ldpc::decoder::guess_tree::guess_tree(guess_tree *parent, uint64_t guess_pos, bool pref_pos, uint64_t max_children) {
    this->parent = parent;
    this->level = (this->parent) ? this->parent->level+1 : 0;
    
    this->guess_pos = guess_pos;
    this->pref_pos = pref_pos;
    
    this->num_children = 0;
    
    if(max_children > 0) {
        this->children = new guess_tree*[max_children];
        
        if(!this->children) {
            fprintf(stderr, "ERROR: allocation for children failed in ldpc::decoder::guess_tree\n");
            exit( EXIT_FAILURE );
        }
        this->num_children_alloc = max_children;
    } else {
        this->num_children_alloc = 0;
    }
    
    this->results[0] = PENDING;
    this->results[1] = PENDING;
    
    this->traverse_counter = 0;
}

ldpc::decoder::guess_tree::~guess_tree(void) {
    for(size_t i=0; i<this->num_children; i++) {
        delete this->children[i];
    }
    delete[] this->children;
}

void ldpc::decoder::guess_tree::add_child(uint64_t guess_pos, bool pref_pos, uint64_t max_children) {
    if(this->num_children >= this->num_children_alloc) {
        guess_tree **new_children = new guess_tree*[this->num_children+1];
        std::memcpy(new_children, this->children, this->num_children*sizeof(guess_tree*));
        
        this->num_children_alloc++;
        
        delete[] this->children;
        this->children = new_children;
    }
    
    this->children[this->num_children] = new guess_tree(this, guess_pos, pref_pos, max_children);
    this->num_children++;
}

ldpc::decoder::guess_tree *ldpc::decoder::guess_tree::traverse(void) {
    if(this->traverse_counter == 0) {
        this->traverse_counter++;
        return this;
    } else if(this->traverse_counter <= this->num_children) {
        const uint64_t child_indx = this->traverse_counter - 1;
        this->traverse_counter++;
        return this->children[child_indx];
    } else if(this->traverse_counter == this->num_children+1) {
        this->traverse_counter++;
        return this->parent;
    } else {
        fprintf(stderr, "ERROR: Sate machine error while traversing ldpc::decoder::guess_tree\n");
        exit( EXIT_FAILURE );
    }
}

void ldpc::decoder::guess_tree::reset_traverse(void) {
    this->traverse_counter = 0;
}

std::string ldpc::decoder::guess_tree::get_str(void) {
    const uint64_t NUM_LEVELS = this->level + 1;
    char **buffer = new char*[NUM_LEVELS];
    guess_tree *tmp=this;
    for(size_t l=0; l<NUM_LEVELS; l++) {
        buffer[l] = new char[7];
        sprintf(buffer[l], "%05lu%1s", tmp->guess_pos, tmp->pref_pos ? "p" : "n");
        tmp = tmp->parent;
    }
    
    std::string tmp_str;
    for(size_t l=0; l<NUM_LEVELS; l++) {
        if(l>0) {
            tmp_str.append(1u, '-');
        }
        tmp_str.append(buffer[NUM_LEVELS-l-1]);
    }
    
    for(size_t l=0; l<NUM_LEVELS; l++) {
        delete[] buffer[l];
    }
    delete[] buffer;
    
    return tmp_str;
}
