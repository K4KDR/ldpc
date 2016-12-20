#include <ldpc/decoder.h>

void test01(void) {
    ldpc::puncturing::conf_t pconf(ldpc::puncturing::NONE, 0, NULL);
    ldpc::decoder::metadata_t meta;
    
    ldpc::decoder d("/home/v1tzl1/work/MOVE/LDPC/code/codes/martin_test.a", ldpc::systematic::FRONT, &pconf);
    
    ldpc::softbit_t buf_in[16];
    ldpc::softbit_t buf_out[16];
    
    //*
    // All zero input (0x00 => 0x0000)
    for(size_t i=0; i<16; i++) {
        buf_in[i] = ldpc::prob2llr(0.0f);
    }
    d.decode(buf_out, buf_in, &meta);
    printf("Decoding %s after %u iterations. %u bits corrected.\n", meta.success ? "SUCCESSFULL" : "FAILED", meta.num_iterations, meta.num_corrected);
    //*/
    
    //*
    // All ones input (0xFF => 0xFF00)
    for(size_t i=0; i<8; i++) {
        buf_in[i] = ldpc::prob2llr(1.0f);
    }
    for(size_t i=8; i<16; i++) {
        buf_in[i] = ldpc::prob2llr(0.0f);
    }
    d.decode(buf_out, buf_in, &meta);
    printf("Decoding %s after %u iterations. %u bits corrected.\n", meta.success ? "SUCCESSFULL" : "FAILED", meta.num_iterations, meta.num_corrected);
    //*/
    
    //*
    // All ones input (0xFF => 0xFF00) with bit error
    for(size_t i=0; i<8; i++) {
        buf_in[i] = ldpc::prob2llr(1.0f);
    }
    for(size_t i=8; i<16; i++) {
        buf_in[i] = ldpc::prob2llr(0.0f);
    }
    buf_in[5] = ldpc::prob2llr(0.0f);
    d.decode(buf_out, buf_in, &meta);
    printf("Decoding %s after %u iterations. %u bits corrected.\n", meta.success ? "SUCCESSFULL" : "FAILED", meta.num_iterations, meta.num_corrected);
    //*/
}

void test02(void) {
    // All ones
    // encoded message:
    // 1024 bits 0xFF (message)
    // 512  bits 0x00
    // 1024 bits 0xFF
    
    ldpc::puncturing::conf_t pconf(ldpc::puncturing::BACK, 512, NULL);
    //ldpc::puncturing::conf_t pconf(ldpc::puncturing::NONE, 0, NULL);
    
    ldpc::decoder::metadata_t meta;
    
    ldpc::decoder d("/home/v1tzl1/work/MOVE/LDPC/code/codes/AR4JA_r12_k1024n.a", ldpc::systematic::FRONT, &pconf);
    
    ldpc::softbit_t buf_in[2560];
    ldpc::softbit_t buf_out[1024];
    
    for(size_t i=0; i<2560; i++) {
        buf_in[i] = (i < 1024 || i >= 1024+512 ) ? -1e20f : 1e20f;
    }
    // Add some bit errors
    /*
    buf_in[50] = 0.0f;
    buf_in[100] = 1e20f;
    buf_in[150] = 1e20f;
    buf_in[200] = 1e20f;
    buf_in[250] = 1e20f;
    buf_in[300] = 1e20f;
    buf_in[350] = 1e20f;
    //*/
    
    d.decode(buf_out, buf_in, &meta, "/tmp/decoding.txt");
    printf("Decoding %s after %u iterations. %u bits corrected.\n", meta.success ? "SUCCESSFULL" : "FAILED", meta.num_iterations, meta.num_corrected);
}

int main(void) {
    
    //test01();
    test02();
    
    printf("Finished.\n");
}
