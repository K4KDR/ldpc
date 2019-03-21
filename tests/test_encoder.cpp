#include <ldpc/encoder.h>

void test01(void) {
    ldpc::puncturing::conf_t *punctconf = new ldpc::puncturing::conf_t(ldpc::puncturing::NONE, 0, NULL);
    ldpc::encoder *e = new ldpc::encoder("/home/v1tzl1/work/MOVE/LDPC/code/codes/gmartin_test.gen", ldpc::systematic::FRONT, punctconf);
    
    uint8_t buf_in[1];
    uint8_t buf_out[2];
    
    
    buf_in[0] = 0x00;
    e->encode(buf_out, buf_in);
    printf("%02X => %02X%02X\n", buf_in[0], buf_out[0], buf_out[1]);
    
    buf_in[0] = 0xFF;
    e->encode(buf_out, buf_in);
    printf("%02X => %02X%02X\n", buf_in[0], buf_out[0], buf_out[1]);
    
    buf_in[0] = 0x03;
    e->encode(buf_out, buf_in);
    printf("%02X => %02X%02X\n", buf_in[0], buf_out[0], buf_out[1]);
    
    buf_in[0] = 0xAA;
    e->encode(buf_out, buf_in);
    printf("%02X => %02X%02X\n", buf_in[0], buf_out[0], buf_out[1]);
    
    delete punctconf;
    delete e;
}

void test02(void) {
    ldpc::puncturing::conf_t punctconf(ldpc::puncturing::NONE, 0, NULL);
    ldpc::encoder e("/home/v1tzl1/work/MOVE/LDPC/code/codes/gAR4JA_r12_k1024n.gen", ldpc::systematic::FRONT, &punctconf);
    
    const uint64_t K_bytes = e.get_num_input();
    const uint64_t M_bytes = e.get_num_output();
    
    uint8_t buf_in[K_bytes];
    uint8_t buf_out[M_bytes];
    
    for(size_t i=0; i<K_bytes; i++) {
        buf_in[i] = 0xFF;
    }
    
    e.encode(buf_out, buf_in);
    printf("AR4JA_r12_k1024n:\n");
    
    for(size_t i=0; i<M_bytes; i++) {
        printf("Bit %4lu - %4lu: %02X\n", i*8, i*8+7, buf_out[i]);
    }
    
}

int main(void) {
    
    printf("martin_test:\n");
    test01();
    
    
    test02();
    
}
