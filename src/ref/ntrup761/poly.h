#ifndef POLYMUL_H
#define POLYMUL_H

#include <stdint.h>

void poly_Rq_mul_small(int16_t *h, const int16_t *f, const int8_t *g);


// TODO: Unsure should we user another macro
int16_t find_primitive_root(int16_t q);
void fft_poly_mul(int16_t *des, const int16_t *src1, const int8_t *src2);

#endif



