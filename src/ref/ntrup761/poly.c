
#include <stdint.h>

#include "params.h"

#include "poly.h"

static int16_t cmod(int32_t a, int16_t mod){
    int16_t t;
    t = a % mod;
    if(t > (mod >> 1)){
        t -= mod;
    }
    if(t < -(mod >> 1)){
        t += mod;
    }
    return t;
}


void poly_Rq_mul_small(int16_t *h, const int16_t *f,const int8_t *g)
{
    int16_t fg[NTRUP_P + NTRUP_P - 1];
    int16_t result;
    int i,j;

    for (i = 0; i < NTRUP_P; i++) {
      result = 0;
      for (j = 0;j <= i;++j) 
          result = cmod(result+f[j]*(int32_t)g[i-j], NTRUP_Q);
      fg[i] = result;
    }
    for (i = NTRUP_P;i < 2 * NTRUP_P - 1;++i) {
      result = 0;
      for (j = i - NTRUP_P + 1; j < NTRUP_P; ++j) 
          result = cmod(result+f[j]*(int32_t)g[i-j], NTRUP_Q);
      fg[i] = result;
    }

    for (i = 2 * NTRUP_P - 2; i >= NTRUP_P; i--) {
      fg[i - NTRUP_P] = cmod(fg[i - NTRUP_P] + fg[i], NTRUP_Q);
      fg[i - NTRUP_P + 1] = cmod(fg[i - NTRUP_P + 1] + fg[i], NTRUP_Q);
    }

    for (i = 0; i < NTRUP_P; ++i) 
        h[i] = fg[i];
}

static int16_t modpow(int16_t a, int16_t b, int16_t m) {
    int16_t res = 1;
    while (b > 0) {
        if (b & 1) {
            res = cmod((int32_t)res * a, m);
        }
        a = cmod((int32_t)a * a, m);
        b >>= 1;
    }
    return res;
}

int16_t find_primitive_root(int16_t q) {
    printf("Finding primitive root for q = %d\n", q);
    int16_t phi = q - 1;
    int16_t n = phi;
    for (int16_t i = 2; i * i <= n; ++i) {
        if (n % i == 0) {
            phi -= phi / i;
            while (n % i == 0) {
                n /= i;
            }
        }
    }
    if (n > 1) {
        phi -= phi / n;
    }

    for (int16_t g = 2; g < q; ++g) {
        if (modpow(g, phi, q) == 1) {
            printf("Primitive root found: %d\n", g);
            return g;
        }
    }

    printf("No primitive root found\n");
    return -1; // No primitive root found
}

static void butterfly(int16_t *a, int16_t *b, int16_t twiddle, int16_t m) {
    int16_t t = cmod((int32_t)*b * twiddle, m);
    *b = cmod((int32_t)*a - t, m);
    *a = cmod((int32_t)*a + t, m);
}

static void fft(int16_t *a, int16_t n, int16_t m) {
    int16_t j = 0;
    for (int16_t i = 1; i < n; i++) {
        int16_t bit = n >> 1;
        for (; j >= bit; bit >>= 1) {
            j -= bit;
        }
        j += bit;
        if (i < j) {
            int16_t temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }
    }

    for (int16_t len = 2; len <= n; len <<= 1) {
        int16_t half_len = len >> 1;
        int16_t twiddle = modpow(NTRUP_ROOT, NTRUP_Q / len, NTRUP_Q);
        for (int16_t i = 0; i < n; i += len) {
            int16_t t = 1;
            for (int16_t j = 0; j < half_len; j++) {
                butterfly(&a[i + j], &a[i + j + half_len], t, m);
                t = cmod((int32_t)t * twiddle, m);
            }
        }
    }
}

static void ifft(int16_t *a, int16_t n, int16_t m) {
    int16_t inv_n = modpow(n, NTRUP_Q - 2, NTRUP_Q);
    fft(a, n, m);
    for (int16_t i = 0; i < n; i++) {
        a[i] = cmod((int32_t)a[i] * inv_n, m);
    }
    for (int16_t i = 0; i < n / 2; i++) {
        int16_t temp = a[i];
        a[i] = a[n - i];
        a[n - i] = temp;
    }
}

void fft_poly_mul(int16_t *des, const int16_t *src1, const int8_t *src2) {
    int16_t n = NTRUP_P;
    int16_t m = NTRUP_Q;

    int16_t *a = (int16_t *)malloc(2 * n * sizeof(int16_t));
    int16_t *b = (int16_t *)malloc(2 * n * sizeof(int16_t));

    for (int16_t i = 0; i < n; i++) {
        a[i] = src1[i];
        b[i] = src2[i];
    }
    for (int16_t i = n; i < 2 * n; i++) {
        a[i] = 0;
        b[i] = 0;
    }
    printf("test\n");
    fft(a, 2 * n, m);
    fft(b, 2 * n, m);

    for (int16_t i = 0; i < 2 * n; i++) {
        a[i] = cmod((int32_t)a[i] * b[i], m);
    }

    ifft(a, 2 * n, m);

    for (int16_t i = 0; i < n; i++) {
        des[i] = a[i];
    }

    free(a);
    free(b);
}
























