
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include <stdlib.h>

/*  
modular reduction of the integer a with respect to the modulus mod.
It ensures that the result lies in the range (-mod/2, mod/2].
*/
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

/*  
modular exponentiation a^b mod m.
square-and-multiply
*/
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

/*  
Cooley-Tukey
It takes two complex numbers a and b, and a twiddle factor twiddle, 
and computes a' = a + twiddle * b and b' = a - twiddle * b, 
where a' and b' are the updated values of a and b, respectively.
*/
static void butterfly(int16_t *a, int16_t *b, int16_t twiddle, int16_t m) {
    int16_t t = cmod((int32_t)*b * twiddle, m);
    *b = cmod((int32_t)*a - t, m);
    *a = cmod((int32_t)*a + t, m);
}

// Function to reverse the bits of a given index
uint16_t bit_reverse(uint16_t num, int16_t len) {
    uint16_t rev_num = 0;
    for (int16_t i = 0; i < len; i++) {
        rev_num = (rev_num << 1) | ((num >> i) & 1);
    }
    return rev_num;
}

void bit_reversal_reorder(int16_t *a, int16_t n, int16_t N_bit) {
    for (int16_t i = 0; i < n; i++) {
        int16_t rev_i = bit_reverse(i, N_bit);
        if (rev_i > i) {
            int16_t temp = a[i];
            a[i] = a[rev_i];
            a[rev_i] = temp;
        }
    }
}

static void fft(int16_t *a, int16_t n, int16_t m) {
    int16_t N_bit = 0;
    for (int16_t temp = n; temp > 0; temp >>= 1)
        N_bit++;
    N_bit--;

    bit_reversal_reorder(a, n, N_bit);

    // Cooley-Tukey FFT
    for (int16_t i = 0; i < N_bit; i++) {
        int16_t points1[n / 2];
        int16_t points2[n / 2];
        int16_t len = 1 << (i + 1);
        int16_t half_len = len >> 1;
        int16_t twiddle = modpow(NTRUP_ROOT, NTRUP_Q / len, NTRUP_Q);

        for (int16_t j = 0; j < n / 2; j++) {
            int16_t shift_bits = N_bit - i;
            int16_t P = (j >> shift_bits) << shift_bits;
            int16_t w_P = modpow(twiddle, P, NTRUP_Q);
            int16_t even = a[j];
            int16_t odd = a[j + half_len] * w_P;
            points1[j] = cmod(even + odd, m);
            points2[j] = cmod(even - odd, m);
        }

        int16_t k = 0;
        for (int16_t j = 0; j < n / 2; j++) {
            a[k++] = points1[j];
            a[k++] = points2[j];
        }
    }
}

static void ifft(int16_t *a, int16_t n, int16_t m) {
    int16_t inv_n = modpow(n, NTRUP_Q - 2, NTRUP_Q);
    int16_t inv_w = modpow(NTRUP_ROOT, NTRUP_Q - 1, NTRUP_Q);

    fft(a, n, m);

    for (int16_t i = 0; i < n; i++) {
        a[i] = cmod((int32_t)a[i] * inv_w, m);
        a[i] = cmod((int32_t)a[i] * inv_n, m);
    }

    bit_reversal_reorder(a, n, n / 2);
}

void fft_poly_mul(int16_t des[NTRUP_P], const int16_t src1[NTRUP_P], const int8_t src2[NTRUP_P]) {
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