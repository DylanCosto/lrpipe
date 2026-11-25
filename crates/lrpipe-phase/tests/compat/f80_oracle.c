// Oracle for compat::f80: random long double operands and the x87 results of
// +, -, *, / and the conversion to double, one line per case:
//   a b a+b a-b a*b a/b (double)a
// as 20-digit hex (sign/exponent, significand). Build and run:
//   gcc -O0 -o f80_oracle tests/compat/f80_oracle.c && ./f80_oracle SEED COUNT > cases.txt
//   cargo run --release --example f80_replay < cases.txt
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static uint64_t s;
static uint64_t next(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
static void pr(long double v) {
    uint64_t m; uint16_t e;
    memcpy(&m, &v, 8); memcpy(&e, (char *)&v + 8, 2);
    printf("%04x%016llx ", e, (unsigned long long)m);
}
static long double rnd(void) {
    long double v;
    uint64_t m = next() | (1ull << 63);
    int kind = next() % 4;
    // exponents near 0 (probabilities), sometimes close to each other
    uint16_t e = 16383 - (kind == 0 ? next() % 4 : next() % 200) + (next() % 3 == 0 ? 0 : 0);
    if (kind == 3) m = (next() % 1000 + 1) << 54 | (1ull << 63);
    memcpy(&v, &m, 8); memcpy((char *)&v + 8, &e, 2);
    if (next() % 8 == 0) v = -v;
    return v;
}
int main(int argc, char **argv) {
    s = strtoull(argv[1], 0, 10) * 2654435761u + 1;
    long n = strtol(argv[2], 0, 10);
    for (long i = 0; i < n; i++) {
        volatile long double a = rnd(), b = rnd();
        if (i % 5 == 0) b = a * (1 + 1e-15L * (long double)(next() % 7));
        pr(a); pr(b); pr(a + b); pr(a - b); pr(a * b); pr(a / b);
        double d = (double)a;
        long double dd = d;
        pr(dd);
        printf("\n");
    }
    return 0;
}
