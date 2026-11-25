// Prints the iteration order of Clair3's khash counters (string KH_COUNTER and int
// KH_INT_COUNTER) for key sequences, and glibc rand() outputs, as fixtures for snv-fullalign.
// Input (khash mode): one case per line, "str k1 k2 ..." or "int k1 k2 ...". Output: per case,
// "key=count" in bucket order. `rand N` prints the first N rand() values of an unseeded process.
// Build: cc -O1 -I<Clair3>/src tools/khash_oracle_fullalign.c <Clair3>/src/medaka_khcounter.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "medaka_khcounter.h"

int main(int argc, char **argv) {
    if (argc == 3 && strcmp(argv[1], "rand") == 0) {
        long n = atol(argv[2]);
        for (long i = 0; i < n; i++) printf("%d\n", rand());
        return 0;
    }
    static char line[1 << 20];
    while (fgets(line, sizeof line, stdin)) {
        char *kind = strtok(line, " \n");
        int first = 1;
        if (strcmp(kind, "int") == 0) {
            khash_t(KH_INT_COUNTER) *h = kh_init(KH_INT_COUNTER);
            for (char *tok = strtok(NULL, " \n"); tok; tok = strtok(NULL, " \n"))
                kh_int_counter_add(h, atoi(tok), 1);
            for (khiter_t k = kh_begin(h); k != kh_end(h); ++k)
                if (kh_exist(h, k)) { printf(first ? "%d=%d" : " %d=%d", kh_key(h, k), kh_val(h, k)); first = 0; }
            kh_int_counter_destroy(h);
        } else {
            khash_t(KH_COUNTER) *h = kh_init(KH_COUNTER);
            for (char *tok = strtok(NULL, " \n"); tok; tok = strtok(NULL, " \n"))
                kh_counter_add(h, tok, 1);
            for (khiter_t k = kh_begin(h); k != kh_end(h); ++k)
                if (kh_exist(h, k)) { printf(first ? "%s=%d" : " %s=%d", kh_key(h, k), kh_val(h, k)); first = 0; }
            kh_counter_destroy(h);
        }
        printf("\n");
    }
    return 0;
}
