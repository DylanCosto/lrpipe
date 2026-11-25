// Prints the iteration order of Clair3's khash string counter for key sequences.
// Input: one case per line, keys separated by spaces. Output: per case, "key=count" in bucket
// order. Build: cc -O1 -I<Clair3>/src tools/khash_oracle.c <Clair3>/src/medaka_khcounter.c
#include <stdio.h>
#include <string.h>
#include "medaka_khcounter.h"

int main(void) {
    static char line[1 << 20];
    while (fgets(line, sizeof line, stdin)) {
        khash_t(KH_COUNTER) *h = kh_init(KH_COUNTER);
        for (char *tok = strtok(line, " \n"); tok; tok = strtok(NULL, " \n"))
            kh_counter_increment(h, tok);
        int first = 1;
        for (khiter_t k = kh_begin(h); k != kh_end(h); ++k)
            if (kh_exist(h, k)) { printf(first ? "%s=%d" : " %s=%d", kh_key(h, k), kh_val(h, k)); first = 0; }
        printf("\n");
        kh_counter_destroy(h);
    }
    return 0;
}
