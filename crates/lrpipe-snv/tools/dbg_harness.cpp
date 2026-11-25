// Prints the haplotypes Clair3's get_consensus (preprocess/realign/debruijn_graph.cpp) returns for
// each "reference<TAB>reads<TAB>low-quality positions" line on stdin; see tools/parity_realign.sh.
#include <iostream>
#include <string>
struct struct_str_arr { int consensus_size; char* consensus[500]; };
extern "C" struct_str_arr* get_consensus(char* reference, char* c_reads, char* c_base_quality, int read_size);
int main() {
  std::string line;
  while (std::getline(std::cin, line)) {
    size_t a = line.find('\t'), b = line.find('\t', a + 1);
    std::string ref = line.substr(0, a), reads = line.substr(a + 1, b - a - 1), bq = line.substr(b + 1);
    struct_str_arr* r = get_consensus(&ref[0], &reads[0], &bq[0], 0);
    for (int i = 0; i < r->consensus_size; i++) std::cout << (i ? "," : "") << r->consensus[i];
    std::cout << '\n';
  }
}
