// Prints what Clair3's SSW (preprocess/realign/ssw_cpp.cpp, default Aligner) returns for each
// "reference<TAB>query" line on stdin; see tools/parity_realign.sh.
#include "ssw_cpp.h"
#include <iostream>
#include <string>
int main() {
  std::string line;
  while (std::getline(std::cin, line)) {
    size_t tab = line.find('\t');
    std::string ref = line.substr(0, tab), query = line.substr(tab + 1);
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment a;
    aligner.SetReferenceSequence(ref.c_str(), ref.length());
    if (!aligner.Align(query.c_str(), filter, &a)) { std::cout << "none\n"; continue; }
    std::cout << a.sw_score << '\t' << a.ref_begin << '\t' << a.ref_end << '\t' << a.query_begin << '\t' << a.query_end << '\t' << a.cigar_string << '\n';
  }
}
