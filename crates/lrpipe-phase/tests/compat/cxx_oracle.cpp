// Oracle for compat::cxx. Reads commands from stdin:
//   H <string>     print std::hash<std::string>
//   N              new (empty) unordered_set<int>
//   I <int>        insert
//   C              clear
//   P              print iteration order and bucket count
//   G <n>          print the bucket counts seen while inserting 0..n into a fresh set
#include <iostream>
#include <string>
#include <unordered_set>
#include <functional>
int main() {
    std::string cmd;
    std::unordered_set<int>* s = new std::unordered_set<int>();
    while (std::cin >> cmd) {
        if (cmd == "H") { std::string x; std::cin >> x; std::cout << std::hash<std::string>()(x) << "\n"; }
        else if (cmd == "N") { delete s; s = new std::unordered_set<int>(); }
        else if (cmd == "I") { int v; std::cin >> v; s->insert(v); }
        else if (cmd == "C") { s->clear(); }
        else if (cmd == "P") { for (int v : *s) std::cout << v << " "; std::cout << "| " << s->bucket_count() << "\n"; }
        else if (cmd == "G") { long n; std::cin >> n; std::unordered_set<long> t; size_t last = 0;
            for (long i = 0; i < n; ++i) { t.insert(i); if (t.bucket_count() != last) { last = t.bucket_count(); std::cout << last << " "; } }
            std::cout << "\n"; }
    }
}
