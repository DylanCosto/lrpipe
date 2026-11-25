//! Emulation of implementation details of CPython and libstdc++ that affect
//! WhatsHap's results (tie-breaking orders). See the module docs.
pub mod cxx;
pub mod f80;
pub mod pyfloat;
pub mod pyset;
pub mod pysort;
pub mod pystr;
