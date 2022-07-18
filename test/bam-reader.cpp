#include "SeqLib/BamReader.h"
#include "catch.hh"

#define INSAM "../htslib/test/range.bam"

using namespace SeqLib;

TEST_CASE( "SAM file is parsed", "[bam-reader]" ) {
    SeqLib::BamReader br;
    br.Open(INSAM);
    SeqLib::BamHeader h = br.Header();
    std::cout << h.AsString() << std::endl;
}
