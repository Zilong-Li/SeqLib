#include "SeqLib/BcfReader.h"
#include "SeqLib/BcfWriter.h"
#include "catch.hh"
#include <iostream>

using namespace SeqLib;
using namespace std;

TEST_CASE("Write VCF", "[bcf-writer]")
{
    BcfWriter bw("test.vcf.gz");
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    bw.WriteHeader(br.header);

}
