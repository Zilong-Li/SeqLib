#include "catch.hh"
#include <iostream>
#include "vcfpp/vcfpp.h"

using namespace std;
using namespace vcfpp;

TEST_CASE("Test bcf-writer", "[bcf-writer]")
{
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    std::cout << br.header << std::endl;

    br.header.addSample("NA11111");
    try {
        auto vec = br.header.getSeqnames();
    } catch (std::runtime_error & e) {
        std::cout << e.what() << "start adding contig 1 in header\n\n";
        br.header.addContig("1");
    }
    std::cout << br.header.asString() << std::endl;
    auto vec = br.header.getSeqnames();
    for (auto& i : vec)
    {
        std::cout << i << std::endl;
    }

    BcfWriter bw("test.vcf");
    bw.writeHeader();

}

TEST_CASE("Write BCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.bcf");
    bw.initalHeader();
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addContig("chr20");
    bw.header.addSample("NA12878");
    bw.writeLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
}

TEST_CASE("Write VCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.vcf");
    bw.initalHeader("VCF4.3");
    bw.header.addFORMAT("GT", "1", "String", "Genotype");
    bw.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header.addContig("chr20");
    bw.header.addSample("NA12878");
    bw.writeLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
}


TEST_CASE("Write VCF by copying header from another VCF", "[bcf-writer]")
{
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    BcfRecord v(br.header);
    BcfWriter bw("test.vcf");
    bw.initalHeader(br.header);
    bw.writeHeader();
    br.getNextVariant(v);
    bw.writeRecord(v);
}
