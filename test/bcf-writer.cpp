#include "SeqLib/BcfReader.h"
#include "SeqLib/BcfWriter.h"
#include "catch.hh"
#include <iostream>

using namespace SeqLib;
using namespace std;

TEST_CASE("Test bcf-writer", "[bcf-writer]")
{
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    std::cout << br.header->AsString() << std::endl;

    br.header->AddSample("NA11111");
    try {
        auto vec = br.header->GetSeqnames();
    } catch (std::runtime_error & e) {
        std::cout << e.what() << "start adding contig 1 in header\n\n";
        br.header->AddContig("1");
    }
    std::cout << br.header->AsString() << std::endl;
    auto vec = br.header->GetSeqnames();
    for (auto& i : vec)
    {
        std::cout << i << std::endl;
    }

    BcfWriter bw("test.vcf");
    bw.WriteHeader(br.header);

}

TEST_CASE("Write BCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.bcf");
    bw.InitalHeader();
    bw.header->AddFormat("GT", "1", "String", "Genotype");
    bw.header->AddInfo("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header->AddContig("chr20");
    bw.header->AddSample("NA12878");
    bw.WriteLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
}

TEST_CASE("Write VCF with custome header and variants", "[bcf-writer]")
{
    BcfWriter bw("test.vcf");
    bw.InitalHeader("VCF4.3");
    bw.header->AddFormat("GT", "1", "String", "Genotype");
    bw.header->AddInfo("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.header->AddContig("chr20");
    bw.header->AddSample("NA12878");
    bw.WriteLine("chr20\t2006060\trs146931526\tG\tC\t100\tPASS\tAF=0.000998403\tGT\t1|0");
}

TEST_CASE("Write VCF by copying header from another VCF", "[bcf-writer]")
{
    BcfWriter bw("test.vcf.gz");
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    bw.WriteHeader(br.header);
}

TEST_CASE("Write VCF with custome header", "[bcf-writer]")
{
    BcfWriter bw("test.vcf");
    bw.InitalHeader();
    bw.header->AddContig("21");
    bw.header->AddFormat("GT", "1", "String", "Genotype");
    bw.header->AddInfo("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bw.WriteHeader();
}
