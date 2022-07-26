#include "SeqLib/BcfReader.h"
#include "catch.hh"
#include <iostream>

using namespace SeqLib;
using namespace std;

TEST_CASE("Parsing VCF with specific tag", "[bcf-reader]")
{
    // BcfReader br("../htslib/test/index.vcf");
    BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    BcfRecord v(br.header);
    vector<int> ad;
    vector<char> gatk;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetFormat(ad, "AD");
        for (auto i : ad) {
            cout << i << ":";
        }
        cout << endl;
        v.GetFormat(gatk, "GATK");
        for (int i = 0; i < v.nsamples; i++) {
            cout << string(gatk.begin()+ i*v.shape1, gatk.begin() + i * v.shape1 + v.shape1 - 1) << ":";
        }
        cout << endl;
        n++;
    }
    REQUIRE(n == 10);
}

TEST_CASE("Parsing VCF", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing VCF with subset samples", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing VCF in target region", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz", "-", "chr20:2006060");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3243);
}

TEST_CASE("Parsing VCF with subset samples in target region", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116", "chr20:2006060-");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3243);
}

TEST_CASE("Parsing BCF", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.bcf.gz");
    BcfRecord v(br.header);
    vector<int> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing BCF in target region", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz", "-", "chr20");
    BcfRecord v(br.header);
    vector<int> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing BCF with subset samples", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.bcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing BCF with subset samples in target region", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.bcf.gz", "HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116", "chr20:2006060-2010000");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.GetNextVariant(v))
    {
        v.GetGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.isAllPhased) n++;
    }
    REQUIRE(n == 129);
}
