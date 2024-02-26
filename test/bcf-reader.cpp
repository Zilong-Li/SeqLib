#include "vcfpp/vcfpp.h"
#include "catch.hh"
#include <iostream>

using namespace std;
using namespace vcfpp;

TEST_CASE("Parsing VCF with specific tag", "[bcf-reader]")
{
    // BcfReader br("../htslib/test/index.vcf");
    // BcfReader br("../htslib/test/test-vcf-hdr-in.vcf");
    BcfReader br("test-vcf-read.vcf");
    BcfRecord v(br.header);
    vector<float> ad;
    vector<int> gq;
    vector<char> gatk;
    vector<bool> gt;
    int n = 0;
    while (br.getNextVariant(v))
    {
        n++;
        v.getGenotypes(gt);
        for (auto i : gt) {
            cout << i << ":";
        }
        cout << endl;

        v.getFORMAT("AD", ad);
        for (auto i : ad) {
            cout << i << ":";
        }
        cout << endl;

        v.getFORMAT("GQ", gq);
        for (auto i : gq) {
            cout << i << ":";
        }
        cout << endl;

        v.getFORMAT("GATK", gatk);
        cout << string(gatk.begin(), gatk.end()) << endl;

    }

    BcfReader br2("../htslib/test/index.vcf");
    BcfRecord v2(br2.header);
    br2.getNextVariant(v2);
    v2.removeINFO("DP");
    br2.header.addINFO("Str", "1", "String", "this is a test for adding string in INFO");
    v2.setINFO("Str", string{"S1S2"});
    v2.setINFO("Str", string{"str"});
    vector<int> dp;
    v2.setINFO("DP", vector<int>{1,2});
    v2.getINFO("DP", dp);
    for (auto i : dp) {
        cout << i << endl;
    }
    v2.setINFO("DP", 2);
    v2.getINFO("DP", dp);
    for (auto i : dp) {
        cout << i << endl;
    }

    cout << br2.header << endl;
    cout << v2.asString() << endl;

    REQUIRE(n == 10);
}

TEST_CASE("Parsing VCF", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.allPhased()) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing VCF with subset samples", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz", "","HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        REQUIRE(gt.size() == 20);
        if(v.allPhased()) n++;
    }
    REQUIRE(n == 3786);
}

TEST_CASE("Parsing VCF in target region", "[bcf-reader]")
{
    BcfReader br("chr20.2000001.2100000.vcf.gz", "chr20:2006060", "-");
    BcfRecord v(br.header);
    vector<bool> gt;
    int n = 0;
    while (br.getNextVariant(v))
    {
        v.getGenotypes(gt);
        REQUIRE(gt.size() == 2504 * 2);
        if(v.allPhased()) n++;
    }
    REQUIRE(n == 3243);
}

