#include "SeqLib/BamReader.h"
#include "SeqLib/SeqPlot.h"
#include "catch.hh"

#define INBAM "../htslib/test/range.bam"

using namespace SeqLib;

using namespace std;

TEST_CASE( "BAM header is parsed", "[bam-reader]" ) {
    BamReader br;
    br.Open(INBAM);
    BamHeader h = br.Header();
    std::string hdr{"@HD\tVN:1.4\tSO:coordinate\n@RG\tID:1\tPL:ILLUMINA\tPU:130410_HS18_09653_A_C1JT2ACXX_4\tLB:7053878\tDT:2013-04-10T00:00:00+0100\tSM:ERS225193\tCN:SC\n@SQ\tSN:CHROMOSOME_I\tLN:1009800\tM5:8ede36131e0dbf3417807e48f77f3ebd\tUR:/\n@SQ\tSN:CHROMOSOME_II\tLN:5000\tM5:8e7993f7a93158587ee897d7287948ec\tUR:/\n@SQ\tSN:CHROMOSOME_III\tLN:5000\tM5:3adcb065e1cf74fafdbba1e8c352b323\tUR:/\n@SQ\tSN:CHROMOSOME_IV\tLN:5000\tM5:251af66a69ee589c9f3757340ec2de6f\tUR:/\n@SQ\tSN:CHROMOSOME_V\tLN:5000\tM5:cf200a65fb754836dcc56b24b3170ee8\tUR:/\n@SQ\tSN:CHROMOSOME_X\tLN:5000\tM5:6f9368fd2192c89c613718399d2d31fc\tUR:/\n@SQ\tSN:CHROMOSOME_MtDNA\tLN:5000\tM5:cd05857ece6411f40257a565ccfe15bb\tUR:/\n@PG\tID:scramble\tPN:scramble\tVN:1.14.7\tCL:scramble -M -I sam -s 50 -r /tmp/ce.fa - /tmp/ERR304769_subset.cram \n"};
    REQUIRE(h.AsString() == hdr);
}

TEST_CASE( "Plot a collection of gapped alignments", "[bam-reader]" ) {
    SeqLib::BamReader br;
    br.Open(INBAM);
    GenomicRegion gr("CHROMOSOME_II:1135-2984", br.Header());
    br.SetRegion(gr);

    SeqPlot s;
    s.SetView(gr);

    BamRecord rec;
    BamRecordVector brv;
    while (br.GetNextRecord(rec))
        if (!rec.CountNBases() && rec.MappedFlag())
            brv.push_back(rec);
    s.SetPadding(5);

    string plotrecs = s.PlotAlignmentRecords(brv);
    size_t n = std::count(plotrecs.begin(), plotrecs.end(), '\n');
    REQUIRE(n == 5);
    cout << plotrecs << endl;
}
