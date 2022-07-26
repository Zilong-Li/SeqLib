#ifndef BCFREADER_H_
#define BCFREADER_H_

#include "SeqLib/BcfRecord.h"

extern "C"
{
#include "htslib/htslib/tbx.h"
#include "htslib/htslib/vcf.h"
}


namespace SeqLib
{
    class BcfReader
    {
    public:
        BcfReader(const std::string& fname_);

        /**
         *  @param samples  samples to include or exclude as a comma-separated string.
         *              LIST        .. select samples in list
         *              ^LIST       .. exclude samples from list
         *              -           .. include all samples
         *              NULL        .. exclude all samples
         */
        BcfReader(const std::string& fname_, const std::string& samples);

        BcfReader(const std::string& fname_, const std::string& samples, const std::string& region);

        ~BcfReader();

        inline int SetThreads(int n)
        {
            return hts_set_threads(fp, n);
        }

        void SetRegion(const std::string& region);

        // reset current region. seek to the first record.
        void Reset();

        bool GetNextVariant(BcfRecord& r);

        std::string fname;
        bool isBcf; // if the input file is bcf or vcf;
        std::shared_ptr<BcfHeader> header = std::make_shared<BcfHeader>(); // bcf header
    private:
        htsFile* fp = NULL;         // hts file
        hts_idx_t* hidx = NULL;     // hts index file
        tbx_t* tidx = NULL;         // .tbi .csi index file for vcf files
        hts_itr_t* itr = NULL;      // hts records iterator
        kstring_t s = {0, 0, NULL}; // kstring
    };

} // namespace SeqLib

#endif // BCFREADER_H_
