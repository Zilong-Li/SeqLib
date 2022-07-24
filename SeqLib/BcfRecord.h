#ifndef BCFRECORD_H_
#define BCFRECORD_H_

#include <vector>
#include <string>

extern "C"
{
#include "htslib/htslib/vcf.h"
}


namespace SeqLib {

    class BcfRecord
    {
        friend class BcfReader;
    public:
        BcfRecord();
        ~BcfRecord();

        void Init(bcf_hdr_t* hdr_, int nsamples_);

        template <typename T>
        bool GetGenotypes(T& gv)
        {
            static_assert(std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<bool>>::value || std::is_same<T, std::vector<int>>::value, "the vector must be of char, bool or int type");

            if (bcf_get_genotypes(hdr, line, &gts, &ngt) <= 0)
                return 0; // gt not present
            gv.resize(ngt);
            /* int nsamples = bcf_hdr_nsamples(hdr); */
            nploidy = ngt / nsamples;
            int i, j, k = 0;
            for (i = 0; i < nsamples; i++)
            {
                for (j = 0; j < nploidy; j++)
                {
                    gv[k++] = bcf_gt_allele(gts[j + i * nploidy]) != 0;
                }
            }
            return 1;
        }

    private:
        bcf1_t* line = bcf_init(); // current bcf record
        bcf_hdr_t* hdr;
        int32_t *gts = NULL, ngt = 0, nploidy = 0;
        int nsamples;
    };

}

#endif // BCFRECORD_H_
