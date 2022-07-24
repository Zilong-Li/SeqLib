#ifndef BCFRECORD_H_
#define BCFRECORD_H_

#include <string>
#include <type_traits>
#include <vector>

extern "C"
{
#include "htslib/htslib/vcf.h"
}


namespace SeqLib
{

    class BcfRecord
    {
        friend class BcfReader;

    public:
        BcfRecord();
        ~BcfRecord();

        void Init(const bcf_hdr_t* hdr_, int nsamples_);

        template <class T>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<bool>>::value || std::is_same<T, std::vector<int>>::value, bool>::type
        GetGenotypes(T& gv)
        {
            ret = bcf_get_genotypes(hdr, line, &gts, &ndst);
            if (ret <= 0)
                return 0; // gt not present
            gv.resize(ret);
            nploidy = ret / nsamples;
            int i, j, k = 0, nphased = 0;
            for (i = 0; i < nsamples; i++)
            {
                for (j = 0; j < nploidy; j++)
                {
                    gv[k++] = bcf_gt_allele(gts[j + i * nploidy]) != 0; // only parse 0 and 1, ie max(nploidy)=2; other values 2,3... will be converted to 1;
                }
                nphased += (gts[1 + i * nploidy] & 1) == 1;
            }
            if (nphased == nsamples)
                isAllPhased = true;
            return 1;
        }

        // return a array for the requested field
        template <typename T>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int32_t>>::value || std::is_same<T, std::vector<float>>::value, bool>::type
        GetFormat(T& v, const std::string& tag)
        {
            fmt = bcf_get_fmt(hdr, line, tag.c_str());
            shape1 = fmt->n;
            if (std::is_same<T, std::vector<int32_t>>::value)
            {
                ret = bcf_get_format_int32(hdr, line, tag.c_str(), &buf, &ndst);
            }
            else if (std::is_same<T, std::vector<char>>::value)
            {
                ret = bcf_get_format_char(hdr, line, tag.c_str(), &buf, &ndst);
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                ret = bcf_get_format_float(hdr, line, tag.c_str(), &buf, &ndst);
            }
            if (ret >= 0)
            {
                v.resize(ret);
                _GetFormat(v);
                return true;
            }
            else
            {
                return false;
            }
        }

        template <typename T>
        void _GetFormat(T& v)
        {
            int i, j, k = 0;
            for (i = 0; i < nsamples; i++)
            {
                for (j = 0; j < fmt->n; j++)
                {
                    // https://stackoverflow.com/questions/16687172/cast-a-value-using-decltype-is-it-possible
                    v[k++] = static_cast<typename T::value_type*>(buf)[j + i * fmt->n];
                }
            }
        }

        bool isAllPhased = false;
        int nploidy = 0;
        int nsamples;
        int shape1;

    private:
        bcf1_t* line = bcf_init(); // current bcf record
        const bcf_hdr_t* hdr;
        bcf_fmt_t* fmt = NULL;
        int32_t* gts = NULL;
        void* buf = NULL;
        int ndst = 0;
        int ret;
    };

} // namespace SeqLib

#endif // BCFRECORD_H_
