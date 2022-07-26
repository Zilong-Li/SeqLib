#ifndef BCFRECORD_H_
#define BCFRECORD_H_

#include <type_traits>

#include "SeqLib/SeqLibUtils.h"

extern "C"
{
#include "htslib/htslib/vcf.h"
}


namespace SeqLib
{
    class BcfHeader
    {
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfHeader()
        {
            hdr = NULL;
        }

        ~BcfHeader()
        {
            bcf_hdr_destroy(hdr);
        }

        // todo : check if the value is valid for vcf specification
        inline void AddInfo(const std::string& id, const std::string& number, const std::string& type, const std::string& description)
        {
            AddLine("##INFO=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description + "\">");
        }

        inline void AddFormat(const std::string& id, const std::string& number, const std::string& type, const std::string& description)
        {
            AddLine("##FORMAT=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description + "\">");
        }

        inline void AddFilter(const std::string& id, const std::string& number, const std::string& type, const std::string& description)
        {
            AddLine("##FILTER=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description + "\">");
        }

        inline void AddContig(const std::string& id)
        {
            AddLine("##contig=<ID=" + id + ">");
        }

        inline void AddLine(const std::string& line)
        {
            ret = bcf_hdr_append(hdr, line.c_str());
            if (ret != 0)
                throw std::runtime_error("could not add " + line + " to header\n");
            ret = bcf_hdr_sync(hdr);
            if (ret != 0)
                throw std::runtime_error("could not add " + line + " to header\n");
        }

        void SetSamples(const std::string& samples)
        {

            ret = bcf_hdr_set_samples(hdr, samples.c_str(), 0);
            if (ret > 0)
            {
                throw std::runtime_error("the " + std::to_string(ret) + "-th sample are not in the VCF.\n");
            }
            else if (ret == -1)
            {
                throw std::runtime_error("couldn't set samples. something wrong.\n");
            }
        }

        inline void AddSample(const std::string& sample)
        {
            bcf_hdr_add_sample(hdr, sample.c_str());
            if (bcf_hdr_sync(hdr) != 0)
            {
                throw std::runtime_error("couldn't add the sample.\n");
            }
        }

        inline int SetVersion(const std::string& version) const
        {
            return bcf_hdr_set_version(hdr, version.c_str());
        }

        inline char* AsString()
        {
            kstring_t s = {0, 0, NULL};          // kstring
            if (bcf_hdr_format(hdr, 0, &s) == 0) // append header string to s.s! append!
                return s.s;
            else
                throw std::runtime_error("failed to convert formatted header to string");
        }

        std::vector<std::string> GetSamples() const
        {
            std::vector<std::string> vec;
            for (int i = 0; i < bcf_hdr_nsamples(hdr); i++)
            {
                vec.emplace_back(std::string(hdr->samples[i]));
            }
            return vec;
        }

        std::vector<std::string> GetSeqnames()
        {
            const char** seqs = bcf_hdr_seqnames(hdr, &ret);
            if (ret == 0)
                printf("there is no contig id in the header!\n");
            std::vector<std::string> vec;
            for (int i = 0; i < ret; i++)
            {
                vec.emplace_back(std::string(seqs[i]));
            }
            // TODO: return uninitialized vec may be undefined.
            return vec;
        }

    private:
        bcf_hdr_t* hdr;          // bcf header
        bcf_hrec_t* hrec = NULL; // populate header
        int ret = 0;
    };

    class BcfRecord
    {
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfRecord()
        {
        }
        ~BcfRecord()
        {
        }

        inline void Init(bcf_hdr_t* hdr_, int nsamples_)
        {
            hdr = hdr_;
            nsamples = nsamples_;
        }

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

        template <class T>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<bool>>::value || std::is_same<T, std::vector<int>>::value, void>::type
        SetGenotypes(T& gv)
        {
            ret = bcf_update_genotypes(hdr, line, &gv[0], gv.size());
            if (ret < 0)
                throw std::runtime_error("couldn't set genotypes correctly.\n");
        }

        // return a array for the requested field
        template <typename T>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value || std::is_same<T, std::vector<float>>::value, bool>::type
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
                int i, j, k = 0;
                for (i = 0; i < nsamples; i++)
                {
                    for (j = 0; j < fmt->n; j++)
                    {
                        // https://stackoverflow.com/questions/16687172/cast-a-value-using-decltype-is-it-possible
                        v[k++] = static_cast<typename T::value_type*>(buf)[j + i * fmt->n];
                    }
                }
                return true;
            }
            else
            {
                return false;
            }
        }

        template <typename T>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value || std::is_same<T, std::vector<float>>::value, void>::type
        SetFormat(T& v, const std::string& tag)
        {
            if (std::is_same<T, std::vector<int32_t>>::value)
            {
                ret = bcf_update_format_int32(hdr, line, tag.c_str(), &v[0], v.size());
            }
            else if (std::is_same<T, std::vector<char>>::value)
            {
                ret = bcf_update_format_char(hdr, line, tag.c_str(), &v[0], v.size());
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                ret = bcf_update_format_float(hdr, line, tag.c_str(), &v[0], v.size());
            }
            if (ret < 0)
                throw std::runtime_error("couldn't set format correctly.\n");
        }

        bool isAllPhased = false;
        int nploidy = 0;
        int nsamples;
        int shape1;

    private:
        bcf1_t* line = bcf_init(); // current bcf record
        bcf_hdr_t* hdr = NULL;
        bcf_fmt_t* fmt = NULL;
        int32_t* gts = NULL;
        void* buf = NULL;
        int ndst = 0;
        int ret;
    };

} // namespace SeqLib

#endif // BCFRECORD_H_
