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
        friend class BcfRecord;
        friend class BcfReader;
        friend class BcfWriter;

    public:
        BcfHeader()
        {
            hdr = NULL;
        }

        ~BcfHeader()
        {
            // bcf_hdr_destroy(hdr); cause double free issue
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

        inline void AddLine(const std::string& str)
        {
            ret = bcf_hdr_append(hdr, str.c_str());
            if (ret != 0)
                throw std::runtime_error("could not add " + str + " to header\n");
            ret = bcf_hdr_sync(hdr);
            if (ret != 0)
                throw std::runtime_error("could not add " + str + " to header\n");
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

        std::string AsString() const
        {
            kstring_t s = {0, 0, NULL};          // kstring
            if (bcf_hdr_format(hdr, 0, &s) == 0) // append header string to s.s! append!
                return std::string(s.s, s.l);
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

        int nsamples = 0;

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
        BcfRecord(BcfHeader& h_) : header(std::make_shared<BcfHeader>(h_))
        {
        }

        ~BcfRecord()
        {
        }

        std::string AsString()
        {
            s.s = NULL;
            s.l = 0;
            s.m = 0;
            if (vcf_format(header->hdr, line, &s) == 0)
                return std::string(s.s, s.l);
            else
                throw std::runtime_error("couldn't format current record into a string.\n");
        }

        template <class T>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<bool>>::value || std::is_same<T, std::vector<int>>::value, bool>::type
        GetGenotypes(T& gv)
        {
            ret = bcf_get_genotypes(header->hdr, line, &gts, &ndst);
            if (ret <= 0)
                return 0; // gt not present
            gv.resize(ret);
            nploidy = ret / header->nsamples;
            int i, j, k = 0, nphased = 0;
            for (i = 0; i < header->nsamples; i++)
            {
                for (j = 0; j < nploidy; j++)
                {
                    gv[k++] = bcf_gt_allele(gts[j + i * nploidy]) != 0; // only parse 0 and 1, ie max(nploidy)=2; other values 2,3... will be converted to 1;
                }
                nphased += (gts[1 + i * nploidy] & 1) == 1;
            }
            if (nphased == header->nsamples)
                isAllPhased = true;
            return 1;
        }

        template <class T>
        typename std::enable_if<
            std::is_same<T, std::vector<bool>>::value || std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value, bool>::type
        SetGenotypes(const T& gv, bool phased = false)
        {
            ret = bcf_get_genotypes(header->hdr, line, &gts, &ndst);
            if (ret <= 0)
                return false; // gt not present
            assert(ret == gv.size());
            nploidy = ret / header->nsamples;
            int i, j, k;
            for (i = 0; i < header->nsamples; i++)
            {
                for (j = 0; j < nploidy; j++)
                {
                    k = i * nploidy + j;
                    if (phased)
                        gts[k] = bcf_gt_phased(gv[k]);
                    else
                        gts[k] = bcf_gt_unphased(gv[k]);
                }
            }
            if (bcf_update_genotypes(header->hdr, line, gts, ret) < 0)
                throw std::runtime_error("couldn't set genotypes correctly.\n");
            else
                return true;
        }

        // return a array for the requested field
        template <typename T, typename S = typename T::value_type>
        typename std::enable_if<
            std::is_same<T, std::vector<char>>::value || std::is_same<T, std::vector<int>>::value || std::is_same<T, std::vector<float>>::value, bool>::type
        GetFormat(T& v, const std::string& tag)
        {
            fmt = bcf_get_fmt(header->hdr, line, tag.c_str());
            shape1 = fmt->n;
            S* buf = NULL;
            if (std::is_same<T, std::vector<int>>::value)
            {
                ret = bcf_get_format_int32(header->hdr, line, tag.c_str(), &buf, &ndst);
            }
            else if (std::is_same<T, std::vector<char>>::value)
            {
                ret = bcf_get_format_char(header->hdr, line, tag.c_str(), &buf, &ndst);
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                ret = bcf_get_format_float(header->hdr, line, tag.c_str(), &buf, &ndst);
            }
            printf("fmt-n %d, nsaples %d \n", fmt->n, header->nsamples);
            if (ret >= 0)
            {
                v.resize(ret);
                int i, j, k = 0;
                for (i = 0; i < header->nsamples; i++)
                {
                    for (j = 0; j < fmt->n; j++)
                    {
                        // https://stackoverflow.com/questions/16687172/cast-a-value-using-decltype-is-it-possible
                        v[k++] = buf[j + i * fmt->n];
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
        SetFormat(const T& v, const std::string& tag)
        {
            if (std::is_same<T, std::vector<int32_t>>::value)
            {
                ret = bcf_update_format_int32(header->hdr, line, tag.c_str(), v.data(), v.size());
            }
            else if (std::is_same<T, std::vector<char>>::value)
            {
                ret = bcf_update_format_char(header->hdr, line, tag.c_str(), v.data(), v.size());
            }
            else if (std::is_same<T, std::vector<float>>::value)
            {
                ret = bcf_update_format_float(header->hdr, line, tag.c_str(), v.data(), v.size());
            }
            if (ret < 0)
                throw std::runtime_error("couldn't set format correctly.\n");
        }

        void AddLineFromString(const std::string& vcfline)
        {
            std::vector<char> str(vcfline.begin(), vcfline.end());
            str.push_back('\0'); // don't forget string has no \0;
            s.s = &str[0];
            s.l = vcfline.length();
            s.m = vcfline.length();
            ret = vcf_parse(&s, header->hdr, line);
            if (ret > 0)
                throw std::runtime_error("error parsing: " + vcfline + "\n");
            if (line->errcode == BCF_ERR_CTG_UNDEF)
            {
                std::string contig(bcf_hdr_id2name(header->hdr, line->rid));
                hdr_d = bcf_hdr_dup(header->hdr);
                header->hrec = bcf_hdr_id2hrec(hdr_d, BCF_DT_CTG, 0, line->rid);
                if (header->hrec == NULL)
                    throw std::runtime_error("contig" + contig + " unknow and not found in the header.\n");
                ret = bcf_hdr_add_hrec(header->hdr, header->hrec);
                printf("bcf_hdr_add_hrec %i \n", ret);
                if (ret < 0)
                    throw std::runtime_error("error adding contig " + contig + " to header.\n");
                ret = bcf_hdr_sync(header->hdr);
            }
        }

        bool isAllPhased = false;
        int nploidy = 0;
        int shape1;
        std::shared_ptr<BcfHeader> header;

    private:
        bcf1_t* line = bcf_init(); // current bcf record
        bcf_hdr_t* hdr_d;          // a dup header by bcf_hdr_dup(header->hdr)
        bcf_fmt_t* fmt = NULL;
        int32_t* gts = NULL;
        int ndst = 0;
        int ret;
        kstring_t s = {0, 0, NULL}; // kstring
    };

} // namespace SeqLib

#endif // BCFRECORD_H_
