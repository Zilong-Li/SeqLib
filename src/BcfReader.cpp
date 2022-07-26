#include "SeqLib/BcfReader.h"


namespace SeqLib
{
    BcfReader::BcfReader(const std::string& fname_) : fname(fname_)
    {
        fp = hts_open(fname.c_str(), "r");
        header->hdr = bcf_hdr_read(fp);
        nsamples = bcf_hdr_nsamples(header->hdr);
    }

    BcfReader::BcfReader(const std::string& fname_, const std::string& samples) : fname(fname_)
    {
        fp = hts_open(fname.c_str(), "r");
        header->hdr = bcf_hdr_read(fp);
        header->SetSamples(samples);
        nsamples = bcf_hdr_nsamples(header->hdr);
    }

    BcfReader::BcfReader(const std::string& fname_, const std::string& samples, const std::string& region) : fname(fname_)
    {
        fp = hts_open(fname.c_str(), "r");
        header->hdr = bcf_hdr_read(fp);
        header->SetSamples(samples);
        nsamples = bcf_hdr_nsamples(header->hdr);
        SetRegion(region);
    }

    BcfReader::~BcfReader()
    {
        if (fp)
            hts_close(fp);
        if (itr)
            hts_itr_destroy(itr);
    }

    // get next variant and automatically update header if something is missing
    bool BcfReader::GetNextVariant(BcfRecord& r)
    {
        int ret;
        if (itr != NULL)
        {
            if (isBcf)
            {
                ret = bcf_itr_next(fp, itr, r.line);
                return (ret >= 0);
            }
            else
            {
                int slen = tbx_itr_next(fp, tidx, itr, &s);
                if (slen > 0)
                {
                    ret = vcf_parse(&s, header->hdr, r.line); // ret > 0, error
                }
                return (ret <= 0) && (slen > 0);
            }
        }
        else
        {
            ret = bcf_read(fp, header->hdr, r.line);
            return (ret == 0);
        }
    }

    // 1. check and load index first
    // 2. query iterval region
    void BcfReader::SetRegion(const std::string& region)
    {
        auto endWith = [](std::string const& s, std::string e)
        {
            if (s.length() >= e.length())
            {
                return (0 == s.compare(s.length() - e.length(), e.length(), e));
            }
            else
            {
                return false;
            }
        };

        if (endWith(fname, "bcf") || endWith(fname, "bcf.gz"))
        {
            isBcf = true;
            hidx = bcf_index_load(fname.c_str());
            itr = bcf_itr_querys(hidx, header->hdr, region.c_str());
        }
        else
        {
            isBcf = false;
            tidx = tbx_index_load(fname.c_str());
            assert(tidx != NULL && "error loading tabix index!");
            itr = tbx_itr_querys(tidx, region.c_str());
            assert(itr != NULL && "no interval region found.failed!");
        }
    }

} // namespace SeqLib
