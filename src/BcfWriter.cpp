#include "SeqLib/BcfWriter.h"

namespace SeqLib
{
    // infer the mode by fname suffix
    BcfWriter::BcfWriter(const std::string& fname_) : fname(fname_)
    {
        std::string mode{"w"};
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
        if (endWith(fname, "bcf.gz"))
            mode += "b";
        if (endWith(fname, "bcf"))
            mode += "bu";
        if (endWith(fname, "vcf.gz"))
            mode += "z";
        fp = hts_open(fname.c_str(), mode.c_str());
    }

    BcfWriter::BcfWriter(const std::string& fname_, const std::string& mode) : fname(fname_)
    {
        fp = hts_open(fname.c_str(), mode.c_str());
    }

    BcfWriter::~BcfWriter()
    {
        hts_close(fp);
        bcf_destroy(b);
    }

    void BcfWriter::InitalHeader(std::string version)
    {
        header = BcfHeader();
        header.hdr = bcf_hdr_init("w");
        ret = header.SetVersion(version);
        if (ret != 0)
            throw std::runtime_error("couldn't set the version correctly.\n");
    }

    void BcfWriter::InitalHeader(const BcfHeader& h)
    {
        header = BcfHeader();
        header.hdr = bcf_hdr_dup(h.hdr); // make a copy of given header
        header.nsamples = bcf_hdr_nsamples(header.hdr);
        if (header.hdr == NULL)
            throw std::runtime_error("couldn't copy the header from another vcf.\n");
    }


    bool BcfWriter::WriteHeader()
    {
        ret = bcf_hdr_write(fp, header.hdr);
        if (ret == 0)
            return isHeaderWritten = true;
        else
            return false;
    }

    void BcfWriter::WriteLine(const std::string& vcfline)
    {
        if (!isHeaderWritten)
            WriteHeader();
        std::vector<char> line(vcfline.begin(), vcfline.end());
        line.push_back('\0'); // don't forget string has no \0;
        s.s = &line[0];
        s.l = vcfline.length();
        s.m = vcfline.length();
        ret = vcf_parse(&s, header.hdr, b);
        if (ret > 0)
            throw std::runtime_error("error parsing: " + vcfline + "\n");
        if (b->errcode == BCF_ERR_CTG_UNDEF)
        {
            printf("contig id=%s is not in the header. please use header->AddContig() to add the contig first!", bcf_hdr_id2name(header.hdr, b->rid));
            throw std::runtime_error("contig id not found in the header. please run header->AddContig() first.\n");
        }
        ret = bcf_write(fp, header.hdr, b);
        if (ret != 0)
            throw std::runtime_error("error writing: " + vcfline + "\n");
    }

} // namespace SeqLib
