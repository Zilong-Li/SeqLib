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
        header->hdr = bcf_hdr_init("w");
        ret = header->SetVersion(version);
        if (ret != 0)
            throw std::runtime_error("couldn't set the version correctly.\n");
    }

    void BcfWriter::WriteLine(const std::string& vcfline)
    {
        std::vector<char> line(vcfline.begin(), vcfline.end());
        line.push_back('\0'); // don't forget string has no \0;
        s.s = &line[0];
        s.l = vcfline.length();
        s.m = vcfline.length();
        ret = vcf_parse(&s, header->hdr, b);
        if (ret > 0)
            throw std::runtime_error("error parsing: " + vcfline + "\n");
        ret = bcf_write(fp, header->hdr, b);
        if (ret != 0)
            throw std::runtime_error("error writing: " + vcfline + "\n");
        if (b->errcode == BCF_ERR_CTG_UNDEF)
        {
            header->AddContig(std::string(bcf_hdr_id2name(header->hdr, b->rid)));
            b->errcode = 0;
        }
    }

} // namespace SeqLib
