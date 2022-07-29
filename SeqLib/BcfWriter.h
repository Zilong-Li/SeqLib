#ifndef BCFWRITER_H_
#define BCFWRITER_H_

#include "SeqLib/BcfRecord.h"


namespace SeqLib
{
    class BcfWriter
    {
    public:
        BcfWriter(const std::string& fname_);

        /*!
          @abstract       Open a sequence data (SAM/BAM/CRAM) or variant data (VCF/BCF)
                          or possibly-compressed textual line-orientated file
          @param fn       The file name or "-" for stdin/stdout. For indexed files
                          with a non-standard naming, the file name can include the
                          name of the index file delimited with HTS_IDX_DELIM
          @param mode     Mode matching / [rwa][bcefFguxz0-9]* /
          @example
              [rw]b  .. compressed BCF, BAM, FAI
              [rw]bu .. uncompressed BCF
              [rw]z  .. compressed VCF
              [rw]   .. uncompressed VCF
        */
        BcfWriter(const std::string& fname_, const std::string& mode);

        ~BcfWriter();

        void InitalHeader(std::string version = "VCF4.1");

        // make a copy of given header
        void InitalHeader(const BcfHeader& h);

        void WriteLine(const std::string& vcfline);

        bool WriteHeader();

        inline bool WriteRecord(BcfRecord& v)
        {
            if (!isHeaderWritten)
                WriteHeader();
            if (bcf_write(fp, v.header->hdr, v.line) < 0)
                return false;
            else
                return true;
        }


        BcfHeader header; // bcf header

    private:
        htsFile* fp = NULL; // hts file
        int ret;
        bcf1_t* b = bcf_init();
        kstring_t s = {0, 0, NULL}; // kstring
        std::string fname;
        bool isHeaderWritten = false;
    };
} // namespace SeqLib

#endif // BCFWRITER_H_
