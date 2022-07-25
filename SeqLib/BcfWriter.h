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
        void InitalHeader(std::string version = "VCF4.3");
        bool WriteRecord(BcfRecord& r);
        void WriteLine(const std::string& vcfline);

        inline bool WriteHeader()
        {
            ret = bcf_hdr_write(fp, header->hdr);
            if (ret == 0)
                return true;
            else
                return false;
        }

        // use copy assignment of shared_ptr
        inline bool WriteHeader(const std::shared_ptr<BcfHeader>& h)
        {
            header = h;
            // printf("address of first header = %p\n", &h);
            // printf("address of second header = %p\n", &header);
            ret = bcf_hdr_write(fp, header->hdr);
            if (ret == 0)
                return true;
            else
                return false;
        }

        std::shared_ptr<BcfHeader> header = std::make_shared<BcfHeader>(); // bcf header

    private:
        htsFile* fp = NULL;                                                // hts file
        int ret;
        bcf1_t *b = bcf_init();
        kstring_t s = {0, 0, NULL}; // kstring
        std::string fname;
    };
} // namespace SeqLib

#endif // BCFWRITER_H_
