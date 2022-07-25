#include "SeqLib/BcfRecord.h"

namespace SeqLib
{
    BcfHeader::BcfHeader()
    {
    }
    BcfHeader::~BcfHeader()
    {
        if (hdr)
            bcf_hdr_destroy(hdr);
    }
    BcfRecord::BcfRecord()
    {
    }

    BcfRecord::~BcfRecord()
    {
    }

} // namespace SeqLib
