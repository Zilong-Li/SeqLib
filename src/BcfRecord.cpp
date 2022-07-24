#include "SeqLib/BcfRecord.h"

namespace SeqLib
{
    BcfRecord::BcfRecord()
    {
    }

    BcfRecord::~BcfRecord()
    {
    }

    // init some necessary members
    void BcfRecord::Init(bcf_hdr_t* hdr_, int nsamples_)
    {
        hdr = hdr_;
        nsamples = nsamples_;
    }
} // namespace SeqLib
