#pragma once
/*
 * PairedRead.hpp
 *
 *  Created on: Jul 26, 2015
 *      Author: nick
 */

#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {

class PairedRead: public readObject {
public:
	PairedRead(const seqInfo & seqFirst, const seqInfo & seqSecond,
			bool processed) ;

	// just for holding seq and qual of paired read, cnt_, frac_, name_,
	// and on_ should be ignored and seqBase_ should be used instead
	seqInfo mateSeqBase_;

	virtual double getQualCheck(uint32_t qualCutOff) const;

	virtual void setBaseCountOnQualCheck(uint32_t qualCheck);
  virtual void setLetterCount();
  virtual void setLetterCount(const std::vector<char> & alph);

  virtual double getAverageErrorRate() const;

  virtual void createCondensedSeq();


  void outFastq(std::ostream & firstOut, std::ostream & secondOut,const readObjectIOOptions & ioOptions )const;

  using size_type = readObject::size_type;
};

template<>
inline PairedRead::size_type len(const PairedRead & read){
	return len(read.seqBase_) + len(read.mateSeqBase_);
}

} /* namespace bibseq */


