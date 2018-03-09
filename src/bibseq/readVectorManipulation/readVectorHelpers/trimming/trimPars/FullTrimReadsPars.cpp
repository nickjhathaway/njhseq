/*
 * FullTrimReadsPars.cpp
 *
 *  Created on: Dec 13, 2017
 *      Author: nick
 */


#include "FullTrimReadsPars.hpp"

namespace bibseq {

FullTrimReadsPars::FullTrimReadsPars(){
	tSeqPars_.includeSequence_ = false;
	tSeqPars_.sequenceToLowerCase_ = false;
	tSeqPars_.removePreviousSameBases_ = false;
	allowableErrors.distances_.query_.coverage_ = .75;
}

void FullTrimReadsPars::initForKSharedTrim(){
  allowableErrors.distances_.query_.coverage_ = .50;
  allowableErrors.hqMismatches_ = 2;
  allowableErrors.lqMismatches_ = 2;
  allowableErrors.oneBaseIndel_ = 2;
  allowableErrors.twoBaseIndel_ = 2;
  allowableErrors.largeBaseIndel_ = 2;
}


}  // namespace bibseq



