/*
 * HmmerHitBase.cpp
 *
 *  Created on: Feb 19, 2022
 *      Author: nick
 */


#include "HmmerHitBase.hpp"

namespace njhseq {

bool HmmerHitBase::isReverseStrand() const{
	return alignTo_ < alignFrom_;
}

uint32_t HmmerHitBase::env0BasedPlusStrandStart() const{
	uint32_t start = !isReverseStrand() ? envFrom_ - 1 : envTo_ - 1;
	return start;
}

uint32_t HmmerHitBase::env0BasedPlusStrandEnd() const{
	uint32_t end = !isReverseStrand() ? envTo_ : envFrom_;
	return end;
}

uint32_t HmmerHitBase::align0BasedPlusStrandStart() const{
	uint32_t start = !isReverseStrand() ? alignFrom_ - 1 : alignTo_ - 1;
	return start;
}

uint32_t HmmerHitBase::align0BasedPlusStrandEnd() const{
	uint32_t end = !isReverseStrand() ? alignTo_ : alignFrom_;
	return end;
}


uint32_t HmmerHitBase::zeroBasedHmmFrom() const{
	//per file specs the positions are 1-based
	return hmmFrom_ - 1;
}

uint32_t HmmerHitBase::envLen() const{
	return env0BasedPlusStrandEnd() - env0BasedPlusStrandStart();
}
uint32_t HmmerHitBase::hmmLen() const{
	return hmmTo_ - zeroBasedHmmFrom();
}
uint32_t HmmerHitBase::alignLen() const{
	return align0BasedPlusStrandEnd() - align0BasedPlusStrandStart();
}

Json::Value HmmerHitBase::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["hmmFrom_"] = njh::json::toJson(hmmFrom_);
	ret["hmmTo_"] = njh::json::toJson(hmmTo_);
	ret["alignFrom_"] = njh::json::toJson(alignFrom_);
	ret["alignTo_"] = njh::json::toJson(alignTo_);
	ret["envFrom_"] = njh::json::toJson(envFrom_);
	ret["envTo_"] = njh::json::toJson(envTo_);
	return ret;
}

}  // namespace njhseq
