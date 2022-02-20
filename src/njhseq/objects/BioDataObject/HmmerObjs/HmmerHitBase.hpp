#pragma once

/*
 * HmmerHitBase.hpp
 *
 *  Created on: Feb 19, 2022
 *      Author: nick
 */

#include "njhseq/common.h"
#include "njhseq/utils.h"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"

namespace njhseq {





class HmmerHitBase{
public:


	uint32_t hmmFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t hmmTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	uint32_t alignFrom_ {std::numeric_limits<uint32_t>::max()};//1-based as per file spec
	uint32_t alignTo_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec

	uint32_t envFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t envTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	bool isReverseStrand() const;

	uint32_t env0BasedPlusStrandStart() const;
	uint32_t env0BasedPlusStrandEnd() const;
	uint32_t align0BasedPlusStrandStart() const;
	uint32_t align0BasedPlusStrandEnd() const;

	uint32_t zeroBasedHmmFrom() const; //per file specs the positions are 1-based

	uint32_t envLen() const;
	uint32_t hmmLen() const;
	uint32_t alignLen() const;

	virtual Json::Value toJson()const;

	virtual ~HmmerHitBase(){

	}
};




}  // namespace njhseq

