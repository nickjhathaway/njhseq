#pragma once
/*
 * HmmerDomainHitTab.hpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */



#include "bibseq/common.h"
#include "bibseq/utils.h"

namespace bibseq {


class HmmerDomainHitTab{
public:
	HmmerDomainHitTab();
	HmmerDomainHitTab(const std::string & line);

	std::string targetName_{""};
	std::string targetAcc_{""};
	uint32_t targetLen_ {std::numeric_limits<uint32_t>::max()};

	std::string queryName_{""};
	std::string queryAcc_{""};
	uint32_t queryLen_ {std::numeric_limits<uint32_t>::max()};

	double seqEvalue_ {std::numeric_limits<double>::max()};
	double seqScore_ {std::numeric_limits<double>::max()};
	double seqBias_ {std::numeric_limits<double>::max()};

	uint32_t domainId_ {std::numeric_limits<uint32_t>::max()};
	uint32_t domainsTotals_ {std::numeric_limits<uint32_t>::max()};
	double domain_c_evalue_ {std::numeric_limits<double>::max()};
	double domain_i_evalue_ {std::numeric_limits<double>::max()};
	double domainScore_ {std::numeric_limits<double>::max()};
	double domainBias_ {std::numeric_limits<double>::max()};

	uint32_t hmmFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t hmmTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	uint32_t alignFrom_ {std::numeric_limits<uint32_t>::max()};//1-based as per file spec
	uint32_t alignTo_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec

	uint32_t envFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t envTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	double acc_{std::numeric_limits<double>::max()};
	std::string targetDesc_{""};


	Json::Value toJson() const;





};


}  // namespace bibseq



