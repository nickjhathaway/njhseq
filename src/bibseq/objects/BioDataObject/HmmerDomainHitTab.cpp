/*
 * HmmerDomainHitTab.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */



#include "HmmerDomainHitTab.hpp"

namespace bibseq {

HmmerDomainHitTab::HmmerDomainHitTab(){

}


HmmerDomainHitTab::HmmerDomainHitTab(const std::string & line) {
	auto toks = tokenizeString(line, "whitespace");
	if (23 != toks.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error there should be 22 items not "
				<< toks.size() << "\n";
		for(const auto & pos : iter::range(toks.size())){
			ss << pos << ": " << toks[pos] << "\n";
		}
		throw std::runtime_error { ss.str() };
	}
	targetName_ = toks[0];
	targetAcc_ = toks[1];
	targetLen_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[2]);

	queryName_ = toks[3];
	queryAcc_ = toks[4];
	queryLen_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[5]);

	seqEvalue_ = bib::StrToNumConverter::stoToNum<double>(toks[6]);
	seqScore_ = bib::StrToNumConverter::stoToNum<double>(toks[7]);
	seqBias_ = bib::StrToNumConverter::stoToNum<double>(toks[8]);

	domainId_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[9]);
	domainsTotals_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[10]);

	domain_c_evalue_ = bib::StrToNumConverter::stoToNum<double>(toks[11]);
	domain_i_evalue_ = bib::StrToNumConverter::stoToNum<double>(toks[12]);
	domainScore_ = bib::StrToNumConverter::stoToNum<double>(toks[13]);
	domainBias_ = bib::StrToNumConverter::stoToNum<double>(toks[14]);

	hmmFrom_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[15]);
	hmmTo_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[16]);

	alignFrom_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[17]);
	alignTo_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[18]);

	envFrom_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[19]);
	envTo_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[20]);

	acc_ =  bib::StrToNumConverter::stoToNum<double>(toks[21]);
	targetDesc_ = toks[22];


}


Json::Value HmmerDomainHitTab::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["targetName_"] = bib::json::toJson(targetName_);
	ret["targetAcc_"] = bib::json::toJson(targetAcc_);
	ret["targetLen_"] = bib::json::toJson(targetLen_);
	ret["queryName_"] = bib::json::toJson(queryName_);
	ret["queryAcc_"] = bib::json::toJson(queryAcc_);
	ret["queryLen_"] = bib::json::toJson(queryLen_);
	ret["seqEvalue_"] = bib::json::toJson(seqEvalue_);
	ret["seqScore_"] = bib::json::toJson(seqScore_);
	ret["seqBias_"] = bib::json::toJson(seqBias_);
	ret["domainId_"] = bib::json::toJson(domainId_);
	ret["domainsTotals_"] = bib::json::toJson(domainsTotals_);
	ret["domain_c_evalue_"] = bib::json::toJson(domain_c_evalue_);
	ret["domain_i_evalue_"] = bib::json::toJson(domain_i_evalue_);
	ret["domainScore_"] = bib::json::toJson(domainScore_);
	ret["domainBias_"] = bib::json::toJson(domainBias_);
	ret["hmmFrom_"] = bib::json::toJson(hmmFrom_);
	ret["hmmTo_"] = bib::json::toJson(hmmTo_);
	ret["alignFrom_"] = bib::json::toJson(alignFrom_);
	ret["alignTo_"] = bib::json::toJson(alignTo_);
	ret["envFrom_"] = bib::json::toJson(envFrom_);
	ret["envTo_"] = bib::json::toJson(envTo_);
	ret["acc_"] = bib::json::toJson(acc_);
	ret["targetDesc_"] = bib::json::toJson(targetDesc_);
	return ret;
}


}  // namespace bibseq



