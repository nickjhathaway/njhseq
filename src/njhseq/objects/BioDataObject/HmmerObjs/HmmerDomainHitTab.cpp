/*
 * HmmerDomainHitTab.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */



#include "HmmerDomainHitTab.hpp"

#include "njhseq/utils/vectorUtils.hpp"
#include "njhseq/objects/Meta/MetaDataInName.hpp"


namespace njhseq {

HmmerDomainHitTab::HmmerDomainHitTab(){

}


HmmerDomainHitTab::HmmerDomainHitTab(const std::string & line) {
	auto toks = tokenizeString(line, "whitespace");
	if (toks.size() < 23) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error there should be at least 23 items not "
				<< toks.size() << "\n";
		for(const auto pos : iter::range(toks.size())){
			ss << pos << ": " << toks[pos] << "\n";
		}
		throw std::runtime_error { ss.str() };
	}
	targetName_ = toks[0];
	targetAcc_ = toks[1];
	targetLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[2]);

	queryName_ = toks[3];
	queryAcc_ = toks[4];
	queryLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);

	seqEvalue_ = njh::StrToNumConverter::stoToNum<double>(toks[6]);
	seqScore_ = njh::StrToNumConverter::stoToNum<double>(toks[7]);
	seqBias_ = njh::StrToNumConverter::stoToNum<double>(toks[8]);

	domainId_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[9]);
	domainsTotals_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[10]);

	domain_c_evalue_ = njh::StrToNumConverter::stoToNum<double>(toks[11]);
	domain_i_evalue_ = njh::StrToNumConverter::stoToNum<double>(toks[12]);
	domainScore_ = njh::StrToNumConverter::stoToNum<double>(toks[13]);
	domainBias_ = njh::StrToNumConverter::stoToNum<double>(toks[14]);

	hmmFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[15]);
	hmmTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[16]);

	alignFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[17]);
	alignTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[18]);

	envFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[19]);
	envTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[20]);

	acc_ =  njh::StrToNumConverter::stoToNum<double>(toks[21]);
	targetDesc_ = toks[22];
	if(toks.size() > 23){
		for(const auto pos : iter::range<uint32_t>(23,toks.size())){
			targetDesc_ += " " + std::string(toks[pos]);
		}
	}
}


uint32_t HmmerDomainHitTab::zeroBasedHmmFrom() const{
	return hmmFrom_ -1;
}

uint32_t HmmerDomainHitTab::zeroBasedAlignFrom() const{
	return alignFrom_ -1;
}

uint32_t HmmerDomainHitTab::zeroBasedEnvFrom() const{
	return envFrom_ - 1;
}

uint32_t HmmerDomainHitTab::envLen() const{
	return envTo_ - zeroBasedEnvFrom();
}

uint32_t HmmerDomainHitTab::hmmLen() const{
	return hmmTo_ - zeroBasedHmmFrom();
}

uint32_t HmmerDomainHitTab::alignLen() const{
	return alignTo_ - zeroBasedAlignFrom();
}


double HmmerDomainHitTab::modelCoverage() const{
	return hmmLen()/static_cast<double>(queryLen_);
}

Bed6RecordCore HmmerDomainHitTab::genBed6_env() const {
	uint32_t start = zeroBasedEnvFrom();
	uint32_t end = envTo_;
	return genBed6(start, end);
}

Bed6RecordCore HmmerDomainHitTab::genBed6_aln() const {

	uint32_t start = zeroBasedAlignFrom();
	uint32_t end = alignTo_;
	return genBed6(start, end);

}

Bed6RecordCore HmmerDomainHitTab::genBed6(uint32_t start, uint32_t end) const{
	Bed6RecordCore out(targetName_, start, end, queryName_, end - start, '+');
	out.name_ = out.genUIDFromCoordsWithStrand();
	MetaDataInName meta;
	meta.addMeta("modelName", queryName_);
	meta.addMeta("hmmTo", hmmTo_);
	meta.addMeta("hmmFrom", zeroBasedHmmFrom());
	meta.addMeta("hmmCovered", modelCoverage());
	meta.addMeta("modelAccuracy", acc_);
	meta.addMeta("score", domainScore_);
	meta.addMeta("evalue", domain_i_evalue_);
	out.extraFields_.emplace_back(meta.createMetaName());
	return out;
}





Json::Value HmmerDomainHitTab::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["targetName_"] = njh::json::toJson(targetName_);
	ret["targetAcc_"] = njh::json::toJson(targetAcc_);
	ret["targetLen_"] = njh::json::toJson(targetLen_);
	ret["queryName_"] = njh::json::toJson(queryName_);
	ret["queryAcc_"] = njh::json::toJson(queryAcc_);
	ret["queryLen_"] = njh::json::toJson(queryLen_);
	ret["seqEvalue_"] = njh::json::toJson(seqEvalue_);
	ret["seqScore_"] = njh::json::toJson(seqScore_);
	ret["seqBias_"] = njh::json::toJson(seqBias_);
	ret["domainId_"] = njh::json::toJson(domainId_);
	ret["domainsTotals_"] = njh::json::toJson(domainsTotals_);
	ret["domain_c_evalue_"] = njh::json::toJson(domain_c_evalue_);
	ret["domain_i_evalue_"] = njh::json::toJson(domain_i_evalue_);
	ret["domainScore_"] = njh::json::toJson(domainScore_);
	ret["domainBias_"] = njh::json::toJson(domainBias_);
	ret["hmmFrom_"] = njh::json::toJson(hmmFrom_);
	ret["hmmTo_"] = njh::json::toJson(hmmTo_);
	ret["alignFrom_"] = njh::json::toJson(alignFrom_);
	ret["alignTo_"] = njh::json::toJson(alignTo_);
	ret["envFrom_"] = njh::json::toJson(envFrom_);
	ret["envTo_"] = njh::json::toJson(envTo_);
	ret["acc_"] = njh::json::toJson(acc_);
	ret["targetDesc_"] = njh::json::toJson(targetDesc_);
	return ret;
}



VecStr HmmerDomainHitTab::toDelimStrHeader (){
	return VecStr{
		"targetName",
		"targetAccession",
		"targetLen",
		"queryName",
		"queryAccession",
		"queryLen",
		"seqEvalue",
		"seqScore",
		"seqBias",
		"domainId",
		"domainsTotals",
		"domain_c_evalue",
		"domain_i_evalue",
		"domainScore",
		"domainBias",
		"hmmFrom",
		"hmmTo",
		"alignFrom",
		"alignTo",
		"envFrom",
		"envTo",
		"modelAccuracy",
		"targetDesc"
	};
}



std::string HmmerDomainHitTab::toDelimStr() const{

	return njh::conToStr(toVecStr(
			targetName_,
			targetAcc_,
			targetLen_,
			queryName_,
			queryAcc_,
			queryLen_,
			seqEvalue_,
			seqScore_,
			seqBias_,
			domainId_,
			domainsTotals_,
			domain_c_evalue_,
			domain_i_evalue_,
			domainScore_,
			domainBias_,
			hmmFrom_,
			hmmTo_,
			alignFrom_,
			alignTo_,
			envFrom_,
			envTo_,
			acc_,
			targetDesc_
		), "\t");
}




}  // namespace njhseq


