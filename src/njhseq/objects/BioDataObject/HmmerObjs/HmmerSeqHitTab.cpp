/*
 * HmmerTableDomainHit.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */



#include "HmmerSeqHitTab.hpp"
#include "njhseq/utils/vectorUtils.hpp"
#include "njhseq/objects/Meta/MetaDataInName.hpp"


namespace njhseq {

HmmerSeqHitTab::HmmerSeqHitTab(){

}

// target name            accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  modlen strand   E-value  score  bias  description of target

HmmerSeqHitTab::HmmerSeqHitTab(const std::string & line) {
	auto toks = tokenizeString(line, "whitespace");
	if (toks.size() < 16) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error there should be at least 16 items not "
				<< toks.size() << "\n";
		ss << line << "\n";
		for(const auto pos : iter::range(toks.size())){
			ss << pos << ": " << toks[pos] << "\n";
		}

		throw std::runtime_error { ss.str() };
	}
	targetName_ = toks[0];
	targetAcc_ = toks[1];
	queryName_ = toks[2];

	queryAcc_ = toks[3];

	hmmFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[4]);
	hmmTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);

	alignFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[6]);
	alignTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[7]);

	envFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[8]);
	envTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[9]);
	modelLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[10]);

	strand_ = toks[11].at(0);
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "toks[12]: " << toks[12] << std::endl;
	modelEvalue_ = njh::StrToNumConverter::stoToNum<double>(toks[12]);
	modelScore_ = njh::StrToNumConverter::stoToNum<double>(toks[13]);
//	if(0 == modelEvalue_ && modelScore_ > 900){
//		modelEvalue_ = std::numeric_limits<double>::min() * 10;
//	}
	modelBias_ = njh::StrToNumConverter::stoToNum<double>(toks[14]);

	targetDesc_ = toks[15];
	if(toks.size() > 16){
		for(const auto pos : iter::range<uint32_t>(16, toks.size())){
			targetDesc_ += " " + std::string(toks[pos]);
		}
	}
}

//
//bool HmmerSeqHitTab::isReverseStrand() const {
//	return '-' == strand_;
//}

double HmmerSeqHitTab::modelCoverage() const{
	double modelLenCovered = uAbsdiff(hmmFrom_, hmmTo_) + 1; //positions are 1 based
	return modelLenCovered/modelLen_;
}

double HmmerSeqHitTab::queryCoverageAln(uint32_t queryLen) const{
	return alignLen()/static_cast<double>(queryLen);
}
double HmmerSeqHitTab::queryCoverageEnv(uint32_t queryLen) const{
	return envLen()/static_cast<double>(queryLen);
}





//
//uint32_t HmmerSeqHitTab::env0BasedPlusStrandStart() const{
//	uint32_t start = !isReverseStrand() ? envFrom_ - 1 : envTo_ - 1;
//	return start;
//}
//
//uint32_t HmmerSeqHitTab::env0BasedPlusStrandEnd() const{
//	uint32_t end = !isReverseStrand() ? envTo_ : envFrom_;
//	return end;
//}
//
//uint32_t HmmerSeqHitTab::align0BasedPlusStrandStart() const{
//	uint32_t start = !isReverseStrand() ? alignFrom_ - 1 : alignTo_ - 1;
//	return start;
//}
//
//uint32_t HmmerSeqHitTab::align0BasedPlusStrandEnd() const{
//	uint32_t end = !isReverseStrand() ? alignTo_ : alignFrom_;
//	return end;
//}
//
//
//uint32_t HmmerSeqHitTab::zeroBasedHmmFrom() const{
//	//per file specs the positions are 1-based
//	return hmmFrom_ - 1;
//}
//
//uint32_t HmmerSeqHitTab::envLen() const{
//	return env0BasedPlusStrandEnd() - env0BasedPlusStrandStart();
//}
//uint32_t HmmerSeqHitTab::hmmLen() const{
//	return hmmTo_ - zeroBasedHmmFrom();
//}
//uint32_t HmmerSeqHitTab::alignLen() const{
//	return align0BasedPlusStrandEnd() - align0BasedPlusStrandStart();
//}



Json::Value HmmerSeqHitTab::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["targetName_"] = njh::json::toJson(targetName_);
	ret["targetAcc_"] = njh::json::toJson(targetAcc_);
	ret["queryName_"] = njh::json::toJson(queryName_);
	ret["queryAcc_"] = njh::json::toJson(queryAcc_);

	ret["hmmFrom_"] = njh::json::toJson(hmmFrom_);
	ret["hmmTo_"] = njh::json::toJson(hmmTo_);
	ret["alignFrom_"] = njh::json::toJson(alignFrom_);
	ret["alignTo_"] = njh::json::toJson(alignTo_);
	ret["envFrom_"] = njh::json::toJson(envFrom_);
	ret["envTo_"] = njh::json::toJson(envTo_);
	ret["modelLen_"] = njh::json::toJson(modelLen_);
	ret["strand_"] = njh::json::toJson(strand_);


	ret["modelEvalue_"] = njh::json::toJson(modelEvalue_);
	ret["modelScore_"] = njh::json::toJson(modelScore_);
	ret["modelBias_"] = njh::json::toJson(modelBias_);
	ret["targetDesc_"] = njh::json::toJson(targetDesc_);
	return ret;
}





VecStr HmmerSeqHitTab::toDelimStrHeader (){
	return VecStr{
		"targetName",
		"targetAcc",
		"queryName",
		"queryAcc",
		"hmmFrom",
		"hmmTo",
		"alignFrom",
		"alignTo",
		"envFrom",
		"envTo",

		"modelLen",
		"strand",

		"modelEvalue",
		"modelScore",
		"modelBias",
		"targetDesc"
	};
}



std::string HmmerSeqHitTab::toDelimStr() const{

	return njh::conToStr(toVecStr(
			targetName_,
			targetAcc_,
			queryName_,
			queryAcc_,
			hmmFrom_,
			hmmTo_,
			alignFrom_,
			alignTo_,
			envFrom_,
			envTo_,
			modelLen_,
			strand_,

			modelEvalue_,
			modelScore_,
			modelBias_,
			targetDesc_
		), "\t");
}



Bed6RecordCore HmmerSeqHitTab::genBed6(uint32_t start, uint32_t end) const{

	Bed6RecordCore out(queryName_, start, end, targetName_, end - start, strand_);
	out.name_ = out.genUIDFromCoordsWithStrand();
	MetaDataInName meta;
	meta.addMeta("modelName", targetName_);
	meta.addMeta("targetDesc", targetDesc_);
	meta.addMeta("hmmTo", hmmTo_);
	meta.addMeta("hmmFrom", hmmFrom_);
	meta.addMeta("hmmCovered", modelCoverage());
	meta.addMeta("score", modelScore_);
	meta.addMeta("evalue", modelEvalue_);
	out.extraFields_.emplace_back(meta.createMetaName());
	return out;
}



Bed6RecordCore HmmerSeqHitTab::genBed6_aln() const {
	uint32_t start = align0BasedPlusStrandStart();
	uint32_t end = align0BasedPlusStrandEnd();
	return genBed6(start, end);
}

Bed6RecordCore HmmerSeqHitTab::genBed6_env() const {
	uint32_t start = env0BasedPlusStrandStart();
	uint32_t end = env0BasedPlusStrandEnd();
	return genBed6(start, end);
}

}  // namespace njhseq



