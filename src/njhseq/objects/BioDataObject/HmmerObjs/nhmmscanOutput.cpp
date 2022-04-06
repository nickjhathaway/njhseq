/*
 * nhmmscanOutput.cpp
 *
 *  Created on: Feb 19, 2022
 *      Author: nick
 */


#include "nhmmscanOutput.hpp"

#include "njhseq/IO/InputStream.hpp"
#include "njhseq/IO/OutputStream.hpp"


namespace njhseq {





Json::Value nhmmscanOutput::Hit::toJson() const {
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

	ret["hmmEdgeInfo_"] = njh::json::toJson(hmmEdgeInfo_);
	ret["aliEdgeInfo_"] = njh::json::toJson(aliEdgeInfo_);
	ret["envEdgeInfo_"] = njh::json::toJson(envEdgeInfo_);
	ret["acc_"] = njh::json::toJson(acc_);
	ret["modelAln_"] = njh::json::toJson(modelAln_);
	ret["alnAgreement_"] = njh::json::toJson(alnAgreement_);
	ret["queryAln_"] = njh::json::toJson(queryAln_);
	ret["aln_posterior_probability_"] = njh::json::toJson(aln_posterior_probability_);
	return ret;
}


Json::Value nhmmscanOutput::QueryResults::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["queryName_"] = njh::json::toJson(queryName_);
	ret["queryLen_"] = njh::json::toJson(queryLen_);
	ret["hits_"] = njh::json::toJson(hits_);
	ret["targetModels_"] = njh::json::toJson(targetModels_);
	ret["targetModNodes_"] = njh::json::toJson(targetModNodes_);
	ret["residuesSearched_"] = njh::json::toJson(residuesSearched_);
	ret["residuesPass_SSV_filter_"] = njh::json::toJson(residuesPass_SSV_filter_);
	ret["residuesPass_bias_filter_"] = njh::json::toJson(residuesPass_bias_filter_);
	ret["residuesPass_Vit_filter_"] = njh::json::toJson(residuesPass_Vit_filter_);
	ret["residuesPass_Fwd_filter_"] = njh::json::toJson(residuesPass_Fwd_filter_);
	ret["cpuRunInfo_"] = njh::json::toJson(cpuRunInfo_);
	ret["Mc_per_sec_"] = njh::json::toJson(Mc_per_sec_);
	ret["hitsFrac_"] = njh::json::toJson(hitsFrac_);
	return ret;
}

void nhmmscanOutput::QueryResults::sortHitsByEvaluesScores(std::vector<Hit> &hits) {
	njh::sort(hits, [](const Hit &hit1, const Hit &hit2) {
		if (hit1.modelEvalue_ == hit2.modelEvalue_) {
			if (hit1.modelScore_ == hit2.modelScore_) {
				return hit1.envLen() > hit2.envLen();
			} else {
				return hit1.modelScore_ > hit2.modelScore_;
			}
		} else {
			return hit1.modelEvalue_ < hit2.modelEvalue_;
		}
	});
}

std::vector<nhmmscanOutput::Hit> nhmmscanOutput::QueryResults::getNonOverlapHits(const std::vector<Hit> &hits ){
	std::vector<Hit> ret;
	for (const auto &region: hits) {
		bool overlap = false;
		for (const auto &outRegion: ret) {
			if (region.overlaps_env(outRegion, 1)) {
				overlap = true;
				break;
			}
		}
		if (!overlap) {
			ret.emplace_back(region);
		}
	}
	return ret;
}

std::vector<nhmmscanOutput::Hit> nhmmscanOutput::QueryResults::getNonOverlapHits(std::vector<Hit> &hits, const std::function<bool(const Hit&, const Hit&)> & sortFunc){
	njh::sort(hits, sortFunc);
	return getNonOverlapHits(hits);
}




Json::Value nhmmscanOutput::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["header_"] = njh::json::toJson(header_);
	ret["parameterInfo_"] = njh::json::toJson(parameterInfo_);
	ret["qResults_"] = njh::json::toJson(qResults_);
	return ret;
}

bool nhmmscanOutput::run_hmmpress_ifNeed(const bfs::path & hmmModelFnp){
	bfs::path hmmModelPressCheckFnp = njh::files::make_path(hmmModelFnp.string() + ".h3f");
	if(!bfs::exists(hmmModelPressCheckFnp) || njh::files::firstFileIsOlder(hmmModelPressCheckFnp, hmmModelFnp)){
		std::stringstream cmdSs;
		cmdSs << "hmmpress " << hmmModelFnp;
		auto cmdOutput = njh::sys::run( {cmdSs.str() });

		if (!cmdOutput.success_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch "
				 << "\n";
			ss << cmdOutput.toJson() << "\n";
			throw std::runtime_error { ss.str() };
		}
		return true;
	}
	return false;
}


nhmmscanOutput nhmmscanOutput::parseRawOutput(const bfs::path & input, const std::unordered_map<uint32_t, std::string> & seqKey){
	auto ret = parseRawOutput(input);
	ret.renameQuery(seqKey);
	return ret;
}

nhmmscanOutput nhmmscanOutput::parseRawOutput(const bfs::path & input){
	nhmmscanOutput ret;
	InputStream in(input);
	std::string line = "";

	{
		//first read header
		bool readingParameters = false;
		while(njh::files::crossPlatGetline(in, line) && njh::beginsWith(line, "#")){
			if(njh::beginsWith(line, "# - - - - -")){
				readingParameters = !readingParameters;
			}else	if(readingParameters){
				auto toks = njh::tokenizeString(line, ":");
				if(2 != toks.size()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "expected two tokens, not: " << toks.size()<< "\n";
					throw std::runtime_error{ss.str()};
				}
				//remove trailing or heading whitespace
				njh::lstrip(toks[0], '#');
				njh::trim(toks[0]);
				njh::trim(toks[1]);
				ret.parameterInfo_[toks[0]] = toks[1];
			}else{
				ret.header_.emplace_back(line);
			}
		}
	}
	auto checkLineTokNumber = [](uint32_t expectedTokNumber, uint32_t observedTokNumber,const std::string & line, const std::string & parsedLineName, const std::string & funcname){
		if(expectedTokNumber != observedTokNumber){
			std::stringstream ss;
			ss << funcname << ", error " << "error in parsing " << parsedLineName<< ": " << line << "\n";
			ss << "Expected " << expectedTokNumber << " toks, have: " << observedTokNumber	<< "\n";
			throw std::runtime_error{ss.str()};
		}
	};

	{
		//read through each query hits
		while(njh::files::crossPlatGetline(in, line)){
			//check start of query
			if(njh::beginsWith(line, "Query:")){
				//process for Query info
				nhmmscanOutput::QueryResults queryRes;
				{
					auto queryLineToks = njh::tokenizeString(line, ":");
					checkLineTokNumber(2, queryLineToks.size(), line, "query line", __PRETTY_FUNCTION__);
					auto lenBegin = queryLineToks[1].rfind("[L=");
					auto lenEnd = queryLineToks[1].rfind("]");
					if(std::string::npos == lenBegin || std::string::npos == lenEnd || lenEnd <= (lenBegin + 3)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "in parsing query line for length: " << line << "\n";
						ss << "lenBegin: " << lenBegin << "\n";
						ss << "lenEnd: " << lenEnd << "\n";
						throw std::runtime_error{ss.str()};
					}
					queryRes.queryName_ = queryLineToks[1].substr(0, lenBegin);
					njh::trim(queryRes.queryName_);
					queryRes.queryLen_ = njh::StrToNumConverter::stoToNum<uint64_t>(queryLineToks[1].substr(lenBegin + 3, lenEnd - (lenBegin + 3)));
				}


				std::string subLine;

				//parse until get
				bool endQueryParse = false;
				while(njh::files::crossPlatGetline(in, subLine) && !endQueryParse){

					if(njh::beginsWith(subLine, ">> ")){
						//process subLine for model name;
						nhmmscanOutput::Hit hit;
						hit.queryName_ = queryRes.queryName_;
						hit.targetName_ = subLine.substr(3);
						njh::trim(hit.targetName_);
						std::string hitline = "";

						//hit header
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						//under header
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						//hit values
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						auto hitValuesToks = tokenizeString(hitline, "whitespace");
						checkLineTokNumber(15, hitValuesToks.size(), hitline, "hit values line", __PRETTY_FUNCTION__);
						//[1]=score, [2]=bias,, [3]=Evalue, [4]=hmmfrom, [5]=hmm to, [6]=hmmflank, [7]=alifrom, [8]=ali to, [9]=aliflank, [10]=envfrom, [11]=env to [12]=envflnak, [13]=mod len, [14]=acc
						hit.modelScore_ = njh::StrToNumConverter::stoToNum<double>(hitValuesToks[1]);
						hit.modelBias_ = njh::StrToNumConverter::stoToNum<double>(hitValuesToks[2]);
						hit.modelEvalue_ = njh::StrToNumConverter::stoToNum<double>(hitValuesToks[3]);
						//hmm
						hit.hmmFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[4]);
						hit.hmmTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[5]);
						hit.hmmEdgeInfo_ = hitValuesToks[6];
						//ali
						hit.alignFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[7]);
						hit.alignTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[8]);
						hit.aliEdgeInfo_ = hitValuesToks[9];
						hit.strand_ = hit.isReverseStrand() ? '-' : '+';
						//env
						hit.envFrom_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[10]);
						hit.envTo_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[11]);
						hit.envEdgeInfo_ = hitValuesToks[12];

						hit.modelLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitValuesToks[13]);
						hit.acc_ = njh::StrToNumConverter::stoToNum<double>(hitValuesToks[14]);

						//blank line
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						//Alignment header
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						//score
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						//model aln
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						{
							auto modelAlnToks = tokenizeString(hitline, "whitespace");
							checkLineTokNumber(4, modelAlnToks.size(), hitline, "model aln line", __PRETTY_FUNCTION__);
							if(hit.targetName_ != modelAlnToks[0]){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing model aln line: " << hitline << "\n";
								ss << "Expected modelAlnToks[0] tp be " << hit.targetName_ << ", not: " << modelAlnToks[0]	<< "\n";
								throw std::runtime_error{ss.str()};
							}
							if(njh::pasteAsStr(hit.hmmFrom_) != modelAlnToks[1]){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing model aln line: " << hitline << "\n";
								ss << "Expected modelAlnToks[1] tp be " << hit.hmmFrom_ << ", not: " << modelAlnToks[1]	<< "\n";
								throw std::runtime_error{ss.str()};
							}
							if(njh::pasteAsStr(hit.hmmTo_) != modelAlnToks[3]){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing model aln line: " << hitline << "\n";
								ss << "Expected modelAlnToks[3] tp be " << hit.hmmTo_ << ", not: " << modelAlnToks[3]	<< "\n";
								throw std::runtime_error{ss.str()};
							}
							hit.modelAln_ = modelAlnToks[2];
						}
						//aln agreement
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						{
							hit.alnAgreement_ = hitline;
						}
						//query aln
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						{
							auto queryAlnToks = tokenizeString(hitline, "whitespace");
							checkLineTokNumber(4, queryAlnToks.size(), hitline, "query aln line", __PRETTY_FUNCTION__);
							if(queryRes.queryName_ != queryAlnToks[0]){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing model aln line: " << hitline << "\n";
								ss << "Expected queryAlnToks[0] tp be " << queryRes.queryName_ << ", not: " << queryAlnToks[0]	<< "\n";
								throw std::runtime_error{ss.str()};
							}
							if(njh::pasteAsStr(hit.alignFrom_) != queryAlnToks[1]){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing model aln line: " << hitline << "\n";
								ss << "Expected modelAlnToks[1] tp be " << hit.alignFrom_ << ", not: " << queryAlnToks[1]	<< "\n";
								throw std::runtime_error{ss.str()};
							}
							if(njh::pasteAsStr(hit.alignTo_) != queryAlnToks[3]){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing model aln line: " << hitline << "\n";
								ss << "Expected modelAlnToks[3] tp be " << hit.alignTo_ << ", not: " << queryAlnToks[3]	<< "\n";
								throw std::runtime_error{ss.str()};
							}
							hit.queryAln_ = queryAlnToks[2];
						}
						//posterior prob
						njh::files::crossPlatGetline(in, hitline);
						njh::trim(hitline);
						auto ppToks = njh::tokenizeString(hitline, "whitespace");
						checkLineTokNumber(2, ppToks.size(), hitline, "posterior probability line", __PRETTY_FUNCTION__);

						if(ppToks[1] != "PP"){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "error in parsing posterior probability line: " << hitline << "\n";
							ss << "Expected ppToks[1] tp be PP, not: " << ppToks[1]	<< "\n";
							throw std::runtime_error{ss.str()};
						}
						hit.aln_posterior_probability_ = ppToks[0];
						queryRes.hits_.emplace_back(hit);
					} else if (njh::beginsWith(subLine, "Internal pipeline statistics summary")){
						//reading summary stats for query
//						uint64_t targetModNodes_{std::numeric_limits<uint64_t>::max()};
//						uint64_t residuesSearched_{std::numeric_limits<uint64_t>::max()};
//						uint64_t residuesPass_SSV_filter_{std::numeric_limits<uint64_t>::max()};
//						uint64_t residuesPass_bias_filter_{std::numeric_limits<uint64_t>::max()};
//						uint64_t residuesPass_Vit_filter_{std::numeric_limits<uint64_t>::max()};
//						uint64_t residuesPass_Fwd_filter_{std::numeric_limits<uint64_t>::max()};
//
//						std::string cpuRunInfo_;
//						double Mc_per_sec_{std::numeric_limits<double>::max()};
//						double hitsFrac_{std::numeric_limits<double>::max()};

						std::string summaryLine = "";
						//all dahses
						njh::files::crossPlatGetline(in, summaryLine);
						//query sequence info
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto queryInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, queryInfoToks.size(), summaryLine, "query info line", __PRETTY_FUNCTION__);
							auto parStart = queryInfoToks[1].find("(");
							auto end = queryInfoToks[1].find(" residues searched");
							if(end <= parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parStart: " << parStart << " should be less than end: " << end << "\n";
								ss << "for line: " << queryInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							queryRes.residuesSearched_ = njh::StrToNumConverter::stoToNum<uint32_t>(queryInfoToks[1].substr(parStart + 1, end -1 - parStart));
						}
						//target model
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto modelsInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, modelsInfoToks.size(), summaryLine, "model info line", __PRETTY_FUNCTION__);
							auto parStart = modelsInfoToks[1].find("(");
							auto end = modelsInfoToks[1].find(" nodes");
							if(end <= parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parStart: " << parStart << " should be less than end: " << end << "\n";
								ss << "for line: " << modelsInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							std::string modelNum = modelsInfoToks[1].substr(0, parStart);
							njh::trim(modelNum);
							queryRes.targetModels_ = njh::StrToNumConverter::stoToNum<uint32_t>(modelNum);
							queryRes.targetModNodes_ = njh::StrToNumConverter::stoToNum<uint32_t>(modelsInfoToks[1].substr(parStart + 1, end -1 - parStart));
						}
						//SSV filter
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto ssvInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, ssvInfoToks.size(), summaryLine, "ssv info line", __PRETTY_FUNCTION__);
							auto parStart = ssvInfoToks[1].find("(");
							if(std::string::npos == parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parStart couldn't be found\n";
								ss << "for line: " << ssvInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							std::string ssvNum = ssvInfoToks[1].substr(0, parStart);
							njh::trim(ssvNum);
							queryRes.residuesPass_SSV_filter_ = njh::StrToNumConverter::stoToNum<uint32_t>(ssvNum);
						}
						//bias filter
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto biasInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, biasInfoToks.size(), summaryLine, "bias info line", __PRETTY_FUNCTION__);
							auto parStart = biasInfoToks[1].find("(");
							if(std::string::npos == parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parStart couldn't be found\n";
								ss << "for line: " << biasInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							std::string biasNum = biasInfoToks[1].substr(0, parStart);
							njh::trim(biasNum);
							queryRes.residuesPass_bias_filter_ = njh::StrToNumConverter::stoToNum<uint32_t>(biasNum);
						}
						//Vit filter
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto vitInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, vitInfoToks.size(), summaryLine, "vit info line", __PRETTY_FUNCTION__);
							auto parStart = vitInfoToks[1].find("(");
							if(std::string::npos == parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parStart couldn't be found\n";
								ss << "for line: " << vitInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							std::string vitNum = vitInfoToks[1].substr(0, parStart);
							njh::trim(vitNum);
							queryRes.residuesPass_Vit_filter_ = njh::StrToNumConverter::stoToNum<uint32_t>(vitNum);
						}
						//Fwd filter
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto fwdInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, fwdInfoToks.size(), summaryLine, "fwd info line", __PRETTY_FUNCTION__);
							auto parStart = fwdInfoToks[1].find("(");
							if(std::string::npos == parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parStart couldn't be found\n";
								ss << "for line: " << fwdInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							std::string fwdNum = fwdInfoToks[1].substr(0, parStart);
							njh::trim(fwdNum);
							queryRes.residuesPass_Fwd_filter_ = njh::StrToNumConverter::stoToNum<uint32_t>(fwdNum);
						}
						//total hits
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto hitsInfoToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, hitsInfoToks.size(), summaryLine, "hits info line", __PRETTY_FUNCTION__);
							auto parStart = hitsInfoToks[1].find("(");
							auto parEnd = hitsInfoToks[1].find(")");
							if(parEnd <= parStart){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "parEnd: " << parEnd << " shouldn't be less than parStart: " << parStart << "\n";
								ss << "for line: " << hitsInfoToks[1] << "\n";
								throw std::runtime_error{ss.str()};
							}
							queryRes.hitsFrac_ = njh::StrToNumConverter::stoToNum<uint32_t>(hitsInfoToks[1].substr(parStart + 1, parEnd -1 - parStart));
						}
						//CPU run info
						njh::files::crossPlatGetline(in, summaryLine);
						queryRes.cpuRunInfo_ = summaryLine;
						//Mc/sec
						njh::files::crossPlatGetline(in, summaryLine);
						{
							auto speedToks = njh::tokenizeString(summaryLine, ":");
							checkLineTokNumber(2, speedToks.size(), summaryLine, "speed line", __PRETTY_FUNCTION__);
							njh::trim(speedToks[1]);
							queryRes.Mc_per_sec_ = njh::StrToNumConverter::stoToNum<double>(speedToks[1]);
						}
						//end
						njh::files::crossPlatGetline(in, summaryLine);
						if("//" != summaryLine){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "error, last line in query should be //, not: " << summaryLine<< "\n";
							throw std::runtime_error{ss.str()};
						}else{
							endQueryParse = true;
							ret.qResults_.emplace_back(queryRes);
							break;
						}
					}
				}
			}
		}
	}

	for(auto & query : ret.qResults_){
		QueryResults::sortHitsByEvaluesScores(query.hits_);
	}
	return ret;
}



uint32_t nhmmscanOutput::Hit::sumOfPosteriorProbability() const{
	uint32_t ret = 0;
	for(const auto c : aln_posterior_probability_){
		if('*' == c){
			ret += 10;
		}else if ('.' == c){
			ret += 0;
		}else{
			ret += njh::StrToNumConverter::stoToNum<uint32_t>(std::string(1, c));
		}
	}
	return ret;
}

double nhmmscanOutput::Hit::averagePP() const{
	return sumOfPosteriorProbability()/static_cast<double>(aln_posterior_probability_.size())	;
}



double nhmmscanOutput::Hit::percentPerfectHit() const{
	uint32_t perfect = njh::count_if(aln_posterior_probability_, [](char c){
		return '*' == c;
	});
	return static_cast<double>(perfect)/aln_posterior_probability_.size();
}

double nhmmscanOutput::Hit::percentGappedHit() const{
	uint32_t aln_modle_gapped = njh::count_if(aln_posterior_probability_, [](char c){
		return '.' == c;
	});
	uint32_t aln_query_gapped = njh::count_if(queryAln_, [](char c){
		return '-' == c;
	});
	return (aln_modle_gapped + aln_query_gapped)/static_cast<double>(aln_posterior_probability_.size());
}

bool nhmmscanOutput::Hit::overlaps_env(const Hit & otherHit, uint32_t minOverlap)const{
	return Bed3RecordCore::getOverlapLen(queryName_, env0BasedPlusStrandStart(), env0BasedPlusStrandEnd(),
																			 otherHit.queryName_, otherHit.env0BasedPlusStrandStart(), otherHit.env0BasedPlusStrandEnd()) >= minOverlap;
}

VecStr nhmmscanOutput::Hit::getOutputDetHeader(){
	return VecStr {"model"
					,"hmm_From"
					,"hmm_to"
					,"hmm_edges"
					,"model_len"

					,"aln_from"
					,"aln_to"
					,"aln_len"
					,"aln_edge"

					,"env_from"
					,"env_to"
					,"env_len"
					,"env_edge"

					,"strand"
					,"evalue"
					,"score"
					,"bias"
					,"acc"
					,"scoreOverLen"
					,"percentGappedHit"
					,"percentPerfectHit"
					,"averagePP"};
}

VecStr nhmmscanOutput::Hit::getOutputDet() const{
	return toVecStr(targetName_
														,hmmFrom_
														,hmmTo_
														,hmmEdgeInfo_
														,modelLen_

														,alignFrom_
														,alignTo_
														,alignLen()
														,aliEdgeInfo_

														,envFrom_
														,envTo_
														,envLen()
														,envEdgeInfo_

														,strand_
														,modelEvalue_
														,modelScore_
														,modelBias_
														,acc_
														,modelScore_/envLen()
														,percentGappedHit()
														,percentPerfectHit()
														,averagePP());
}


void nhmmscanOutput::renameQuery(const std::unordered_map<uint32_t, std::string> & seqKey){
	for(auto & query : qResults_){
		query.queryName_ = seqKey.at(njh::StrToNumConverter::stoToNum<uint32_t>(query.queryName_));
		for(auto & hit : query.hits_){
			hit.queryName_ = query.queryName_;
		}
	}
}

void nhmmscanOutput::outputCustomHitsTable(const OutOptions &outOpts) const {
	OutputStream out(outOpts);

	out << "query"
			<< "\t" << "queryLen";
	out << "\t" << njh::conToStr(Hit::getOutputDetHeader(), "\t");
	out << std::endl;

	for (const auto &res: qResults_) {
		for (const auto &hit: res.hits_) {
			out << res.queryName_
					<< "\t" << res.queryLen_;
			out << "\t" << njh::conToStr(hit.getOutputDet(), "\t");
			out << std::endl;
		}
	}
}


nhmmscanOutput::PostProcessHitsRes nhmmscanOutput::postProcessHits(const PostProcessHitsPars & pars) {
	PostProcessHitsRes ret;
	for(auto & query : qResults_){
		for(const auto & hit : query.hits_){
			auto loc = hit.genBed6_env();
			uint32_t hmmFromBack = hit.modelLen_ - hit.hmmTo_;
			uint32_t hmmFromFront = hit.zeroBasedHmmFrom();
			if(
							hmmFromBack <= pars.hmmStartFilter &&
							hmmFromFront <= pars.hmmStartFilter &&
							loc.length() >= pars.minLength &&
							hit.acc_ >= pars.hardAccCutOff &&
							hit.modelScore_ >= pars.hardScoreCutOff &&
							( hit.acc_ >=  pars.accCutOff || hit.modelScore_ >=  pars.scoreCutOff)){
				ret.filteredHitsByQuery_[hit.queryName_].emplace_back(hit);
			}
		}
	}
	for(auto & filteredHits : ret.filteredHitsByQuery_){
		nhmmscanOutput::QueryResults::sortHitsByEvaluesScores(filteredHits.second);
		ret.filteredNonOverlapHitsByQuery_[filteredHits.first] = nhmmscanOutput::QueryResults::getNonOverlapHits(filteredHits.second);
	}
	return ret;
}

void nhmmscanOutput::writeInfoFiles(const PostProcessHitsRes & postProcessResults, const bfs::path & outputDir) const {
	OutputStream allHitsBedOut(njh::files::make_path(outputDir, "all_hits.bed"));
	for(auto & query : qResults_){
		for(const auto & hit : query.hits_){
			auto loc = hit.genBed6_env();
			allHitsBedOut << loc.toDelimStrWithExtra() << std::endl;
		}
	}
	OutOptions hitTableOutOpts = njh::files::make_path(outputDir, "nhmmscan_hits_table.tab.txt");
	outputCustomHitsTable(hitTableOutOpts);

	OutputStream hitFilteredTableOut(njh::files::make_path(outputDir, "nhmmscan_hits_filtered_table.tab.txt"));
	OutputStream hitFilteredBedOut(njh::files::make_path(outputDir, "nhmmscan_hits_filtered.bed"));

	OutputStream hitNonOverlapFilteredTableOut(njh::files::make_path(outputDir, "nhmmscan_hits_nonOverlap_filtered_table.tab.txt"));
	OutputStream hitNonOverlapFilteredBedOut(njh::files::make_path(outputDir, "nhmmscan_hits_nonOverlap_filtered.bed"));

	hitFilteredTableOut << "queryName\t" << njh::conToStr(nhmmscanOutput::Hit::getOutputDetHeader(), "\t") << std::endl;
	hitNonOverlapFilteredTableOut << "queryName\t" << njh::conToStr(nhmmscanOutput::Hit::getOutputDetHeader(), "\t") << std::endl;
	for(const auto & filteredHits : postProcessResults.filteredHitsByQuery_){
		for(const auto & hit : filteredHits.second){
			hitFilteredTableOut << filteredHits.first << "\t" << njh::conToStr(hit.getOutputDet(), "\t") << std::endl;
			hitFilteredBedOut << hit.genBed6_env().toDelimStrWithExtra() << std::endl;
		}
		for(const auto & hit : postProcessResults.filteredNonOverlapHitsByQuery_.at(filteredHits.first)){
			hitNonOverlapFilteredTableOut << filteredHits.first << "\t" << njh::conToStr(hit.getOutputDet(), "\t") << std::endl;
			hitNonOverlapFilteredBedOut << hit.genBed6_env().toDelimStrWithExtra() << std::endl;
		}
	}
}



}  // namespace njhseq


