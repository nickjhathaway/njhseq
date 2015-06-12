/*
 * MidDeterminator.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */

#include "MidDeterminator.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {

midPos::midPos() :
		midName_("unrecognized"), midPos_(midPos::npos), barcodeSize_(0) {
}
midPos::midPos(const std::string & midName, uint64_t midPos, uint64_t barcodeSize) :
		midName_(midName), midPos_(midPos), barcodeSize_(barcodeSize) {
}

Json::Value midPos::toJson()const{
	Json::Value ret;
	ret["midName_"] = bib::json::toJson(midName_);
	ret["midPos_"] = bib::json::toJson(midPos_);
	ret["barcodeSize_"] = bib::json::toJson(barcodeSize_);
	return ret;
}

midPos::operator bool() const {
	return npos != midPos_;
}

MidDeterminator::MidDeterminator(const table & mids) {
	if (!bib::in(std::string("id"), mids.columnNames_)
			|| !bib::in(std::string("barcode"), mids.columnNames_)) {
		throw std::runtime_error {
				"Error in creating MidDeterminator, need to have at "
						"least two columns, id and barcode, only have "
						+ bib::conToStr(mids.columnNames_, ",") };
	}
	bool addRComp = true;;
	for (const auto & row : mids.content_) {
		addBarcode(row[mids.getColPos("id")],row[mids.getColPos("barcode")], addRComp);
	}
}

bool MidDeterminator::containsMidByName(const std::string & name)const{
	return bib::in(name, mids_);
}
bool MidDeterminator::containsMidByBarcode(const std::string & barcode, bool checkComp)const{
	for(const auto & mid : mids_){
		if(mid.second.motifOriginal_ == barcode){
			return true;
		}
	}
	if(checkComp){
		for(const auto & mid : compMids_){
			if(mid.second.motifOriginal_ == seqUtil::reverseComplement(barcode, "DNA")){
				return true;
			}
		}
	}
	return false;
}

std::string MidDeterminator::getMidName(const std::string & barcode, bool checkComp)const{
	for(const auto & mid : mids_){
		if(mid.second.motifOriginal_ == barcode){
			return mid.first;
		}
	}
	if(checkComp){
		for(const auto & mid : compMids_){
			if(mid.second.motifOriginal_ == seqUtil::reverseComplement(barcode, "DNA")){
				return mid.first;
			}
		}
	}
	return "no_name_for_barcode:" + barcode;
}

void MidDeterminator::addBarcode(const std::string & name, const std::string & barcode,
		bool addToComp) {
	if (containsMidByName(name)) {
		std::stringstream ss;
		ss << "Error, MidDeterminator already contains mid by name: " << name
				<< "\n";
		ss << "original barcode: " << mids_.at(name).motifOriginal_
				<< ", adding barcode: " << barcode << "\n";
		throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
	}
	if (containsMidByBarcode(barcode, addToComp)) {
		std::stringstream ss;
		ss << "Error, MidDeterminator already contains mid by barcode: "
				<< barcode << "\n";
		ss << "original name: " << getMidName(barcode, addToComp)
				<< ", adding name: " << name << "\n";
		throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
	}
	mids_.emplace(name, motif { barcode });
	if (addToComp) {
		compMids_.emplace(name,
				motif { seqUtil::reverseComplement(barcode, "DNA") });
	}
}




midPos MidDeterminator::determineMidSimple(const std::string & seq) {
	for (auto & mid : mids_) {
		if (mid.second.passMotifParameter(
				seq.substr(0, mid.second.motifOriginal_.size()),
				mid.second.motifOriginal_.size() - 0)) {
			return {mid.first,0, mid.second.motifOriginal_.size()};
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidSimpleBack(const std::string & seq) {
	for (auto & mid : mids_) {
		if(seq.size() < mid.second.motifOriginal_.size()){
			continue;
		}
		if (mid.second.passMotifParameter(
				seq.substr(seq.size() - mid.second.motifOriginal_.size(), mid.second.motifOriginal_.size()),
				mid.second.motifOriginal_.size() - 0)) {
			return {mid.first, seq.size() - mid.second.motifOriginal_.size(), mid.second.motifOriginal_.size()};
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidSimpleComp(const std::string & seq) {
	for (auto & mid : compMids_) {
		if(seq.size() < mid.second.motifOriginal_.size()){
			continue;
		}
		if (mid.second.passMotifParameter(
				seq.substr(seq.size() - mid.second.motifOriginal_.size(), mid.second.motifOriginal_.size()),
				mid.second.motifOriginal_.size() - 0)) {
			return {mid.first,seq.size() - mid.second.motifOriginal_.size(), mid.second.motifOriginal_.size()};
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidSimpleCompFront(const std::string & seq) {
	for (auto & mid : compMids_) {
		if (mid.second.passMotifParameter(
				seq.substr(0, mid.second.motifOriginal_.size()), mid.second.motifOriginal_.size() - 0)) {
			return {mid.first,0, mid.second.motifOriginal_.size()};
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidPosVarStart(const std::string & seq, uint32_t varStop) {
	for (auto & mid : mids_) {
		varStop = std::min<uint32_t>(varStop, seq.size());
		auto positions = mid.second.findPositionsFull(seq, 0, 0,varStop);
		if (!positions.empty() && positions.front() <= varStop) {
			return midPos { mid.first, positions.front(),mid.second.motifOriginal_.size() };
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidPosVarStartBack(const std::string & seq, uint32_t varStop) {
	for (auto & mid : mids_) {
		if(varStop > seq.size()){
			varStop = seq.size();
		}
		auto positions = mid.second.findPositionsFull(seq, 0, seq.size() - varStop, seq.size());
		if (!positions.empty() && positions.front() >= (seq.size() - varStop)) {
			return midPos { mid.first, positions.front(), mid.second.motifOriginal_.size() };
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidPosVarStartComp(const std::string & seq, uint32_t varStop) {
	for (auto & mid : compMids_) {
		if(varStop > seq.size()){
			varStop = seq.size();
		}
		auto positions = mid.second.findPositionsFull(seq, 0, seq.size() - varStop, seq.size());
		if (!positions.empty() && positions.front() >= (seq.size() - varStop)) {
			return midPos { mid.first, positions.front(), mid.second.motifOriginal_.size() };
		}
	}
	return midPos {};
}

midPos MidDeterminator::determineMidPosVarStartCompFront(const std::string & seq, uint32_t varStop) {
	for (auto & mid : compMids_) {
		varStop = std::min<uint32_t>(varStop, seq.size());
		auto positions = mid.second.findPositionsFull(seq, 0, 0,varStop);
		if (!positions.empty() && positions.front() <= varStop) {
			return midPos { mid.first, positions.front(),mid.second.motifOriginal_.size() };
		}
	}
	return midPos {};
}

midPos MidDeterminator::fullDetermine(seqInfo & info, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds ){
  midPos ret;
  if(variableStart){
  	auto detMid = determineMidPosVarStart(info.seq_, variableStop);
  	if(detMid){
  		ret = detMid;
  	}else if(checkComplement){
  		auto compDetMid = determineMidPosVarStartComp(info.seq_, variableStop);
  		if(compDetMid){
  			compDetMid.midPos_ = info.seq_.size() - compDetMid.midPos_ - compDetMid.barcodeSize_ ;
  			info.reverseComplementRead(true, true);
  			ret = compDetMid;
  		}
  	}
  }else{
  	auto detMid = determineMidSimple(info.seq_);
  	if(detMid){
  		ret = detMid;
  	}else if(checkComplement){
  		auto compDetMid = determineMidSimpleComp(info.seq_);
  		if(compDetMid){
  			compDetMid.midPos_ = info.seq_.size() - compDetMid.midPos_ - compDetMid.barcodeSize_ ;
  			info.reverseComplementRead(true, true);
  			ret = compDetMid;
  		}
  	}
  }
  if(!ret){
    if(barcodesBothEnds){
      if(variableStart){
      	auto detMid = determineMidPosVarStartBack(info.seq_, variableStop);
      	if(detMid){
      		ret = detMid;
      	}else if(checkComplement){
      		auto compDetMid = determineMidPosVarStartCompFront(info.seq_, variableStop);
      		if(compDetMid){
      			compDetMid.midPos_ = info.seq_.size() - compDetMid.midPos_ - compDetMid.barcodeSize_ ;
      			info.reverseComplementRead(true, true);
      			ret = compDetMid;
      		}
      	}
      }else{
      	auto detMid = determineMidSimpleBack(info.seq_);
      	if(detMid){
      		ret = detMid;
      	}else if(checkComplement){
      		auto compDetMid = determineMidSimpleCompFront(info.seq_);
      		if(compDetMid){
      			compDetMid.midPos_ = info.seq_.size() - compDetMid.midPos_ - compDetMid.barcodeSize_ ;
      			info.reverseComplementRead(true, true);
      			ret = compDetMid;
      		}
      	}
      }
    }

    if(ret){
    	//remove barcodes
    	info.trimBack(ret.midPos_);
    	//remove the same size other side, should add a better check for the other side barcode
    	info.trimFront(ret.barcodeSize_);
    }
  }else{
    //remove barcodes
    info.trimFront(ret.midPos_ + ret.barcodeSize_);
    if(barcodesBothEnds){
    	info.trimBack(len(info) - ret.barcodeSize_);
    }
  }
  return ret;
}


} /* namespace bibseq */
