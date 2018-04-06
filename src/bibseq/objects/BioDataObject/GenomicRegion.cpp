/*
 * GenomicRegion.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include "GenomicRegion.hpp"
#include "bibseq/objects/BioDataObject/reading.hpp"
#include "bibseq/helpers/seqUtil.hpp"


namespace bibseq {

GenomicRegion::GenomicRegion() :
		uid_(""), chrom_(""), start_(0), end_(0) {
}

GenomicRegion::GenomicRegion(const std::string & uid, const std::string & chrom,
		size_t start, size_t end, bool reverseSrand) :
		uid_(uid), alternateUids_ { uid }, chrom_(chrom), start_(start), end_(end), reverseSrand_(
				reverseSrand) {
	if (end_ <= start_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", end: " << end_ << " shouldn't be less than start: " << start_ << " for region: " << uid_ << "\n";
		throw std::runtime_error{ss.str()};
	}
}

GenomicRegion::GenomicRegion(const Bed6RecordCore & bed) :
		GenomicRegion(bed.name_, bed.chrom_, bed.chromStart_, bed.chromEnd_,
				'-' == bed.strand_) {

}

GenomicRegion::GenomicRegion(const Bed3RecordCore & bed) :
		GenomicRegion(
				bib::pasteAsStr(bed.chrom_, "-", bed.chromStart_, "-", bed.chromEnd_),
				bed.chrom_, bed.chromStart_, bed.chromEnd_, false) {

}


GenomicRegion::GenomicRegion(const TandemRepeatFinderRecord & trfRecord):GenomicRegion(
		bib::pasteAsStr(trfRecord.repeatPatSeq_, "_x", trfRecord.numberOfAlignedRepeats_),
		trfRecord.seqName_, trfRecord.start_ - 1, trfRecord.end_, false) {

}


GenomicRegion::GenomicRegion(const GFFCore & gff) :
		GenomicRegion(
				bib::pasteAsStr(gff.seqid_, "_", gff.start_ - 1, "_", gff.end_),
				gff.seqid_, gff.start_ - 1, gff.end_, '-' == gff.strand_ ) {
	if (gff.hasAttr("ID")) {
		uid_.append("_" + urldecode(gff.getAttr("ID")));
	}

	if(gff.hasAttr("Name")){
		alternateUids_.emplace_back(uid_);
		uid_ = gff.getAttr("Name");
	}
	alternateUids_.emplace_back(uid_);
}

GenomicRegion::GenomicRegion(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData):
				uid_(bAln.Name),
				alternateUids_(VecStr{ bAln.Name }),
				chrom_(-1 == bAln.RefID ? std::string("*") : refData.at(bAln.RefID).RefName),
				start_(bAln.Position),
				end_(bAln.GetEndPosition()),
				reverseSrand_(bAln.IsReverseStrand()) {
}


void GenomicRegion::checkRegion(const BamTools::BamReader & bReader,
		const BamTools::RefVector & refs, const VecStr & twoBitRefs){
	if (!bib::in(chrom_, twoBitRefs)) {
		std::stringstream ss;
		ss << "Error, ref name: " << chrom_
				<< " not in TwoBit file, options are: " << std::endl;
		ss << bib::conToStr(twoBitRefs, ",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	if (!bib::has(refs, chrom_,
			[](const BamTools::RefData & rData, const std::string & refName) {
					return rData.RefName == refName;})) {
		std::stringstream ss;
		ss << "Error, ref name: " << chrom_
				<< " not in Bam file, options are: " << std::endl;
		for(const auto & refD : iter::enumerate(refs)){
			if(0 != refD.index){
				ss << ",";
			}
			ss << refD.element.RefName;
		}
		ss << std::endl;
		throw std::runtime_error{ss.str()};
	}
	for (const auto & refN : refs) {
		if (refN.RefName == chrom_) {
			if(refN.RefLength < 0){
				std::stringstream ss;
				ss << "Error, ref length for " << chrom_ << " is less than 0: " << refN.RefLength << std::endl;
				throw std::runtime_error{ss.str()};
			}else	if (start_ > static_cast<size_t>(refN.RefLength)) {
				std::stringstream ss;
				ss << "Error, start: " << start_ << " is greater than " << chrom_
						<< " length: " << refN.RefLength << std::endl;
				throw std::runtime_error{ss.str()};
			}
			break;
		}
	}
	if(end_ <=start_){
		std::stringstream ss;
		ss << "Error, end should be greater than start" << std::endl;
		ss << " Start: " << start_ << "End:" << end_ << std::endl;
		throw std::runtime_error{ss.str()};
	}
	auto refId = bReader.GetReferenceID(chrom_);
	if (end_ > static_cast<size_t>(refs[refId].RefLength)) {
		std::cerr << "End was set to greater than " << chrom_ << ", length: "
				<< refs[refId].RefLength << " so it is being set to the end length"
				<< std::endl;
		end_ = refs[refId].RefLength;
	}
	if(std::numeric_limits<size_t>::max() == end_){
		end_ = refs[refId].RefLength;
	}
}


bool GenomicRegion::sameRegion(const GenomicRegion & otherRegion)const{
	return otherRegion.chrom_ == chrom_ &&
			otherRegion.start_ == start_ &&
			otherRegion.end_ == end_;
}
bool GenomicRegion::sameRegion(const GFFCore & gff)const{
	//have to take into account that posiitions in gff are 1 based
	return gff.seqid_ == chrom_ &&
			gff.start_ - 1 == start_ &&
			gff.end_ == end_;
}

bool GenomicRegion::overlaps(const GenomicRegion & otherRegion,
		const size_t overlapMin) const {
	if(sameRegion(otherRegion)){
		return true;
	}
	if(getOverlapLen(otherRegion) >=overlapMin){
		return true;
	}
	return false;
}

bool GenomicRegion::overlaps(const GFFCore & gff,
		const size_t overlapMin) const {
	if(sameRegion(gff)){
		return true;
	}
	if(getOverlapLen(gff) >=overlapMin){
		return true;
	}
	return false;
}



bool GenomicRegion::endsInThisRegion(
		const GenomicRegion & otherRegion) const {
	if(otherRegion.chrom_ != chrom_){
		return false;
	}
	return otherRegion.end_ > start_ && otherRegion.end_ <= end_;
}

bool GenomicRegion::startsInThisRegion(
		const GenomicRegion & otherRegion) const {
	if(otherRegion.chrom_ != chrom_){
		return false;
	}
	return otherRegion.start_ >= start_ && otherRegion.start_ < end_;
}


size_t GenomicRegion::getOverlapLen(const GenomicRegion & otherRegion) const {
	if(otherRegion.chrom_ != chrom_){
		return 0;
	}
	if( (otherRegion.start_ > start_ && otherRegion.start_ < end_) ||
			(otherRegion.end_ > start_ && otherRegion.end_ < end_) ||
			(start_ > otherRegion.start_ &&  start_ <  otherRegion.end_) ||
			(end_ > otherRegion.start_ && end_  < otherRegion.end_)  ) {
		auto overlapStart = std::max(otherRegion.start_, start_);
		auto overlapEnd = std::min(otherRegion.end_, end_);
		return overlapEnd - overlapStart;
	}
	return 0;
}

size_t GenomicRegion::getOverlapLen(const GFFCore & gff) const {
	auto chrom = gff.seqid_;
	auto start = gff.start_ -1 ;
	auto end = gff.end_;
	if(chrom != chrom_){
		return 0;
	}
	if( (start > start_ && start < end_) ||
			(end > start_ && end < end_) ||
			(start_ > start &&  start_ <  end) ||
			(end_ > start && end_  < end)  ) {
		auto overlapStart = std::max(start, start_);
		auto overlapEnd = std::min(end, end_);
		return overlapEnd - overlapStart;
	}
	return 0;
}



void GenomicRegion::addAltUid(const std::string & altUid) {
	alternateUids_.emplace_back(altUid);
}

Bed6RecordCore GenomicRegion::genBedRecordCore() const {
	return Bed6RecordCore(chrom_, start_, end_, uid_, getLen(),
			reverseSrand_ ? '-' : '+');
}

Bed3RecordCore GenomicRegion::genBed3RecordCore() const{
	return Bed3RecordCore(chrom_, start_, end_);
}

void GenomicRegion::setLongestUid(){
	std::string longestUid = "";
	for(const auto & uid : alternateUids_){
		if(uid.length() > longestUid.length()){
			longestUid = uid;
		}
	}
	uid_ = longestUid;
}


uint32_t GenomicRegion::getLen() const{
	return end_ - start_;
}


std::string GenomicRegion::createUidFromCoords() const {
//	MetaDataInName meta;
//	meta.addMeta("chrom", chrom_);
//	meta.addMeta("start_", start_);
//	meta.addMeta("end_", end_);
//	return meta.createMetaName();
	//return bib::pasteAsStr(chrom_, ":", start_, "-", end_);
	return bib::pasteAsStr(chrom_, "-", start_, "-", end_);
}
std::string GenomicRegion::createUidFromCoordsStrand() const {
//	MetaDataInName meta;
//	meta.addMeta("chrom", chrom_);
//	meta.addMeta("start_", start_);
//	meta.addMeta("end_", end_);
//	return meta.createMetaName();
	//return bib::pasteAsStr(chrom_, ":", start_, "-", end_);
	return bib::pasteAsStr(chrom_, "-", start_, "-", end_, "-", reverseSrand_ ? "rev" : "for");
}

double GenomicRegion::getPercInRegion(const BamTools::BamAlignment & bAln,
			const std::string & chrom) const{
	if (chrom_ != chrom) {
		return 0;
	}
	if(0 == bAln.AlignedBases.size()){
		return 0;
	}

	if(bAln.GetEndPosition() < start_){
		return 0;
	}
	if(bAln.Position >= end_){
		return 0;
	}

	if((bAln.Position >= start_ && bAln.Position < end_) ||
			(bAln.GetEndPosition() > start_ && bAln.GetEndPosition() <= end_)){

		auto start = std::max<size_t>(bAln.Position, start_);
		auto end = std::min<size_t>(bAln.GetEndPosition(), end_);
		double basesInRegion = end - start;
		return basesInRegion / bAln.AlignedBases.size();
	}else{
		return 0;
	}


}


size_t GenomicRegion::getRelativePositionFromStartStrandAware(
		size_t strandAwarePositionFromStart) const {
	if (strandAwarePositionFromStart >= getLen()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error strandAwareStart: "
				<< strandAwarePositionFromStart << " must be less than length: " << getLen()
				<< " for " << genBed3RecordCore().toDelimStr() << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(reverseSrand_){
		return end_ - 1 - strandAwarePositionFromStart;
	}else{
		return start_ + strandAwarePositionFromStart;
	}
}


seqInfo GenomicRegion::extractSeq(TwoBit::TwoBitFile & twobitReader) const{
	std::string buffer = "";
	twobitReader[chrom_]->getSequence(buffer, start_, end_, false);
	if(reverseSrand_){
		buffer = seqUtil::reverseComplement(buffer, "DNA");
	}
	std::string bName = bfs::basename(twobitReader.getFilename());
	return seqInfo(bName + "-" + uid_, buffer);
}

double GenomicRegion::getPercInRegion(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData) const {
	if(bAln.RefID < 0){
		return 0;
	}
	return getPercInRegion(bAln, refData.at(bAln.RefID).RefName);
}

void GenomicRegion::setUidWtihCoords(){
	uid_ = createUidFromCoords();
}

void GenomicRegion::setUidWtihCoordsStrand(){
	uid_ = createUidFromCoordsStrand();
}


Json::Value GenomicRegion::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["uid_"] = bib::json::toJson(uid_);
	ret["alternateUids_"] = bib::json::toJson(alternateUids_);
	ret["chrom_"] = bib::json::toJson(chrom_);
	ret["start_"] = bib::json::toJson(start_);
	ret["end_"] = bib::json::toJson(end_);
	ret["reverseSrand_"] = bib::json::toJson(reverseSrand_);
	ret["off_"] = bib::json::toJson(off_);
	return ret;
}

std::vector<GenomicRegion> gatherRegions(const std::string & bedFile,
		const std::string & gffFile, bool verbose){
	std::vector<GenomicRegion> allRegions;
	if ("" != bedFile) {
		auto beds = getBeds(bedFile);
		for(const auto & bed : beds){
			allRegions.emplace_back(*bed);
		}
	}
	if("" != gffFile){
		auto gffs = getGFFs(gffFile);
		for(const auto & gff : gffs){
			allRegions.emplace_back(*gff);
		}
	}
	if("" != bedFile || "" != gffFile){
		//remove duplicate regions
		for (const auto regionPos : iter::range(allRegions.size() - 1)) {
			if(allRegions[regionPos].off_){
				continue;
			}
			for(const auto secondRegionPos : iter::range(regionPos + 1,allRegions.size())){
				if(allRegions[secondRegionPos].off_){
					continue;
				}
				if(allRegions[regionPos].sameRegion(allRegions[secondRegionPos])){
					if(verbose){
						std::cout << "Found regions with same coordinates for "
								<< allRegions[regionPos].chrom_
								<< " " << allRegions[regionPos].start_
								<< " " << allRegions[regionPos].end_
								<< " , collapsing " << std::endl;
					}
					allRegions[regionPos].addAltUid(allRegions[secondRegionPos].uid_);
					allRegions[secondRegionPos].off_ = true;
				}
			}
		}
		allRegions.erase(
				std::remove_if(allRegions.begin(), allRegions.end(),
						[](const GenomicRegion & region) {return region.off_;}),
				allRegions.end());
		for(auto & region : allRegions){
			region.setLongestUid();
		}
	}
	return allRegions;
}



std::vector<GenomicRegion> bedPtrsToGenomicRegs(
		const std::vector<std::shared_ptr<Bed6RecordCore>> & beds) {
	return bib::convert<std::shared_ptr<Bed6RecordCore>, GenomicRegion>(beds,
			[](const std::shared_ptr<Bed6RecordCore> & bed)->GenomicRegion {return GenomicRegion(*bed);});
}

std::vector<GenomicRegion> bed3PtrsToGenomicRegs(
		const std::vector<std::shared_ptr<Bed3RecordCore>> & beds){
	return bib::convert<std::shared_ptr<Bed3RecordCore>, GenomicRegion>(beds,
			[](const std::shared_ptr<Bed3RecordCore> & bed)->GenomicRegion {return GenomicRegion(*bed);});
}

std::vector<GenomicRegion> gffPtrsToGenomicRegs(
		const std::vector<std::shared_ptr<GFFCore>> & gffs){
	return bib::convert<std::shared_ptr<GFFCore>, GenomicRegion>(gffs,
			[](const std::shared_ptr<GFFCore> & gff)->GenomicRegion {return GenomicRegion(*gff);});
}

void sortGRegionsByStart(std::vector<GenomicRegion> & regions){
	bib::sort(regions, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
		if(reg1.chrom_ == reg2.chrom_){
			return reg1.start_ < reg2.start_;
		}else{
			return reg1.chrom_ < reg2.chrom_;
		}
	});
}

}  // namespace bibseq
