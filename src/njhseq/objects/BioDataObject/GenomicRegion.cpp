/*
 * GenomicRegion.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: nick
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include "GenomicRegion.hpp"
#include "njhseq/objects/BioDataObject/reading.hpp"
#include "njhseq/helpers/seqUtil.hpp"


namespace njhseq {

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
	if(bed.extraFields_.size() > 0){
		uint32_t extraFieldCount = 0;
		for(const auto & extraField : bed.extraFields_){
			if(MetaDataInName::nameHasMetaData(extraField)){
				meta_.addMeta(MetaDataInName(extraField), false);
			} else {
				meta_.addMeta(njh::pasteAsStr("extraField", leftPadNumStr<uint32_t>(extraFieldCount, bed.extraFields_.size())), extraField, false);
			}
			++extraFieldCount;
		}
	}
}

GenomicRegion::GenomicRegion(const Bed3RecordCore & bed) :
		GenomicRegion(
				njh::pasteAsStr(bed.chrom_, "-", bed.chromStart_, "-", bed.chromEnd_),
				bed.chrom_, bed.chromStart_, bed.chromEnd_, false) {
	//add name if extract fields present, might not always be name though
	if(bed.extraFields_.size() > 0){
		uid_ = bed.extraFields_[0];
	}
}


GenomicRegion::GenomicRegion(const TandemRepeatFinderRecord & trfRecord):GenomicRegion(
		njh::pasteAsStr(trfRecord.seqName_, "-",trfRecord.start_ - 1, "-",trfRecord.end_,
				"__",trfRecord.repeatPatSeq_, "_x", (trfRecord.end_ + 1 - trfRecord.start_)/static_cast<double>(trfRecord.periodSize_)),
		trfRecord.seqName_,
		trfRecord.start_ - 1,
		trfRecord.end_,
		false) {

}


GenomicRegion::GenomicRegion(const GFFCore & gff) :
		GenomicRegion(
				njh::pasteAsStr(gff.seqid_, "_", gff.start_ - 1, "_", gff.end_),
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
	if (!njh::in(chrom_, twoBitRefs)) {
		std::stringstream ss;
		ss << "Error, ref name: " << chrom_
				<< " not in TwoBit file, options are: " << std::endl;
		ss << njh::conToStr(twoBitRefs, ",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	if (!njh::has(refs, chrom_,
			[](const BamTools::RefData & rData, const std::string & refName) {
					return rData.RefName == refName;})) {
		std::stringstream ss;
		ss << "Error, ref name: " << chrom_
				<< " not in Bam file, options are: " << std::endl;
		for(const auto refD : iter::enumerate(refs)){
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


bool GenomicRegion::sameRegion(const GenomicRegion & otherRegion) const {
	return otherRegion.chrom_ == chrom_ && otherRegion.start_ == start_
			&& otherRegion.end_ == end_;
}

bool GenomicRegion::sameRegion(const Bed3RecordCore & otherRegion) const {
	return otherRegion.chrom_ == chrom_ && otherRegion.chromStart_ == start_
			&& otherRegion.chromEnd_ == end_;
}

bool GenomicRegion::sameRegion(const GFFCore & gff) const {
	//have to take into account that posiitions in gff are 1 based
	return gff.seqid_ == chrom_ && gff.start_ - 1 == start_ && gff.end_ == end_;
}

bool GenomicRegion::overlaps(const GenomicRegion & otherRegion,
		const size_t overlapMin) const {
	return getOverlapLen(otherRegion) >=overlapMin;
}

bool GenomicRegion::overlaps(const Bed3RecordCore & otherRegion,
		const size_t overlapMin) const{
	return getOverlapLen(otherRegion) >=overlapMin;
}

bool GenomicRegion::overlaps(const GFFCore & gff,
		const size_t overlapMin) const {
	return getOverlapLen(gff) >=overlapMin;
}

bool GenomicRegion::overlaps(const std::string & chrom, const size_t start, const size_t end,
		const size_t overlapMin) const{
	return getOverlapLen(chrom, start, end) >=overlapMin;
}

size_t GenomicRegion::distBetweenRegions(const GenomicRegion & otherRegion){
	return distBetweenRegions(otherRegion.chrom_, otherRegion.start_, otherRegion.end_);
}
size_t GenomicRegion::distBetweenRegions(const GFFCore & otherRegion){
	return distBetweenRegions(otherRegion.seqid_, otherRegion.start_, otherRegion.end_);

}
size_t GenomicRegion::distBetweenRegions(const Bed3RecordCore & otherRegion){
	return distBetweenRegions(otherRegion.chrom_, otherRegion.chromStart_, otherRegion.chromEnd_);

}

size_t GenomicRegion::distBetweenRegions(const std::string & otherChrom, const size_t otherStart, const size_t otherend){
	if(chrom_ != otherChrom){
		return std::numeric_limits<size_t>::max();
	}
	if(overlaps(otherChrom, otherStart, otherend)){
		return 0;
	}
	if(start_ < otherStart){
		return otherStart - end_;
	}else{
		return start_ - otherend;
	}
}


bool GenomicRegion::startsInThisRegion(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData) const {
	return startsInThisRegion(refData[bAln.RefID].RefName, bAln.Position);
}

bool GenomicRegion::endsInThisRegion(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData) const {
	return endsInThisRegion(refData[bAln.RefID].RefName, bAln.GetEndPosition());

}

bool GenomicRegion::fallsInThisRegion(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData) const {

	return fallsInThisRegion(refData[bAln.RefID].RefName, bAln.Position, bAln.GetEndPosition());
}

bool GenomicRegion::startsInThisRegion(const std::string & chrom,
		uint32_t start) const { //no check for if start is less than end
	if (chrom != chrom_) {
		return false;
	}
	return start >= start_ && start < end_;
}

bool GenomicRegion::endsInThisRegion(const std::string & chrom,
		uint32_t end) const { //no check for if start is less than end
	if (chrom != chrom_) {
		return false;
	}
	return end > start_ && end <= end_;
}

bool GenomicRegion::fallsInThisRegion(const std::string & chrom, uint32_t start,
		uint32_t end) const { //no check for if start is less than end
	if (chrom != chrom_) {
		return false;
	}
	return overlaps(chrom, start, end);
}




bool GenomicRegion::endsInThisRegion(
		const GenomicRegion & otherRegion) const {
	return endsInThisRegion(otherRegion.chrom_, otherRegion.end_);
}

bool GenomicRegion::startsInThisRegion(
		const GenomicRegion & otherRegion) const {
	return startsInThisRegion(otherRegion.chrom_, otherRegion.start_);
}

bool GenomicRegion::fallsInThisRegion(
		const GenomicRegion & otherRegion) const {
	return fallsInThisRegion(otherRegion.chrom_, otherRegion.start_, otherRegion.end_);
}


size_t GenomicRegion::getOverlapLen(const GenomicRegion & otherRegion) const {
	//std::cout << "otherRegion.uid_: " << otherRegion.uid_ << std::endl;
	return getOverlapLen(otherRegion.chrom_, otherRegion.start_, otherRegion.end_);
}

size_t GenomicRegion::getOverlapLen(const GFFCore & gff) const {
	return getOverlapLen(gff.seqid_, gff.start_ -1, gff.end_);
}

size_t GenomicRegion::getOverlapLen(const Bed3RecordCore & otherRegion) const {
	return getOverlapLen(otherRegion.chrom_, otherRegion.chromStart_, otherRegion.chromEnd_);
}

size_t GenomicRegion::getOverlapLen(const std::string & otherChrom,
		const size_t otherStart, const size_t otherEnd) const {
	if (otherChrom != chrom_) {
		return 0;
	}
//	if(1724816	== otherStart && 1726997 == otherEnd ){
//		std::cout << "uid_: " << uid_ << std::endl;
//		std::cout << "otherStart: " << otherStart << std::endl;
//		std::cout << "start_:     " << start_ << std::endl;
//		std::cout << "otherEnd:   " << otherEnd << std::endl;
//		std::cout << "end_:       " << end_ << std::endl;
//		std::cout << "otherStart >= start_  " << njh::colorBool(otherStart >= start_) << std::endl;
//		std::cout << "otherStart <  end_:   " << njh::colorBool(otherStart <  end_) << std::endl;
//		std::cout << "otherEnd    > start_: " << njh::colorBool(otherEnd    > start_) << std::endl;
//		std::cout << "otherEnd   <= end_:   " << njh::colorBool(otherEnd   <= end_) << std::endl;
//		std::cout << "otherStart >= start_ && otherStart <  end_: " << njh::colorBool(otherStart >= start_ && otherStart <  end_) << std::endl;
//		std::cout << "otherEnd    > start_ && otherEnd   <= end_: " << njh::colorBool(otherEnd    > start_ && otherEnd   <= end_) << std::endl;
//		std::cout << "start_ >= otherStart && start_ <  otherEnd: " << njh::colorBool(start_ >= otherStart && start_ <  otherEnd) << std::endl;
//		std::cout << "end_    > otherStart && end_   <= otherEnd: " << njh::colorBool(end_    > otherStart && end_   <= otherEnd) << std::endl;
//
//	}

	if ((otherStart >= start_ && otherStart <  end_) ||
			(otherEnd    > start_ && otherEnd   <= end_) ||
			(start_ >= otherStart && start_ <  otherEnd) ||
			(end_    > otherStart && end_   <= otherEnd)) {
		auto overlapStart = std::max(otherStart, start_);
		auto overlapEnd =   std::min(otherEnd,   end_);
		return overlapEnd - overlapStart;
	}
	return 0;
}



void GenomicRegion::addAltUid(const std::string & altUid) {
	alternateUids_.emplace_back(altUid);
}

Bed6RecordCore GenomicRegion::genBedRecordCore() const {
	Bed6RecordCore ret (chrom_, start_, end_, uid_, getLen(),
			reverseSrand_ ? '-' : '+');
	if("" == uid_){
		ret.name_ = createUidFromCoordsStrand();
	}
	if(!meta_.meta_.empty()){
		ret.extraFields_.emplace_back(meta_.createMetaName());
	}
	return ret;
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
	//return njh::pasteAsStr(chrom_, ":", start_, "-", end_);
	return njh::pasteAsStr(chrom_, "-", start_, "-", end_);
}
std::string GenomicRegion::createUidFromCoordsStrand() const {
//	MetaDataInName meta;
//	meta.addMeta("chrom", chrom_);
//	meta.addMeta("start_", start_);
//	meta.addMeta("end_", end_);
//	return meta.createMetaName();
	//return njh::pasteAsStr(chrom_, ":", start_, "-", end_);
	return njh::pasteAsStr(chrom_, "-", start_, "-", end_, "-", reverseSrand_ ? "rev" : "for");
}

double GenomicRegion::getPercInRegion(const BamTools::BamAlignment & bAln,
			const std::string & chrom) const{
	if (chrom_ != chrom) {
		return 0;
	}
	if(0 == bAln.AlignedBases.size()){
		return 0;
	}

	if(bAln.GetEndPosition() < static_cast<int64_t>(start_)){
		return 0;
	}
	if(bAln.Position >= static_cast<int64_t>(end_)){
		return 0;
	}

	//get percent overlap if
	//1) bAln.Position falls within the region
	//2) baln.GetEndPosition() falls within the region
	//3) if the bAln covers the whole region and beyond
//	std::cout << njh::bashCT::red;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if((     bAln.Position >= static_cast<int64_t>(start_)        && bAln.Position < static_cast<int64_t>(end_)) ||
			(bAln.GetEndPosition() > static_cast<int64_t>(start_) && bAln.GetEndPosition() <= static_cast<int64_t>(end_)) ||
			(bAln.Position <= static_cast<int64_t>(start_)        && bAln.GetEndPosition() >= static_cast<int64_t>(end_))){
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "bAln.Position: " << bAln.Position << std::endl;
//		std::cout << "bAln.GetEndPosition(): " << bAln.GetEndPosition() << std::endl;
//		std::cout << "start_: " << start_ << std::endl;
//		std::cout << "end_: " << end_ << std::endl;

		auto start = std::max<size_t>(bAln.Position, start_);
		auto end = std::min<size_t>(bAln.GetEndPosition(), end_);
//		std::cout << "start: " << start << std::endl;
//		std::cout << "end: " << end << std::endl;

		double basesInRegion = end - start;
//		std::cout << "basesInRegion: " << basesInRegion << std::endl;
//		std::cout << "std::min<uint64_t>(bAln.AlignedBases.size(), getLen()): " << std::min<uint64_t>(bAln.AlignedBases.size(), getLen()) << std::endl;
//		std::cout << njh::bashCT::reset;

		return basesInRegion / std::min<uint64_t>(bAln.AlignedBases.size(), getLen());
	}else{
//		std::cout << njh::bashCT::reset;

		return 0;
	}


}


size_t GenomicRegion::getRelativePositionFromStartStrandAware(
		size_t strandAwarePositionFromStart) const {
	if (strandAwarePositionFromStart > getLen()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error strandAwareStart: "
				<< strandAwarePositionFromStart << " must be less than length: " << getLen()
				<< " for " << genBed3RecordCore().toDelimStr() << "\n";
		ss << "If trying to get the end position, do end minus 1 and add one to the end" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(reverseSrand_){
		if(end_ == strandAwarePositionFromStart){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't get relative positive from reverse strand for " << strandAwarePositionFromStart << " because "<< "\n";
			throw std::runtime_error{ss.str()};
		}
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
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["uid_"] = njh::json::toJson(uid_);
	ret["alternateUids_"] = njh::json::toJson(alternateUids_);
	ret["chrom_"] = njh::json::toJson(chrom_);
	ret["start_"] = njh::json::toJson(start_);
	ret["end_"] = njh::json::toJson(end_);
	ret["reverseSrand_"] = njh::json::toJson(reverseSrand_);
	ret["off_"] = njh::json::toJson(off_);
	return ret;
}

Json::Value GenomicRegion::toJsonLocationOnly() const{
  Json::Value ret;
  ret["chrom"] = njh::json::toJson(chrom_);
  ret["start"] = njh::json::toJson(start_);
  ret["end"] = njh::json::toJson(end_);
  ret["strand"] = njh::json::toJson(reverseSrand_ ? '-' : '+');
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
	return njh::convert<std::shared_ptr<Bed6RecordCore>, GenomicRegion>(beds,
			[](const std::shared_ptr<Bed6RecordCore> & bed)->GenomicRegion {
		return GenomicRegion(*bed);
	});
}

std::vector<GenomicRegion> bed3PtrsToGenomicRegs(
		const std::vector<std::shared_ptr<Bed3RecordCore>> & beds){
	return njh::convert<std::shared_ptr<Bed3RecordCore>, GenomicRegion>(beds,
			[](const std::shared_ptr<Bed3RecordCore> & bed)->GenomicRegion {return GenomicRegion(*bed);});
}

std::vector<GenomicRegion> gffPtrsToGenomicRegs(
		const std::vector<std::shared_ptr<GFFCore>> & gffs){
	return njh::convert<std::shared_ptr<GFFCore>, GenomicRegion>(gffs,
			[](const std::shared_ptr<GFFCore> & gff)->GenomicRegion {return GenomicRegion(*gff);});
}

void sortGRegionsByStart(std::vector<GenomicRegion> & regions){
	njh::sort(regions, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
		if(reg1.chrom_ == reg2.chrom_){
			if(reg1.start_ == reg2.start_) {
				return reg1.end_ < reg2.end_;
			} else {
				return reg1.start_ < reg2.start_;
			}
		}else{
			return reg1.chrom_ < reg2.chrom_;
		}
	});
}

}  // namespace njhseq
