/*
 * GenomicRegion.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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


namespace bibseq {

GenomicRegion::GenomicRegion() :
		uid_(""), chrom_(""), start_(0), end_(0) {
}

GenomicRegion::GenomicRegion(const std::string & uid, const std::string & chrom,
		size_t start, size_t end, bool reverseSrand) :
		uid_(uid),alternateUids_{uid}, chrom_(chrom), start_(start), end_(end), reverseSrand_(reverseSrand) {
}

GenomicRegion::GenomicRegion(const BedRecordCore & bed) :
		GenomicRegion(bed.name_, bed.chrom_, bed.chromStart_, bed.chromEnd_, '-' == bed.strand_) {

}

GenomicRegion::GenomicRegion(const GFFCore & gff) :
		GenomicRegion(
				bib::pasteAsStr(gff.seqName_, "_", gff.start_ - 1, "_", gff.end_),
				gff.seqName_, gff.start_ - 1, gff.end_, '-' == gff.strand_) {
	if (gff.hasAttr("ID")) {
		uid_.append("_" + urldecode(gff.getAttr("ID")));
	}
	alternateUids_.emplace_back(uid_);
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


void GenomicRegion::addAltUid(const std::string & altUid){
	alternateUids_.emplace_back(altUid);
}
BedRecordCore GenomicRegion::genBedRecordCore()const{
	return BedRecordCore(chrom_, start_, end_, uid_, 255, reverseSrand_ ? '-' : '+');
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


}  // namespace bibseq
