#pragma once
/*
 * GenomicRegion.hpp
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
#include <api/BamReader.h>
#include <TwoBit.h>

#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"
#include "njhseq/objects/BioDataObject/Bed3RecordCore.hpp"
#include "njhseq/objects/BioDataObject/GFFCore.hpp"
#include "njhseq/objects/BioDataObject/TandemRepeatFinderRecord.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

namespace njhseq {

class GenomicRegion {
public:

	GenomicRegion();
	GenomicRegion(const std::string & uid, const std::string & chrom,
			size_t start, size_t end, bool reverseStrand);
	GenomicRegion(const Bed6RecordCore & bed);
	GenomicRegion(const Bed3RecordCore & bed);
	GenomicRegion(const GFFCore & gff);
	GenomicRegion(const TandemRepeatFinderRecord & trfRecord);
	GenomicRegion(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData);
	std::string uid_;
	VecStr alternateUids_;
	std::string chrom_;
	size_t start_;
	size_t end_;
	bool reverseSrand_ = false;
	bool off_ = false;

	MetaDataInName meta_;

	uint32_t getLen() const;
	Bed6RecordCore genBedRecordCore() const;
	Bed3RecordCore genBed3RecordCore() const;
	void addAltUid(const std::string & altUid);
	void setLongestUid();
	bool sameRegion(const GenomicRegion & otherRegion) const;
	bool sameRegion(const GFFCore & gff) const;
	bool sameRegion(const Bed3RecordCore & otherRegion) const;

	void checkRegion(const BamTools::BamReader & bReader,
			const BamTools::RefVector & refs, const VecStr & twoBitRefs);

	bool overlaps(const GenomicRegion & otherRegion,
			const size_t overlapMin = 1) const;
	bool overlaps(const Bed3RecordCore & otherRegion,
			const size_t overlapMin = 1) const;
	bool overlaps(const GFFCore & gff,
			const size_t overlapMin = 1) const;
	bool overlaps(const std::string & chrom, const size_t start, const size_t end,
			const size_t overlapMin = 1) const;


	size_t getOverlapLen(const GenomicRegion & otherRegion) const;
	size_t getOverlapLen(const GFFCore & gff) const;
	size_t getOverlapLen(const Bed3RecordCore & otherRegion) const;
	size_t getOverlapLen(const std::string & chrom, const size_t start, const size_t end) const;


	size_t distBetweenRegions(const GenomicRegion & otherRegion);
	size_t distBetweenRegions(const GFFCore & otherRegion);
	size_t distBetweenRegions(const Bed3RecordCore & otherRegion);
	/**@brief Can the distance between two different regions, if on different chroms returns std::numeric<uint32_t>::max(), overlaps return 0
	 *
	 * @param chrom the chromosome of the other region
	 * @param start the start on the chromsome
	 * @param end the end of the chromosome
	 * @return the distance between the regions
	 */
	size_t distBetweenRegions(const std::string & chrom, const size_t start, const size_t end);


	bool startsInThisRegion(const GenomicRegion & otherRegion) const;
	bool endsInThisRegion(const GenomicRegion & otherRegion) const;
	bool fallsInThisRegion(const GenomicRegion & otherRegion) const;

	bool startsInThisRegion(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData) const;
	bool endsInThisRegion(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData) const;
	bool fallsInThisRegion(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData) const;

	bool startsInThisRegion(const std::string & chrom, uint32_t start) const;
	bool endsInThisRegion(const std::string & chrom, uint32_t end) const;
	bool fallsInThisRegion(const std::string & chrom, uint32_t start,
			uint32_t end) const; //no check for if start is less than end


	Json::Value toJson() const;
  Json::Value toJsonLocationOnly() const;

	std::string createUidFromCoords() const;
	std::string createUidFromCoordsStrand() const;

	void setUidWtihCoords();
	void setUidWtihCoordsStrand();

	double getPercInRegion(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData) const;

	double getPercInRegion(const BamTools::BamAlignment & bAln,
				const std::string & chrom) const;

	seqInfo extractSeq(TwoBit::TwoBitFile & twobitReader) const;

	size_t getRelativePositionFromStartStrandAware(size_t strandAwarePositionFromStart) const;

};

std::vector<GenomicRegion> gatherRegions(const std::string & bedFile,
		const std::string & gffFile, bool verbose);

std::vector<GenomicRegion> bedPtrsToGenomicRegs(
		const std::vector<std::shared_ptr<Bed6RecordCore>> & beds);
std::vector<GenomicRegion> bed3PtrsToGenomicRegs(
		const std::vector<std::shared_ptr<Bed3RecordCore>> & beds);

std::vector<GenomicRegion> gffPtrsToGenomicRegs(
		const std::vector<std::shared_ptr<GFFCore>> & gffs);

void sortGRegionsByStart(std::vector<GenomicRegion> & regions);


template<typename IN, typename OUT>
std::vector<OUT> convertGenomeRegions(const std::vector<IN> & regions, std::function<OUT(const IN &)> unaryConvertor){
	std::vector<OUT> ret;
	for(const auto & reg : regions){
		ret.emplace_back(std::forward<OUT>(unaryConvertor(reg)));
	}
	return ret;
}



}  // namespace njhseq

