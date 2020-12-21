#pragma once
/*
 * BedUtility.hpp
 *
 *  Created on: Jan 17, 2018
 *      Author: nick
 */


// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//


#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"

namespace njhseq {


class BedUtility {
public:

	static Bed3RecordCore & extendLeft(Bed3RecordCore & b, uint32_t extendLeft);
	static Bed3RecordCore & extendRight(Bed3RecordCore & b, uint32_t extendRight, uint32_t chromLength);
	static Bed3RecordCore & extendRight(Bed3RecordCore & b, uint32_t extendRight);
	static Bed3RecordCore & extendLeftRight(Bed3RecordCore & b, uint32_t extendLeft, uint32_t extendRight, uint32_t chromLength);
	static Bed3RecordCore & extendLeftRight(Bed3RecordCore & b, uint32_t extendLeft, uint32_t extendRight);

	static GenomicRegion & extendLeft(GenomicRegion & b, uint32_t extendLeft);
	static GenomicRegion & extendRight(GenomicRegion & b, uint32_t extendRight, uint32_t chromLength);
	static GenomicRegion & extendRight(GenomicRegion & b, uint32_t extendRight);
	static GenomicRegion & extendLeftRight(GenomicRegion & b, uint32_t extendLeft, uint32_t extendRight, uint32_t chromLength);
	static GenomicRegion & extendLeftRight(GenomicRegion & b, uint32_t extendLeft, uint32_t extendRight);


	template<typename T>
	static void coordSort(std::vector<T> & beds, bool decending = false){
		auto bedCoordSorterFunc =
				[](const T & reg1In, const T & reg2In) {
					const auto & reg1 = getRef(reg1In);
					const auto & reg2 = getRef(reg2In);
					if(reg1.chrom_ == reg2.chrom_) {
						if(reg1.chromStart_ == reg2.chromStart_) {
							return reg1.chromEnd_ < reg2.chromEnd_;
						} else {
							return reg1.chromStart_ < reg2.chromStart_;
						}
					} else {
						return reg1.chrom_ < reg2.chrom_;
					}
				};
		if(decending){
			std::sort(beds.rbegin(), beds.rend(), bedCoordSorterFunc);
		}else{
			njh::sort(beds, bedCoordSorterFunc);
		}
	}

	template<typename BED>
	static uint32_t getPlotIDForBed(const BED & record,
			std::unordered_map<std::string,
					std::unordered_map<uint32_t, std::vector<uint32_t>>> & alreadyTakenIds) {
		uint32_t id = 0;
		std::set<uint32_t> alreadyTaken;
		for (const auto pos : iter::range(getRef(record).chromStart_, getRef(record).chromEnd_)) {
			for (const auto & otherId : alreadyTakenIds[getRef(record).chrom_][pos]) {
				alreadyTaken.emplace(otherId);
			}
		}
		while (njh::in(id, alreadyTaken)) {
			++id;
		}
		for (const auto pos : iter::range(getRef(record).chromStart_, getRef(record).chromEnd_)) {
			alreadyTakenIds[getRef(record).chrom_][pos].emplace_back(id);
		}
		return id;
	}


	struct genSubRegionCombosPars{
		bool includeFromFront = false;
		bool includeFromBack = false;
		uint32_t includeFromSubRegionSize = std::numeric_limits<uint32_t>::max();
		uint32_t minLen = 1;
		uint32_t maxLen = std::numeric_limits<uint32_t>::max();
		uint32_t minBlockRegionLen = 1;
		bool justToNextRegion = false;
	};


	struct SubRegionCombo {
		SubRegionCombo(
				const Bed6RecordCore & startRegion,
				const Bed6RecordCore & endRegion);
		Bed6RecordCore startRegion_;
		Bed6RecordCore endRegion_;

		Bed6RecordCore genOut(uint32_t maximumToInclude = std::numeric_limits<uint32_t>::max()) const ;
		Bed6RecordCore genSubStart(uint32_t maximumToInclude = std::numeric_limits<uint32_t>::max()) const ;

		Bed6RecordCore genSubEnd(uint32_t maximumToInclude = std::numeric_limits<uint32_t>::max()) const;

	};

	template<typename BED6>
	static std::vector<SubRegionCombo> genSubRegionCombos(std::vector<BED6> regionsRaw,
			const std::unordered_map<std::string, uint32_t> & genomeLen,
			const genSubRegionCombosPars & pars){
		std::vector<SubRegionCombo> ret;
		BedUtility::coordSort(regionsRaw, false);
		std::vector<std::shared_ptr<Bed6RecordCore>> regions;
		regions.emplace_back(std::make_shared<Bed6RecordCore>(getRef(regionsRaw.front())));

		for (const auto regPos : iter::range<uint32_t>(1, regionsRaw.size())) {
			if (regions.back()->overlaps(getRef(regionsRaw[regPos]), 1) ||
					(regions.back()->chrom_ == getRef(regionsRaw[regPos]).chrom_ && regions.back()->chromEnd_ == getRef(regionsRaw[regPos]).chromStart_) ) {
				regions.back()->chromEnd_ = std::max(
						regions.back()->chromEnd_,
						getRef(regionsRaw[regPos]).chromEnd_);
			} else {
				regions.emplace_back(
						std::make_shared<Bed6RecordCore>(getRef(regionsRaw[regPos])));
			}
		}


		std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> regionsByChrom;
		for(const auto & region : regions){
			regionsByChrom[region->chrom_].emplace_back(region);
		}

		if(pars.includeFromFront){
			for(auto & byChrom : regionsByChrom){
				if(byChrom.second.front()->chromStart_ > 10){
					byChrom.second.emplace_back(std::make_shared<Bed6RecordCore>(byChrom.first, 0, 1, "start", 1, '+'));
					BedUtility::coordSort(byChrom.second, false);
				}
			}
		}
		if(pars.includeFromBack){
			for(auto & byChrom : regionsByChrom){
				if(!njh::in(byChrom.first, genomeLen)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "don't have chromosome " << byChrom.first << "\n";
					ss << "has: " << njh::conToStr(njh::getVecOfMapKeys(genomeLen), ",") << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(genomeLen.at(byChrom.first) - byChrom.second.back()->chromEnd_ > 10){
					byChrom.second.emplace_back(std::make_shared<Bed6RecordCore>(byChrom.first, genomeLen.at(byChrom.first) -1, genomeLen.at(byChrom.first), "back", 1, '+'));
				}
			}
		}


		for(auto & byChrom : regionsByChrom){
			if(byChrom.second.size() > 1){
				uint32_t lastDownStreamPos = 0;
				for(const auto pos : iter::range(byChrom.second.size() -1 )){
					Bed6RecordCore startReg = *byChrom.second[pos];
					if(startReg.length() < pars.minBlockRegionLen){
						continue;
					}
					if(pars.justToNextRegion && pos < lastDownStreamPos){
						continue;
					}
					for(const auto downStreamPos : iter::range(pos + 1, byChrom.second.size())){
						Bed6RecordCore endReg   = *byChrom.second[downStreamPos];
						if(endReg.length() < pars.minBlockRegionLen){
							continue;
						}

						SubRegionCombo combinedRegion(startReg, endReg);
						Bed6RecordCore outRegion = combinedRegion.genOut(pars.includeFromSubRegionSize);
						bool passMinLen = outRegion.length() >= pars.minLen;
						if(!passMinLen){
							if(startReg.length() + endReg.length() > pars.maxLen){
								passMinLen = true;
							}else if(0 == pos && endReg.length() > startReg.length() + endReg.length() > pars.maxLen){
								passMinLen = true;
							}else if(startReg.length() + endReg.length() > pars.maxLen && downStreamPos == (byChrom.second.size() -1)){
								passMinLen = true;
							}else if(0 == pos && downStreamPos == (byChrom.second.size() -1)){
								passMinLen = true;
							}
						}
						if(passMinLen && outRegion.length() <= pars.maxLen){
							ret.emplace_back(SubRegionCombo(startReg, endReg));
						}
						if(pars.justToNextRegion && passMinLen && outRegion.length() <= pars.maxLen){
							lastDownStreamPos = downStreamPos;
							break;
						}
					}
				}
			}
		}
		return ret;
	}

};

} /* namespace njhseq */

