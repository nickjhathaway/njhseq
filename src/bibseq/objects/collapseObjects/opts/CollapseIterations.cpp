/*
 * CollapseIterations.cpp
 *
 *  Created on: May 6, 2016
 *      Author: nick
 */
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
#include "CollapseIterations.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"

namespace bibseq {
CollapseIterations::CollapseIterations(const CollapseIterations& other) :
		iters_(other.iters_), onPerId_(other.onPerId_) {
	//std::cout << "copy constructor of CollapseIterations\n";

}

CollapseIterations& CollapseIterations::operator=(const CollapseIterations & other) {
	//std::cout << "copy assignment of CollapseIterations\n";
	onPerId_ = other.onPerId_;
	iters_ = other.iters_;
	return *this;
}

CollapseIterations& CollapseIterations::operator=(CollapseIterations&& other){
	//std::cout << "move assignment of CollapseIterations\n";
	onPerId_ = other.onPerId_;
	//std::cout << __PRETTY_FUNCTION__ << " : 1" << std::endl;
	itersNums_ =  std::move(other.itersNums_);
	//std::cout << __PRETTY_FUNCTION__ << " : 2" << std::endl;
	//iters_.clear();
	//std::cout << __PRETTY_FUNCTION__ << " : 2.5" << std::endl;
	//std::swap(iters_, other.iters_);
	iters_ = other.iters_;
	//std::cout << __PRETTY_FUNCTION__ << " : 3" << std::endl;
	//std::map<uint32_t, IterPar> temp(std::move(other.iters_));
	//iters_ =  std::move(other.iters_);
	//std::cout << __PRETTY_FUNCTION__ << " : 4" << std::endl;
	return *this;
}




CollapseIterations::CollapseIterations(){

}

CollapseIterations::CollapseIterations(bool onPerId): onPerId_(onPerId){

}

CollapseIterations::CollapseIterations(const std::string & parametersFilename,
		bool onPerId):onPerId_(onPerId) {
	table inTab(parametersFilename, ":");
	uint32_t iters = 1;
	for (const auto& row : inTab.content_) {
		if (row.empty() || row.front().front() == 's' || row.front().front() == '#'
				|| row.front().front() == 'S') {
			continue;
		}
		std::vector<double> tempVect;
		for (const auto & colPos : iter::range(row.size())) {
			if (colPos == 0) {
				if (stringToLowerReturn(row[colPos]) == "all") {
					tempVect.emplace_back(std::numeric_limits<uint32_t>::max());
				} else {
					tempVect.emplace_back(std::stod(row[colPos]));
				}
			} else {
				if (row[colPos] == "all") {
					tempVect.emplace_back(std::numeric_limits<uint32_t>::max());
				} else {
					tempVect.emplace_back(std::stod(row[colPos]));
				}
			}
		}
		addIteration(iters, tempVect);
		++iters;
	}
}

CollapseIterations::CollapseIterations(std::map<uint32_t, IterPar> pars,
		bool onPerId): iters_(pars), onPerId_(onPerId) {

}

void CollapseIterations::writePars(std::ostream & out) const {
	table outTab(IterPar::outputIterInfoHeader(onPerId_));
	for (const auto & iter : iters_) {
		outTab.content_.emplace_back(iter.second.outputIterInfo());
	}
	outTab.outPutContents(out, ":");
}


void CollapseIterations::addIteration(uint32_t iterNum, std::vector<double> pars){
	IterPar ipar(pars, iterNum, IterPar::PerIdPars{onPerId_, onHqPerId_});
	addIteration(iterNum, ipar);
}

void CollapseIterations::addIteration(uint32_t iterNum, const IterPar & par) {
	if (bib::in(iterNum, iters_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " : " << "Iteration number: "
				<< bib::bashCT::boldRed(estd::to_string(iterNum))
				<< " already in iterations ";
		throw std::runtime_error { ss.str() };
	}
	itersNums_[iterNum] = iterNum;
	iters_[iterNum] = par;
}

void CollapseIterations::addIteration(const IterPar & par) {
	uint32_t iterNum = iters_.size() + 1;
	addIteration(iterNum, par);
}

void CollapseIterations::addIterations(std::vector<std::vector<double>> pars) {
	for(const auto pos : iter::range(pars.size())){
		addIteration(pos + 1, pars[pos]);
	}
}

CollapseIterations CollapseIterations::gen454ItDefaultPars(uint32_t stopCheck) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl,3,1,0,0,0,0,1},
			{stopCheckDbl,3,2,0,0,0,0,1},
			{stopCheckDbl,3,3,0,0,0,1,1},
			{stopCheckDbl,3,4,0,.99,0,2,1},
			{stopCheckDbl,3,5,0,.99,0,3,1},
			{stopCheckDbl,3,6,0,.99,0,4,1},
			{stopCheckDbl,3,7,1,.99,0,5,1},
			{stopCheckDbl,3,7,2,.99,0,5,1},
			{stopCheckDbl,3,7,3,.99,0,5,1},
			{stopCheckDbl,3,7,4,.99,0,5,2},
			{stopCheckDbl,3,7,5,.99,0,5,2},
			{stopCheckDbl,0,1,0,0,0,0,1},
			{stopCheckDbl,0,2,0,0,0,0,1},
			{stopCheckDbl,0,3,0,0,0,1,1},
			{stopCheckDbl,0,4,0,.99,0,2,1},
			{stopCheckDbl,0,5,0,.99,0,3,1},
			{stopCheckDbl,0,6,0,.99,0,4,1},
			{stopCheckDbl,0,7,1,.99,0,5,1},
			{stopCheckDbl,0,7,2,.99,0,5,1},
			{stopCheckDbl,0,7,3,.99,0,5,1},
			{stopCheckDbl,0,7,4,.99,0,5,2},
			{stopCheckDbl,0,7,5,.99,0,5,2},
			{stopCheckDbl,0,0,0,0,0,0,0}
		};
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genIlluminaDefaultPars(uint32_t stopCheck) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl,3,1,0,0,0,0,1},
			{stopCheckDbl,3,2,0,0,0,0,1},
			{stopCheckDbl,3,2,0,0,0,1,1},
			{stopCheckDbl,3,2,0,0,0,2,1},
			{stopCheckDbl,3,2,0,0,0,8,1},
			{stopCheckDbl,3,2,0,0,0,10,2},
			{stopCheckDbl,3,2,0,0,0,12,2},
			{stopCheckDbl,3,2,0,0,0,14,2},
			{stopCheckDbl,3,2,0,0,0,14,2},
			{stopCheckDbl,3,2,0,0,0,14,2},
			{stopCheckDbl,0,1,0,0,0,0,1},
			{stopCheckDbl,0,2,0,0,0,0,1},
			{stopCheckDbl,0,2,0,0,0,1,1},
			{stopCheckDbl,0,2,0,0,0,2,1},
			{stopCheckDbl,0,2,0,0,0,8,1},
			{stopCheckDbl,0,2,0,0,0,10,2},
			{stopCheckDbl,0,2,0,0,0,12,2},
			{stopCheckDbl,0,2,0,0,0,14,2},
			{stopCheckDbl,0,2,0,0,0,14,2},
			{stopCheckDbl,0,2,0,0,0,14,2},
			{stopCheckDbl,0,0,0,0,0,0,0}
	};
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genStrictNoErrorsDefaultPars(uint32_t stopCheck) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl, 3,0,0,0,0,0,0},
			{stopCheckDbl, 3,0,0,0,0,0,0},
			{stopCheckDbl, 0,0,0,0,0,0,0},
			{stopCheckDbl, 0,0,0,0,0,0,0}
	};
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genStrictDefaultPars(uint32_t stopCheck) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl, 3,0,0,0,0,0,0},
			{stopCheckDbl, 3,1,1,0,0,1,0},
			{stopCheckDbl, 0,0,0,0,0,0,0},
			{stopCheckDbl, 0,1,1,0,0,1,0},
			{stopCheckDbl, 0,0,0,0,0,0,0}
	};
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}


CollapseIterations CollapseIterations::gen454ItDefaultParsWithHqs(
		uint32_t stopCheck, uint32_t hqMismatches) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl,3,1,0,0,0,0,1},
			{stopCheckDbl,3,2,0,0,0,0,1},
			{stopCheckDbl,3,3,0,0,0,1,1},
			{stopCheckDbl,3,4,0,.99,0,2,1},
			{stopCheckDbl,3,5,0,.99,0,3,1},
			{stopCheckDbl,3,6,0,.99,0,4,1},
			{stopCheckDbl,3,7,1,.99,0,5,1},
			{stopCheckDbl,3,7,2,.99,0,5,1},
			{stopCheckDbl,3,7,3,.99,0,5,1},
			{stopCheckDbl,3,7,4,.99,0,5,2},
			{stopCheckDbl,3,7,5,.99,0,5,2}};
	for(double hq : iter::range<uint32_t>(1, hqMismatches + 1) ){
		iters.emplace_back(std::vector<double>{stopCheckDbl,3,7,5,.99,hq,5,2});
		iters.emplace_back(std::vector<double>{stopCheckDbl,3,7,5,.99,hq,5,2});
	}
	addOtherVec(iters,std::vector<std::vector<double>>{{stopCheckDbl,0,1,0,0,0,0,1},
			{stopCheckDbl,0,2,0,0,0,0,1},
			{stopCheckDbl,0,3,0,0,0,1,1},
			{stopCheckDbl,0,4,0,.99,0,2,1},
			{stopCheckDbl,0,5,0,.99,0,3,1},
			{stopCheckDbl,0,6,0,.99,0,4,1},
			{stopCheckDbl,0,7,1,.99,0,5,1},
			{stopCheckDbl,0,7,2,.99,0,5,1},
			{stopCheckDbl,0,7,3,.99,0,5,1},
			{stopCheckDbl,0,7,4,.99,0,5,2},
			{stopCheckDbl,0,7,5,.99,0,5,2}
		});
	for(double hq : iter::range<uint32_t>(1, hqMismatches + 1) ){
		iters.emplace_back(std::vector<double>{stopCheckDbl,0,7,5,.99,hq,5,2});
		iters.emplace_back(std::vector<double>{stopCheckDbl,0,7,5,.99,hq,5,2});
	}
	iters.emplace_back(std::vector<double>{stopCheckDbl,0,0,0,0,0,0,0});
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genIlluminaDefaultParsWithHqs(
		uint32_t stopCheck, uint32_t hqMismatches) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl,3,1,0,0,0,0,1},
			{stopCheckDbl,3,2,0,0,0,0,1},
			{stopCheckDbl,3,2,0,0,0,1,1},
			{stopCheckDbl,3,2,0,0,0,2,1},
			{stopCheckDbl,3,2,0,0,0,8,1},
			{stopCheckDbl,3,2,0,0,0,10,2},
			{stopCheckDbl,3,2,0,0,0,12,2},
			{stopCheckDbl,3,2,0,0,0,14,2},
			{stopCheckDbl,3,2,0,0,0,14,2},
			{stopCheckDbl,3,2,0,0,0,14,2} };
	for(double hq : iter::range<uint32_t>(1, hqMismatches + 1) ){
		iters.emplace_back(std::vector<double>{stopCheckDbl,3,2,0,0,hq,14,2});
		iters.emplace_back(std::vector<double>{stopCheckDbl,3,2,0,0,hq,14,2});
	}
	iters.emplace_back(std::vector<double>{stopCheckDbl,0,0,0,0,0,0,0});
	addOtherVec(iters,std::vector<std::vector<double>>{
			{stopCheckDbl,0,1,0,0,0,0,1},
			{stopCheckDbl,0,2,0,0,0,0,1},
			{stopCheckDbl,0,2,0,0,0,1,1},
			{stopCheckDbl,0,2,0,0,0,2,1},
			{stopCheckDbl,0,2,0,0,0,8,1},
			{stopCheckDbl,0,2,0,0,0,10,2},
			{stopCheckDbl,0,2,0,0,0,12,2},
			{stopCheckDbl,0,2,0,0,0,14,2},
			{stopCheckDbl,0,2,0,0,0,14,2},
			{stopCheckDbl,0,2,0,0,0,14,2}
	});

	for(double hq : iter::range<uint32_t>(1, hqMismatches + 1) ){
		iters.emplace_back(std::vector<double>{stopCheckDbl,0,2,0,0,hq,14,2});
		iters.emplace_back(std::vector<double>{stopCheckDbl,0,2,0,0,hq,14,2});
	}
	iters.emplace_back(std::vector<double>{stopCheckDbl,0,0,0,0,0,0,0});
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genStrictDefaultParsWithHqs(uint32_t stopCheck, uint32_t hqMismatches){
	double stopCheckDbl = stopCheck;
	std::vector<double> noErrorsAll = {stopCheckDbl, 0,0,0,0,0,0,0};
	std::vector<double> noErrorsStopSize3 = {stopCheckDbl, 3,0,0,0,0,0,0};
	std::vector<std::vector<double>> iters;
	iters.push_back(noErrorsStopSize3);
	iters.push_back({stopCheckDbl, 3,1,1,0,0,1,0});
	for(uint32_t hq : iter::range(hqMismatches)){
		iters.push_back({stopCheckDbl, 3,1,1,0,hq + 1.0,hq + 1.0,0});
		iters.push_back({stopCheckDbl, 3,1,1,0,hq + 1.0,hq + 1.0,0});
	}
	iters.push_back(noErrorsStopSize3);
	iters.push_back(noErrorsAll);
	iters.push_back({stopCheckDbl, 0,1,1,0,0,1,0});
	for(uint32_t hq : iter::range(hqMismatches)){
		iters.push_back({stopCheckDbl, 0,1,1,0,hq + 1.0,hq + 1.0,0});
		iters.push_back({stopCheckDbl, 0,1,1,0,hq + 1.0,hq + 1.0,0});
	}
	iters.push_back(noErrorsAll);
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genStrictDefaultParsHq1(uint32_t stopCheck) {
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl, 3,0,0,0,0,0,0},
			{stopCheckDbl, 3,1,1,0,0,1,0},
			{stopCheckDbl, 3,1,1,0,1,1,0},
			{stopCheckDbl, 0,0,0,0,0,0,0},
			{stopCheckDbl, 0,1,1,0,0,1,0},
			{stopCheckDbl, 0,1,1,0,1,1,0},
			{stopCheckDbl, 0,0,0,0,0,0,0}
	};
	CollapseIterations ret;
	ret.addIterations(iters);
	return ret;
}



CollapseIterations CollapseIterations::genOtuPars(uint32_t stopCheck, double perId, bool onHqPerId){
	double stopCheckDbl = stopCheck;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl,3,perId},
			{stopCheckDbl,3,perId},
			{stopCheckDbl,3,perId},
			{stopCheckDbl,0,perId},
			{stopCheckDbl,0,perId},
			{stopCheckDbl,0,perId}
	};
	CollapseIterations ret(true);
	ret.onHqPerId_ = onHqPerId;
	ret.addIterations(iters);
	return ret;
}

CollapseIterations CollapseIterations::genOtu99To97(uint32_t stopCheck, bool onHqPerId){
	double stopCheckDbl = stopCheck;
	double smallestClusToClusSize = 0;
	std::vector<std::vector<double>> iters = {
			{stopCheckDbl,3,.99},
			{stopCheckDbl,3,.99},
			{stopCheckDbl,3,.99},
			{stopCheckDbl,3,.98},
			{stopCheckDbl,3,.98},
			{stopCheckDbl,3,.98},
			{stopCheckDbl,3,.97},
			{stopCheckDbl,3,.97},
			{stopCheckDbl,3,.97},
			{stopCheckDbl,smallestClusToClusSize,.99},
			{stopCheckDbl,smallestClusToClusSize,.99},
			{stopCheckDbl,smallestClusToClusSize,.99},
			{stopCheckDbl,smallestClusToClusSize,.98},
			{stopCheckDbl,smallestClusToClusSize,.98},
			{stopCheckDbl,smallestClusToClusSize,.98},
			{stopCheckDbl,smallestClusToClusSize,.97},
			{stopCheckDbl,smallestClusToClusSize,.97},
			{stopCheckDbl,smallestClusToClusSize,.97}
	};
	CollapseIterations ret(true);
	ret.onHqPerId_ = onHqPerId;
	ret.addIterations(iters);
	return ret;
}


CollapseIterations CollapseIterations::genOtu99(uint32_t stopCheck, bool onHqPerId){
	return genOtuPars(stopCheck, .99, onHqPerId);
}



}  // namespace bibseq

