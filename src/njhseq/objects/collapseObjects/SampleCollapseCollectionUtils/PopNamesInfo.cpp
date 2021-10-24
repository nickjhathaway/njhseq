/*
 * PopNamesInfo.cpp
 *
 *  Created on: Jul 31, 2016
 *      Author: nick
 */
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

#include "PopNamesInfo.hpp"

namespace njhseq {
PopNamesInfo::PopNamesInfo(std::string populationName,
		std::set<std::string> samples, std::set<std::string> controlSamples) :
		populationName_(populationName), samples_(samples), controlSamples_(controlSamples){
	checkPopNameThrow();

}
PopNamesInfo::PopNamesInfo(std::string populationName, VecStr samples, VecStr controlSamples) :
		populationName_(populationName) {
	checkPopNameThrow();
	auto counts = countVec(samples);
	VecStr warnings;
	for (const auto & count : counts) {
		if (count.second > 1) {
			warnings.emplace_back(count.first + ":" + estd::to_string(count.second));
		}
	}
	VecStr controlWarnings;

	for(const auto & controlSamp : controlSamples){
		if(!njh::in(controlSamp, samples)){
			controlWarnings.emplace_back("Control samples should also be in samples, did not find " + controlSamp);
		}
	}
	auto controlCounts = countVec(controlSamples);
	for (const auto & count : controlCounts) {
		if (count.second > 1) {
			controlWarnings.emplace_back("Found multiple of sample: "  + count.first + ":" + estd::to_string(count.second));
		}
	}
	if (!warnings.empty() || !controlWarnings.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": Error, the following samples were entered more than once\n";
		for (const auto & warn : warnings) {
			ss << warn << std::endl;
		}
		if(!controlWarnings.empty()){
			ss << "Control Samples Warnings" << std::endl;
		}
		for(const auto & warn : controlWarnings){
			ss << warn << std::endl;
		}
		throw std::runtime_error { ss.str() };
	}
	samples_ = std::set<std::string> { samples.begin(), samples.end() };
	controlSamples_ = std::set<std::string> { controlSamples.begin(), controlSamples.end() };
}

void PopNamesInfo::checkPopNameThrow()const{
	if(njh::containsSubString(populationName_, ".")){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, population name can't have a period in it, " << populationName_ << std::endl;
		throw std::runtime_error{ss.str()};
	}
}

bool PopNamesInfo::hasSample(const std::string & sample) const {
	return njh::has(samples_, sample);
}


Json::Value PopNamesInfo::toJson() const{
	Json::Value ret;
	ret["class"] =           njh::json::toJson(njh::getTypeName(*this));
	ret["samples_"] =        njh::json::toJson(samples_);
	ret["controlSamples_"] =        njh::json::toJson(controlSamples_);
	ret["populationName_"] = njh::json::toJson(populationName_);
	return ret;
}


}  // namespace njhseq

