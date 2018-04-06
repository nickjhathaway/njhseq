/*
 * PopNamesInfo.cpp
 *
 *  Created on: Jul 31, 2016
 *      Author: nick
 */
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

#include "PopNamesInfo.hpp"

namespace bibseq {
PopNamesInfo::PopNamesInfo(std::string populationName,
		std::set<std::string> samples) :
		populationName_(populationName), samples_(samples) {
	checkPopNameThrow();

}
PopNamesInfo::PopNamesInfo(std::string populationName, VecStr samples) :
		populationName_(populationName) {
	checkPopNameThrow();
	auto counts = countVec(samples);
	VecStr warnings;
	for (const auto & count : counts) {
		if (count.second > 1) {
			warnings.emplace_back(count.first + ":" + estd::to_string(count.second));
		}
	}
	if (!warnings.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": Error, the following samples were entered more than once\n";
		for (const auto & warn : warnings) {
			ss << warn << std::endl;
		}
		throw std::runtime_error { ss.str() };
	}
	samples_ = std::set<std::string> { samples.begin(), samples.end() };
}

void PopNamesInfo::checkPopNameThrow()const{
	if(bib::containsSubString(populationName_, ".")){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, population name can't have a period in it, " << populationName_ << std::endl;
		throw std::runtime_error{ss.str()};
	}
}

bool PopNamesInfo::hasSample(const std::string & sample) const {
	return bib::has(samples_, sample);
}


Json::Value PopNamesInfo::toJson() const{
	Json::Value ret;
	ret["class"] =           bib::json::toJson(bib::getTypeName(*this));
	ret["samples_"] =        bib::json::toJson(samples_);
	ret["populationName_"] = bib::json::toJson(populationName_);
	return ret;
}


}  // namespace bibseq

