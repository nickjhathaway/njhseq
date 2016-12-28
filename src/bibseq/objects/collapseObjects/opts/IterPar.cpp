#include "IterPar.hpp"

#include "bibseq/objects/dataContainers/tables/table.hpp"
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
#include <bibcpp/bashUtils.h>

namespace bibseq {

bool CompareIDProfile::passErrors(const comparison & threshold,
		const comparison & generated) {
	return threshold.passIdThreshold(generated);
}

CompareIDProfile::~CompareIDProfile() {
}

bool CompareIDHqProfile::passErrors(const comparison & threshold,
		const comparison & generated) {
	return threshold.passIdThresholdHq(generated);
}

CompareIDHqProfile::~CompareIDHqProfile() {

}



bool CompareErrorProfile::passErrors(const comparison & threshold,
		const comparison & generated) {
	return threshold.passErrorProfile(generated);
}
CompareErrorProfile::~CompareErrorProfile() {
}

bool CompareIDAndErrorProfile::passErrors(const comparison & threshold,
		const comparison & generated) {
	return threshold.passIdAndErrorThreshold(generated);
}
CompareIDAndErrorProfile::~CompareIDAndErrorProfile() {
}

namespace cto = bib::bashCT;

IterPar& IterPar::operator=(const IterPar & other) {
	//std::cout << "copy constructor of IterPar\n";
	stopCheck_ = other.stopCheck_;
	smallCheckStop_ = other.smallCheckStop_;
	iterNumber_ = other.iterNumber_;
	onPerId_ = other.onPerId_;
	onHqPerId_ = other.onHqPerId_;
	errors_ = other.errors_;

	setCompFunc();
	return *this;
}

IterPar& IterPar::operator=(IterPar&& other) {
	//std::cout << "move assignment of IterPar\n";
	stopCheck_ = other.stopCheck_;
	smallCheckStop_ = other.smallCheckStop_;
	iterNumber_ = other.iterNumber_;
	onPerId_ = other.onPerId_;
	onHqPerId_ = other.onHqPerId_;
	errors_ = other.errors_;
	setCompFunc();
	return *this;
}
/*
 * 	// procedure parameters
	uint32_t stopCheck_;
	uint32_t smallCheckStop_;
	uint32_t iterNumber_;
	bool onPerId_ = false;
	// error parameters
	comparison errors_;
 */

IterPar::IterPar(const IterPar& other) :
		stopCheck_(other.stopCheck_), smallCheckStop_(other.smallCheckStop_), iterNumber_(
				other.iterNumber_), onPerId_(other.onPerId_), onHqPerId_(
				other.onHqPerId_), errors_(other.errors_) {
	//std::cout << "copy assignment of IterPar\n";
	setCompFunc();
}




void IterPar::setCompFunc() {
	if (onPerId_) {
		if (onHqPerId_) {
			comp_ = std::make_shared<CompareIDHqProfile>();
		} else {
			comp_ = std::make_shared<CompareIDProfile>();
		}
	} else {
		comp_ = std::make_shared<CompareErrorProfile>();
	}
}

bool IterPar::passErrorCheck(const comparison & generated) const{
	return comp_->passErrors(errors_, generated);
}

IterPar::IterPar() :
		stopCheck_(100), smallCheckStop_(0), iterNumber_(0) {
	setCompFunc();
}


IterPar::IterPar(const std::vector<double>& parameter,
		uint32_t iterNumber,
		const PerIdPars & perIdPars):iterNumber_(iterNumber),
				onPerId_(perIdPars.onPerId_),
				onHqPerId_(perIdPars.onHqPerId_) {
	//check size
	if (onPerId_) {
		if (parameter.size() != 3) {
			std::stringstream ss;
			ss
					<< cto::boldRed(
							"Error in constructing runningParameters, size should be 3, not ")
					<< parameter.size() << " when constructing a percent identity parameters" << std::endl;
			throw std::runtime_error { ss.str() };
		}
	} else {
		if (parameter.size() != 8) {
			std::stringstream ss;
			ss
					<< cto::boldRed(
							"Error in constructing runningParameters, size should be 8, not ")
					<< parameter.size() << " check input parameters file or check to see if"
							" --onPerId should be on if the input parameters file is for otu clustering"<< std::endl;
			throw std::runtime_error { ss.str() };
		}
	}
  stopCheck_ = parameter[0];
  smallCheckStop_ = parameter[1];

	if (onPerId_) {
	  // error parameters
	  errors_.distances_.eventBasedIdentity_ = parameter[2];
	  errors_.distances_.eventBasedIdentityHq_ = parameter[2];
	} else {
	  // error parameters
	  errors_.oneBaseIndel_ = parameter[2];
	  errors_.twoBaseIndel_ = parameter[3];
	  errors_.largeBaseIndel_ = parameter[4];

	  errors_.hqMismatches_ = parameter[5];
	  errors_.lqMismatches_ = parameter[6];
	  errors_.lowKmerMismatches_ = parameter[7];
	}
	setCompFunc();
}



std::string centerText(const std::string& text, uint32_t maxWidth) {
	uint32_t spaceWidth = maxWidth - text.size();
	uint32_t front = round(spaceWidth/2.0);
	uint32_t back = spaceWidth - front;
  return std::string(front, ' ') + text + std::string(back, ' ');
}

void IterPar::printIterInfo(std::ostream & out, bool colorFormat) const{
	if(onPerId_){
		std::cout << cto::bold
				<< cto::addBGColor(145) << " Iteration "
				<< cto::addBGColor(188) << " StopAfter "
				<< cto::addBGColor(145) << " SizeCutOff "
				<< cto::addBGColor(188);
		if(onHqPerId_){
			std::cout << " HqPerIdCutOff "<< cto::reset << std::endl;
		}else{
			std::cout << " PerIdCutOff "<< cto::reset << std::endl;
		}

		std::cout
		<< cto::addBGColor(145) << " " << centerText(estd::to_string(iterNumber_), 9) << " "
		<< cto::addBGColor(188) << " " << centerText(estd::to_string(stopCheck_), 9) << " "
		<< cto::addBGColor(145) << " " << centerText(estd::to_string(smallCheckStop_), 10) << " ";
		if(onHqPerId_){
			std::cout
			<< cto::addBGColor(188) << " " << centerText(estd::to_string(errors_.distances_.eventBasedIdentity_ * 100) + "%", 13) << " ";
		}else{
			std::cout
			<< cto::addBGColor(188) << " " << centerText(estd::to_string(errors_.distances_.eventBasedIdentity_ * 100) + "%", 11) << " ";
		}
		std::cout << cto::reset << std::endl;
	}else{
		std::cout << cto::bold
				<< cto::addBGColor(145) << " Iteration "
				<< cto::addBGColor(188) << " StopAfter "
				<< cto::addBGColor(145) << " SizeCutOff "
				<< cto::addBGColor(188) << " 1BaseIndels "
				<< cto::addBGColor(145) << " 2BaseIndels "
				<< cto::addBGColor(188) << " >2BaseIndels "
				<< cto::addBGColor(145) << " HQMismatches "
				<< cto::addBGColor(188) << " LQMismatches "
				<< cto::addBGColor(145) << " LKMismatches "  << cto::reset << std::endl;
		std::cout
		<< cto::addBGColor(145) << " " << centerText(estd::to_string(iterNumber_), 9) << " "
		<< cto::addBGColor(188) << " " << centerText(estd::to_string(stopCheck_), 9) << " "
		<< cto::addBGColor(145) << " " << centerText(estd::to_string(smallCheckStop_), 10) << " "
		<< cto::addBGColor(188) << " " << centerText(errors_.oneBaseIndel_ == std::numeric_limits<uint32_t>::max() ? "all" : estd::to_string(errors_.oneBaseIndel_), 11) << " "
		<< cto::addBGColor(145) << " " << centerText(errors_.twoBaseIndel_ == std::numeric_limits<uint32_t>::max() ? "all" : estd::to_string(errors_.twoBaseIndel_), 11) << " "
		<< cto::addBGColor(188) << " " << centerText(errors_.largeBaseIndel_ == std::numeric_limits<uint32_t>::max() ? "all" : estd::to_string(errors_.largeBaseIndel_), 12) << " "
		<< cto::addBGColor(145) << " " << centerText(errors_.hqMismatches_ == std::numeric_limits<uint32_t>::max() ? "all" : estd::to_string(errors_.hqMismatches_), 12) << " "
		<< cto::addBGColor(188) << " " << centerText(errors_.lqMismatches_ == std::numeric_limits<uint32_t>::max() ? "all" : estd::to_string(errors_.lqMismatches_), 12) << " "
		<< cto::addBGColor(145) << " " << centerText(errors_.lowKmerMismatches_ == std::numeric_limits<uint32_t>::max() ? "all" : estd::to_string(errors_.lowKmerMismatches_), 12) << " "
		<< cto::reset << std::endl;
	}

}

VecStr IterPar::outputIterInfoHeader(bool onPerId) {
	if (onPerId) {
		return VecStr { "stopCheck", "smallCutoff", "id" };
	} else {
		return VecStr { "stopCheck", "smallCutoff", "1baseIndel", "2baseIndel",
				">2baseIndel", "HQMismatches", "LQMismatches", "LKMismatches" };
	}
}

VecStr IterPar::outputIterInfo() const {
	if (onPerId_) {
		return toVecStr(stopCheck_, smallCheckStop_,
				errors_.distances_.eventBasedIdentity_);
	} else {
		return toVecStr(stopCheck_, smallCheckStop_, errors_.oneBaseIndel_,
				errors_.twoBaseIndel_, errors_.largeBaseIndel_, errors_.hqMismatches_,
				errors_.lqMismatches_, errors_.lowKmerMismatches_);
	}
}

}  // namespace bib
