//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "readChecker.hpp"



namespace bibseq {

bool readChecker::checkReadLenWithin(seqInfo & info, uint32_t basesWithin,
		double givenLen, bool mark) {
	std::string markWith = "_lenMoreThan" + std::to_string(basesWithin) +
      "_From" + std::to_string(givenLen);
	if (std::abs(givenLen - info.seq_.length()) > basesWithin) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadLenBellow(seqInfo & info, uint32_t maxLen,
		bool mark) {
	std::string markWith = "_length>" + std::to_string(maxLen);
	if (info.seq_.length() > maxLen) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadLenAbove(seqInfo & info, uint32_t minLen,
		bool mark) {
	std::string markWith = "_length<" + std::to_string(minLen);
	if (info.seq_.length() < minLen) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadLenBetween(seqInfo & info, uint32_t maxLen,
		uint32_t minLen, bool mark) {
	if(checkReadLenBellow(info,maxLen, mark) && checkReadLenAbove(info, minLen, mark)){
		return true;
	}
	return false;
}


bool readChecker::checkReadQualCheck(seqInfo & info, uint32_t qualCutOff,
		double qualFracCutOff, bool mark) {

	std::string markWith = "_q" + to_string(qualCutOff) + "<"
			+ std::to_string(qualFracCutOff);

	if (info.getQualCheck(qualCutOff) < qualFracCutOff) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnCount(seqInfo & info, double countCutOff,
		bool mark) {
	std::string markWith = "_cnt<=" + estd::to_string(countCutOff);
	if (info.cnt_ <= countCutOff) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnFrac(seqInfo & info, double fracCutOff,
		bool mark) {
	std::string markWith = "_frac<=" + estd::to_string(fracCutOff);
	if (info.frac_ <= fracCutOff) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnNucComp(seqInfo & info, const charCounterArray & counter, double fracDiff,
		bool mark) {
	std::string markWith = "_nucCompFrac>" + estd::to_string(fracDiff);
	charCounterArray count(counter.alphabet_);
	count.increaseCountByString(info.seq_);
	count.setFractions();
	if (counter.getFracDifference(count, counter.alphabet_) > fracDiff) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnSeqContaining(seqInfo & info, const std::string& str,
		int occurences, bool mark) {
	std::string markWith = "_seqContains" + estd::to_string(str);
	std::string lowerCaseSearch = stringToLowerReturn(str);
	if (countOccurences(info.seq_, str)
			+ countOccurences(info.seq_, lowerCaseSearch) >= occurences) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnNameContaining(seqInfo & info, const std::string& str,
		bool mark) {
	std::string markWith = "_nameContains" + estd::to_string(str);
	if (info.name_.find(str) != std::string::npos) {
		info.on_ = false;
		if(mark) info.name_.append(markWith);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnQualityWindow(seqInfo & info, int qualityWindowLength,
		int qualityWindowStep, int qualityWindowThres, bool mark) {
	std::stringstream window;
	window << "_failedWindowOf_wl:" << qualityWindowLength << ",ws:"
			<< qualityWindowStep << ",wt:" << qualityWindowThres;
	if (!seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
			qualityWindowStep, info.qual_)) {
		info.on_ = false;
		if(mark) info.name_.append(window.str());
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnQualityWindowTrim(seqInfo & info, int qualityWindowLength,
    int qualityWindowStep, int qualityWindowThres, uint32_t minLen, bool mark) {
	std::stringstream window;
	window << "_failedWindowOf_wl:" << qualityWindowLength << ",ws:"
			<< qualityWindowStep << ",wt:" << qualityWindowThres;
	auto pos = seqUtil::checkQualityWindowPos(
      qualityWindowLength, qualityWindowThres, qualityWindowStep,
      info.qual_);
	if (pos <= minLen) {
		info.on_ = false;
		if(mark) info.name_.append(window.str());
		return false;
	} else {
		info.setClip(0, pos - 1);
		info.on_ = true;
		return true;
	}
}

bool readChecker::checkReadOnNs(seqInfo & info, bool mark){
	return checkReadOnSeqContaining(info, "N", 1, mark);
}

bool readChecker::checkReadOnKmerComp(seqInfo & info,kmerInfo compareInfo, uint32_t kLength, double kmerCutoff, bool mark) {
	std::stringstream window;
	window << "_failedKmerCutOff:" << kmerCutoff << "_kLen_" << kLength;
	kmerInfo currentInfo(info.seq_, kLength);
	auto dist = compareInfo.compareKmers(currentInfo);
	if(dist.second < kmerCutoff){
		info.on_ = false;
		if(mark) info.name_.append(window.str());
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

}  // namespace bibseq
