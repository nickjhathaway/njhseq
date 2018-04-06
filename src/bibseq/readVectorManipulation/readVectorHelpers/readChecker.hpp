#pragma once
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

#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/seqToolsUtils/seqToolsUtils.hpp"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/objects/seqObjects/seqKmers/seqWithKmerInfo.hpp"
#include "bibseq/objects/seqObjects/Paired/PairedRead.hpp"

namespace bibseq {

class ReadChecker {

public:
	ReadChecker(std::string markWith, bool mark);
	const std::string markWith_;
	bool mark_ = true;
	void markName(seqInfo & info) const;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;
	virtual ~ReadChecker();
};

class ReadCheckerLenWithin: public ReadChecker {
public:
	ReadCheckerLenWithin(uint32_t basesWithin, double givenLen, bool mark = true);
	const uint32_t basesWithin_;
	const double givenLen_;
	virtual bool checkRead(seqInfo & info) const;
	virtual ~ReadCheckerLenWithin();
};

class ReadCheckerLenBelow: public ReadChecker {
public:
	ReadCheckerLenBelow(uint32_t maxLen, bool mark = true);
	const uint32_t maxLen_;
	virtual bool checkRead(seqInfo & info) const;

	virtual ~ReadCheckerLenBelow();
};

class ReadCheckerLenAbove: public ReadChecker {
public:
	ReadCheckerLenAbove(uint32_t minLen, bool mark = true);
	const uint32_t minLen_;
	virtual bool checkRead(seqInfo & info) const;
	virtual ~ReadCheckerLenAbove();
};

class ReadCheckerLenBetween: public ReadChecker {
public:
	ReadCheckerLenBetween(uint32_t maxLen, uint32_t minLen, bool mark = true);
	const uint32_t minLen_;
	const uint32_t maxLen_;
	virtual bool checkRead(seqInfo & info) const;

	virtual ~ReadCheckerLenBetween();
};

class ReadCheckerQualCheck: public ReadChecker {
public:
	ReadCheckerQualCheck(uint32_t qualCutOff, double qualFracCutOff, bool mark =
			true);
	const uint32_t qualCutOff_;
	const double qualFracCutOff_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & info) const;

	virtual ~ReadCheckerQualCheck();
};

class ReadCheckerOnCount: public ReadChecker {
public:
	ReadCheckerOnCount(double countCutOff, bool mark = true);
	const double countCutOff_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;
	virtual ~ReadCheckerOnCount();
};

class ReadCheckerOnFrac: public ReadChecker {
public:
	ReadCheckerOnFrac(double fracCutOff, bool mark = true);
	const double fracCutOff_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;
	virtual ~ReadCheckerOnFrac();
};

class ReadCheckerOnNucComp: public ReadChecker {
public:
	ReadCheckerOnNucComp(charCounter counter, double fracDiff, bool mark = true);
	const charCounter counter_;
	const double fracDiff_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;

	virtual ~ReadCheckerOnNucComp();
};

class ReadCheckerOnSeqContaining: public ReadChecker {
public:
	ReadCheckerOnSeqContaining(std::string str, uint32_t occurences, bool mark =
			true);
	const std::string str_;
	const uint32_t occurences_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;
	virtual ~ReadCheckerOnSeqContaining();
};

class ReadCheckerOnNameContaining: public ReadChecker {
public:
	ReadCheckerOnNameContaining(std::string str, bool mark = true);
	const std::string str_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;
	virtual ~ReadCheckerOnNameContaining();
};

class ReadCheckerOnQualityWindow: public ReadChecker {
public:
	ReadCheckerOnQualityWindow(uint32_t qualityWindowLength,
			uint32_t qualityWindowStep, uint32_t qualityWindowThres,
			bool mark = true);
	const uint32_t qualityWindowLength_;
	const uint32_t qualityWindowStep_;
	const uint32_t qualityWindowThres_;
	virtual bool checkRead(seqInfo & info) const;

	virtual ~ReadCheckerOnQualityWindow();
};

class ReadCheckerOnQualityWindowTrim: public ReadChecker {
public:
	ReadCheckerOnQualityWindowTrim(uint32_t qualityWindowLength,
			uint32_t qualityWindowStep, uint32_t qualityWindowThres, uint32_t minLen,
			bool mark = true);
	const uint32_t qualityWindowLength_;
	const uint32_t qualityWindowStep_;
	const uint32_t qualityWindowThres_;
	const uint32_t minLen_;

	virtual bool checkRead(seqInfo & info) const;

	virtual ~ReadCheckerOnQualityWindowTrim();
};

class ReadCheckerOnNs: public ReadCheckerOnSeqContaining {
public:
	ReadCheckerOnNs(bool mark = true);

	virtual ~ReadCheckerOnNs();
};

class ReadCheckerOnKmerComp: public ReadChecker {
public:
	ReadCheckerOnKmerComp(kmerInfo compareInfo, uint32_t kLength,
			double kmerCutoff, bool mark = true);

	const kmerInfo compareInfo_;
	const uint32_t kLength_;
	const double kmerCutoff_;
	virtual bool checkRead(seqInfo & info) const;
	virtual bool checkRead(PairedRead & seq) const;
	virtual ~ReadCheckerOnKmerComp();
};

}  // namespace bibseq

