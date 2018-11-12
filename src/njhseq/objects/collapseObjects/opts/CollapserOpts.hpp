#pragma once
/*
 * collapserOpts.hpp
 *
 *  Created on: Jul 30, 2015
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
#include "njhseq/common.h"
#include "njhseq/alignment/alignerUtils/comparison.hpp"

namespace njhseq {

struct ChimeraOpts {
	ChimeraOpts();
	double parentFreqs_ = 2;
	uint32_t runCutOff_ = 1;
	comparison chiOverlap_;
	uint32_t overLapSizeCutoff_ = 5;
	bool checkChimeras_ = false;
	bool keepLowQaulityMismatches_ = false;

};


struct NucCompBinOpts {
	NucCompBinOpts(bool useNucComp, bool useMinLenNucComp, bool findBestNuc,
			const std::vector<double> & diffCutOffVec);
	NucCompBinOpts();
	bool useNucComp_ = false;
	bool useMinLenNucComp_ = true;
	bool findBestNuc_ = true;
	std::vector<double> diffCutOffVec_;
};

struct KmerBinOpts {
	KmerBinOpts(bool useKmerBinning, uint32_t kCompareLen, double kmerCutOff);
	KmerBinOpts();
	bool useKmerBinning_ = false;
	uint32_t kCompareLen_ = 10;
	double kmerCutOff_ = 0.80;

};

struct KmerOpts {
	KmerOpts(bool checkKmers, bool kmersByPosition,const std::string & runCutOffString, uint32_t runCutOff,
			uint32_t kLength);
	KmerOpts();
	bool checkKmers_ = true;
	bool kmersByPosition_ = true;
	std::string runCutOffString_ = "1";
	uint32_t runCutOff_ = 1;
	uint32_t kLength_ = 9;
};

struct BestMatchOpts {
	BestMatchOpts(bool findingBestMatch, uint32_t bestMatchCheck);
	BestMatchOpts();
	bool findingBestMatch_ = true;
	uint32_t bestMatchCheck_ = 10;

};

struct ITOpts {
	ITOpts(bool adjustHomopolyerRuns, bool weighHomopolyer,
			bool removeLowQualityBases, uint32_t lowQualityBaseTrim);
	ITOpts();
	bool adjustHomopolyerRuns_ = false;
	bool weighHomopolyer_ = false;
	bool removeLowQualityBases_ = false;
	uint32_t lowQualityBaseTrim_ = 3;
};

struct SkipOpts {
	SkipOpts(bool skipOnLetterCounterDifference, double fractionDifferenceCutOff,
			bool useReadLen, uint32_t readLenDiff);
	SkipOpts();
	bool skipOnLetterCounterDifference_ = false;
	double fractionDifferenceCutOff_ = 0.05;
	bool useReadLen_ = false;
	uint32_t readLenDiff_ = 15;
};

struct VerboseOpts {
	VerboseOpts(bool verbose, bool debug) ;
	VerboseOpts();
	bool verbose_ = false;
	bool debug_ = false;
};

struct AlignOpts {
	AlignOpts(bool eventBased, bool countEndGaps,  bool noAlign);
	AlignOpts();
	bool eventBased_ = true;
	bool countEndGaps_ = false;
	bool noAlign_ = false;
};

struct ClusteringOpts{
	ClusteringOpts(bool converge);
	ClusteringOpts();
	bool converge_ = false;
};


class CollapserOpts {
public:
	CollapserOpts(KmerOpts kmerOpts, BestMatchOpts bestMatchOpts, ITOpts iTOpts,
			SkipOpts skipOpts, VerboseOpts verboseOpts, AlignOpts alignOpts,
			NucCompBinOpts nucCompBinOpts, KmerBinOpts kmerBinOpts, ClusteringOpts clusOpts);
	CollapserOpts();


	KmerOpts kmerOpts_;
	BestMatchOpts bestMatchOpts_;
	ITOpts iTOpts_;
	SkipOpts skipOpts_;
	VerboseOpts verboseOpts_;
	AlignOpts alignOpts_;
	NucCompBinOpts nucCompBinOpts_;
	KmerBinOpts kmerBinOpts_;
	ClusteringOpts clusOpts_;
};

} /* namespace njhseq */


