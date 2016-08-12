/*
 * collapserOpts.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: nick
 */

#include "CollapserOpts.hpp"

namespace bibseq {

ChimeraOpts::ChimeraOpts(){
	chiOverlap_.oneBaseIndel_ = 2;
	chiOverlap_.twoBaseIndel_ = 1;
	chiOverlap_.lowKmerMismatches_ = 1;
}

NucCompBinOpts::NucCompBinOpts(){

}
KmerBinOpts::KmerBinOpts(){

}


NucCompBinOpts::NucCompBinOpts(bool useNucComp, bool useMinLenNucComp,
		bool findBestNuc, const std::vector<double> & diffCutOffVec) :
		useNucComp_(useNucComp), useMinLenNucComp_(useMinLenNucComp), findBestNuc_(
				findBestNuc), diffCutOffVec_(diffCutOffVec) {
}

KmerBinOpts::KmerBinOpts(bool useKmerBinning, uint32_t kCompareLen,
		double kmerCutOff) :
		useKmerBinning_(useKmerBinning), kCompareLen_(kCompareLen), kmerCutOff_(
				kmerCutOff) {
}

KmerOpts::KmerOpts(bool checkKmers, bool kmersByPosition,
		const std::string & runCutOffString, uint32_t runCutOff, uint32_t kLength) :
		checkKmers_(checkKmers), kmersByPosition_(kmersByPosition), runCutOffString_(
				runCutOffString), runCutOff_(runCutOff), kLength_(kLength) {
}

KmerOpts::KmerOpts(){}

BestMatchOpts::BestMatchOpts(bool findingBestMatch, uint32_t bestMatchCheck) :
		findingBestMatch_(findingBestMatch), bestMatchCheck_(bestMatchCheck) {
}

BestMatchOpts::BestMatchOpts(){}

ITOpts::ITOpts(bool adjustHomopolyerRuns, bool weighHomopolyer,
		bool removeLowQualityBases, uint32_t lowQualityBaseTrim) :
		adjustHomopolyerRuns_(adjustHomopolyerRuns), weighHomopolyer_(
				weighHomopolyer), removeLowQualityBases_(removeLowQualityBases), lowQualityBaseTrim_(
				lowQualityBaseTrim) {
}

ITOpts::ITOpts(){};


SkipOpts::SkipOpts(bool skipOnLetterCounterDifference, double fractionDifferenceCutOff,
		bool useReadLen, uint32_t readLenDiff) :
		skipOnLetterCounterDifference_(skipOnLetterCounterDifference), fractionDifferenceCutOff_(
				fractionDifferenceCutOff), useReadLen_(useReadLen), readLenDiff_(
				readLenDiff) {
}

SkipOpts::SkipOpts(){}

VerboseOpts::VerboseOpts(bool verbose, bool debug) :
		verbose_(verbose), debug_(debug) {
}

VerboseOpts::VerboseOpts() {
}

AlignOpts::AlignOpts(bool eventBased, bool countEndGaps, bool noAlign) :
		eventBased_(eventBased), countEndGaps_(countEndGaps), noAlign_(noAlign) {
}

AlignOpts::AlignOpts() {}

CollapserOpts::CollapserOpts(KmerOpts kmerOpts, BestMatchOpts bestMatchOpts, ITOpts iTOpts,
		SkipOpts skipOpts, VerboseOpts verboseOpts, AlignOpts alignOpts, NucCompBinOpts nucCompBinOpts, KmerBinOpts kmerBinOpts) :
		kmerOpts_(kmerOpts), bestMatchOpts_(bestMatchOpts), iTOpts_(iTOpts), skipOpts_(
				skipOpts), verboseOpts_(verboseOpts), alignOpts_(alignOpts), nucCompBinOpts_(nucCompBinOpts),
				kmerBinOpts_(kmerBinOpts){
}

CollapserOpts::CollapserOpts(){}

} /* namespace bibseq */
