#pragma once
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


#include "njhseq/objects/seqObjects/BaseObjects/baseReadObject.hpp"
#include "njhseq/objects/kmer/kmerMap.hpp"
#include "njhseq/objects/helperObjects/tandemRepeat.hpp"
#include "njhseq/alignment/alignerUtils.h"
#include "njhseq/alignment/alnCache/alnInfoHolder.hpp"
#include "njhseq/alignment/aligner/alnParts.hpp"

namespace njhseq {

/*! \brief Aligner Class
 *
 *
 *  This class can do local or global alignment using simple scoring or using a
 *  provide scoring matrix
 */
class aligner {

 public:
  // constructors
	/**@brief Default aligner, can handle alignments of up to 400 bps
	 */
  aligner();

	aligner(uint64_t maxSize, const gapScoringParameters& gapPars);

	aligner(uint64_t maxSize, const gapScoringParameters& gapPars,
          const substituteMatrix& scoreMatrix);

	aligner(uint64_t maxSize, const gapScoringParameters& gapPars,
          const substituteMatrix& scoreMatrix, bool countEndGaps);

	aligner(uint64_t maxSize, const gapScoringParameters & gapPars,
			const substituteMatrix& subMatrix, const KmerMaps& kmaps,
			QualScorePars qScorePars, bool countEndGaps,
			bool weighHomopolymers);



	void setGapScoring(const gapScoringParameters & gapPars);
	//void setMatchScoring(const substituteMatrix& subMatrix);

  // to hold the sequence alignments
  baseReadObject alignObjectA_;
  baseReadObject alignObjectB_;

  alnParts parts_;
  uint32_t inputAlignmentBlockSize_{100};
  uint32_t inputAlignmentBlockWalkbackSize_{50};

  alnInfoMasterHolder alnHolder_;

  uint32_t numberOfAlingmentsDone_ = 0;

  comparison comp_;
  KmerMaps kMaps_;

  QualScorePars qScorePars_;
  bool countEndGaps_;
  bool weighHomopolymers_;

  void resetAlnCache();

  void processAlnInfoInput(const std::string& alnInfoDirName, bool verbose = false);
	void processAlnInfoOutput(const std::string& outAlnInfoDirName, bool verbose);

  void processAlnInfoInputNoCheck(const std::string& alnInfoDirName, bool verbose);
	void processAlnInfoOutputNoCheck(const std::string& outAlnInfoDirName, bool verbose);


  // Aligner
	void alignScoreLocal(const std::string& firstSeq, const std::string& secondSeq);
	void alignScoreCacheLocal(const std::string& firstSeq,
			const std::string& secondSeq);
	void alignScoreGlobal(const std::string& firstSeq, const std::string& secondSeq);
	void alignScoreGlobalDiag(const std::string& firstSeq, const std::string& secondSeq);
	void alignScoreGlobalNoInternalGaps(const std::string& firstSeq, const std::string& secondSeq);

	void alignScoreCacheGlobal(const std::string& firstSeq,
			const std::string& secondSeq);
	void alignScoreCacheGlobalDiag(const std::string& firstSeq,
			const std::string& secondSeq);

	void alignScore(const std::string& firstSeq, const std::string& secondSeq,
			bool local);

	void alignScoreCache(const std::string& firstSeq,
			const std::string& secondSeq, bool local);

	void alignCacheLocal(const seqInfo & ref, const seqInfo & read);

	template<typename READ1, typename READ2>
	void alignCacheLocal(const READ1 & ref, const READ2 & read){
		alignCacheLocal(getSeqBase(ref), getSeqBase(read));
	}

	void alignCacheGlobal(const seqInfo & ref, const seqInfo & read);

	template<typename READ1, typename READ2>
	void alignCacheGlobal(const READ1 & ref, const READ2 & read){
		alignCacheGlobal(getSeqBase(ref), getSeqBase(read));
	}

	void alignCache(const seqInfo & ref, const seqInfo & read, bool local);

	template<typename READ1, typename READ2>
	void alignCache(const READ1 & ref, const READ2 & read, bool local){
		alignCache(getSeqBase(ref), getSeqBase(read), local);
	}

	template<typename READ1, typename READ2>
	void alignReg(const READ1 & ref, const READ2 & read, bool local){
		alignReg(getSeqBase(ref), getSeqBase(read), local);
	}
	void alignReg(const seqInfo & ref, const seqInfo & read, bool local);

	template<typename READ1, typename READ2>
	void alignRegGlobal(const READ1 & ref, const READ2 & read){
		alignRegGlobal(getSeqBase(ref), getSeqBase(read));
	}
	void alignRegGlobal(const seqInfo & ref, const seqInfo & read);



	void alignCacheGlobalDiag(const std::string & ref, const std::string & read);
	void alignCacheGlobalDiag(const seqInfo & ref, const seqInfo & read);
	template<typename READ1, typename READ2>
	void alignCacheGlobalDiag(const READ1 & ref, const READ2 & read){
		alignCacheGlobalDiag(getSeqBase(ref), getSeqBase(read));
	}

	void alignRegGlobalDiag(const std::string & ref, const std::string & read);
	void alignRegGlobalDiag(const seqInfo & ref, const seqInfo & read);
	template<typename READ1, typename READ2>
	void alignRegGlobalDiag(const READ1 & ref, const READ2 & read){
		alignRegGlobalDiag(getSeqBase(ref), getSeqBase(read));
	}





	template<typename READ1, typename READ2>
	void alignRegGlobalNoInternalGaps(const READ1 & ref, const READ2 & read){
		alignRegGlobalNoInternalGaps(getSeqBase(ref), getSeqBase(read));
	}
	void alignRegGlobalNoInternalGaps(const seqInfo & ref, const seqInfo & read);


	template<typename READ1, typename READ2>
	void alignRegLocal(const READ1 & ref, const READ2 & read){
		alignRegLocal(getSeqBase(ref), getSeqBase(read));
	}
	void alignRegLocal(const seqInfo & ref, const seqInfo & read);

	std::pair<uint32_t, uint32_t> findReversePrimer(const std::string& read,
			const std::string& primer);
	std::pair<uint32_t, uint32_t> findReversePrimer(const baseReadObject& read,
			const baseReadObject& primer);

	void rearrangeSeq(const std::string& firstRead, const std::string& secondRead,
			bool local);

	void rearrangeObjsLocal(const seqInfo& firstRead, const seqInfo& secondRead);
	template<typename READ1, typename READ2>
	void rearrangeObjsLocal(const READ1 & ref, const READ2 & read){
		rearrangeObjsLocal(getSeqBase(ref), getSeqBase(read));
	}

	void rearrangeObjsGlobal(const std::string & firstRead, const std::string& secondRead);
	void rearrangeObjsGlobal(const seqInfo& firstRead, const seqInfo& secondRead);
	template<typename READ1, typename READ2>
	void rearrangeObjsGlobal(const READ1 & ref, const READ2 & read){
		rearrangeObjsGlobal(getSeqBase(ref), getSeqBase(read));
	}

	void rearrangeObjs(const seqInfo& firstRead, const seqInfo& secondRead,bool local);
	template<typename READ1, typename READ2>
	void rearrangeObjs(const READ1 & ref, const READ2 & read, bool local){
		rearrangeObjs(getSeqBase(ref), getSeqBase(read), local);
	}

  void noAlignSetAndScore(const seqInfo& objectA,
                          const seqInfo& objectB);
	template<typename READ1, typename READ2>
	void noAlignSetAndScore(const READ1 & ref, const READ2 & read){
		noAlignSetAndScore(getSeqBase(ref), getSeqBase(read));
	}

  void scoreAlignment(bool editTheSame);
  bool CountEndGaps() { return countEndGaps_; }
  void setQual(QualScorePars qScorePars);

  // Outputting
  void outPutParameterInfo(std::ostream& out) const;


  void handleGapCountingInA(gap& currentGap);
  void handleGapCountingInB(gap& currentGap);

	const comparison & profileAlignment(const seqInfo& objectA,
			const seqInfo& objectB, bool checkKmer, bool usingQuality,
			bool doingMatchQuality, uint32_t start = 0, uint32_t stop = 0);

	template<typename READ1, typename READ2>
	const comparison & profileAlignment(const READ1& objectA,
			const READ2& objectB, bool checkKmer, bool usingQuality,
			bool doingMatchQuality, uint32_t start = 0, uint32_t stop = 0) {
		return profileAlignment(getSeqBase(objectA), getSeqBase(objectB), checkKmer,
				usingQuality, doingMatchQuality, start, stop);
	}


  const comparison & profilePrimerAlignment(const std::string& objectA,
			const std::string& objectB);
  const comparison & profilePrimerAlignment(const seqInfo& objectA,
			const seqInfo& objectB);
	template<typename READ1, typename READ2>
	const comparison & profilePrimerAlignment(const READ1& objectA,
			const READ2& objectB) {
		return profilePrimerAlignment(getSeqBase(objectA), getSeqBase(objectB));
	}

	comparison compareAlignment(const seqInfo& objectA, const seqInfo& objectB,
			bool checkKmers);

	template<typename READ1, typename READ2>
	comparison compareAlignment(const READ1& objectA, const READ2& objectB,
			bool checkKmers) {
		return compareAlignment(getSeqBase(objectA), getSeqBase(objectB),
				checkKmers);
	}





  void setDefaultQualities();
	void resetCounts();
	void resetAlignmentInfo();

  struct checkForTandemRepeatGapPars {
	  checkForTandemRepeatGapPars() {

	  }
	  int32_t match{2};
	  int32_t mismatch{-2};
	  int32_t gap{-7};
	  int32_t minimumAlignScore{16};
  };

	bool checkForTandemRepeatGap(const checkForTandemRepeatGapPars & pars = checkForTandemRepeatGapPars() ) const;
	std::map<uint32_t, gap> getTandemRepeatGapsForCurrentAlignment(const checkForTandemRepeatGapPars & pars = checkForTandemRepeatGapPars() ) const;

	static bool checkTwoEqualSeqs(const std::string& seq1,
			const std::string& seq2, int allowableMismatches);
	// finding tandem repeats
	static std::vector<TandemRepeat> findTandemRepeatsInSequence(
			const std::string& str, int match = 2, int mismatch = -2, int gap = -7,
			int minimumAlignScore = 50);
	static TandemRepeat findTandemRepeatOfStrInSequence(const std::string & str,
			const std::string & tandem, int match = 2, int mismatch = -2, int gap = -7,
			int minimumAlignScore = 50);
	static TandemRepeat findTandemRepeatOfStrInSequenceDegen(std::string str,
			std::string tandem, int match = 2, int mismatch = -2, int gap = -7,
			int minimumAlignScore = 50);
	static bool checkTwoStringsDegen(const std::string& str1,
			const std::string& str2, int allowableMismatches,
			const substituteMatrix& scoringArray);
	size_t getAlignPosForSeqAPos(size_t seqAPos);
	size_t getAlignPosForSeqBPos(size_t seqBPos);
	size_t getSeqPosForAlnAPos(size_t alnAPos);
	size_t getSeqPosForAlnBPos(size_t alnBPos);

 public:
  void setGeneralScorring(int32_t generalMatch, int32_t generalMismatch);

};



}  // namespace njhseq


