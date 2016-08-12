#pragma once
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


#include "bibseq/objects/seqObjects/BaseObjects/baseReadObject.hpp"
#include "bibseq/objects/kmer/kmerMap.hpp"
#include "bibseq/objects/helperObjects/tandemRepeat.hpp"
#include "bibseq/alignment/alignerUtils.h"
#include "bibseq/alignment/alnCache/alnInfoHolder.hpp"
#include "bibseq/alignment/aligner/alnParts.hpp"

namespace bibseq {

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

	aligner(uint64_t maxSize, const gapScoringParameters& gapPars,
          const substituteMatrix& scoreMatrix);

	aligner(uint64_t maxSize, const gapScoringParameters& gapPars,
          const substituteMatrix& scoreMatrix, bool countEndGaps);

	aligner(uint64_t maxSize, const gapScoringParameters & gapPars,
			const substituteMatrix& subMatrix, const KmerMaps& kmaps,
			QualScorePars qScorePars, bool countEndGaps,
			bool weighHomopolymers);



  // to hold the sequence alignments
  baseReadObject alignObjectA_;
  baseReadObject alignObjectB_;

  alnParts parts_;
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
	void alignScoreCacheGlobal(const std::string& firstSeq,
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


	void alignReg(const baseReadObject & ref, const baseReadObject & read,
			bool local);
	void alignReg(const seqInfo & ref, const seqInfo & read, bool local);

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

	void rearrangeObjsGlobal(const seqInfo& firstRead, const seqInfo& secondRead);
	template<typename READ1, typename READ2>
	void rearrangeObjsGlobal(const READ1 & ref, const READ2 & read){
		rearrangeObjsGlobal(getSeqBase(ref), getSeqBase(read));
	}

	void rearrangeObjs(const seqInfo& firstRead, const seqInfo& secondRead,
			bool local);
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
  bool checkForTandemRepeatGap();

  static bool checkTwoEqualSeqs(const std::string& seq1,
                                const std::string& seq2,
                                int allowableMismatches);
  // finding tandem repeats
  std::vector<tandemRepeat> findTandemRepeatsInSequence(
      const std::string& str, int match = 2, int mismatch = -2, int gap = -7,
      int minimumAlignScore = 50);
  tandemRepeat findTandemRepeatOfStrInSequence(std::string str,
                                               std::string tandem,
                                               int match = 2, int mismatch = -2,
                                               int gap = -7,
                                               int minimumAlignScore = 50);
  tandemRepeat findTandemRepeatOfStrInSequenceDegen(
      std::string str, std::string tandem, int match = 2, int mismatch = -2,
      int gap = -7, int minimumAlignScore = 50);
  static bool checkTwoStringsDegen(
      const std::string& str1, const std::string& str2, int allowableMismatches,
      const substituteMatrix& scoringArray);
  size_t getAlignPosForSeqAPos(size_t seqAPos);
  size_t getAlignPosForSeqBPos(size_t seqBPos);
  size_t getSeqPosForAlnAPos(size_t alnAPos);
  size_t getSeqPosForAlnBPos(size_t alnBPos);

 public:
  void setGeneralScorring(int32_t generalMatch, int32_t generalMismatch);

};



}  // namespace bibseq


