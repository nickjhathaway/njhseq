#pragma once
/*
 * MidDeterminator.hpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */
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
#include "bibseq/objects/helperObjects/motif.hpp"
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/objects/seqObjects/Paired/PairedRead.hpp"

#include "bibseq/objects/dataContainers/tables/table.hpp"


namespace bibseq {





class MidDeterminator {
public:
	struct midPos {
		midPos();
		midPos(const std::string & midName, uint64_t midPos, uint64_t barcodeSize, uint32_t barcodeScore);
		std::string midName_;
		uint64_t midPos_;
		uint64_t barcodeSize_;
		uint32_t barcodeScore_;
		bool inRevComp_ = false;
		bool failure_ = false;
		std::string altName_ = "";

		enum class FailureCase{
			NONE,
			NO_MATCHING,
			PARTIAL,
			TOO_MANY_MATCHING,
			MISMATCHING_MIDS,
			MISMATCHING_DIRECTION
		};

		FailureCase fCase_ = FailureCase::NONE;

		static std::string getFailureCaseName(FailureCase fCase);
		static VecStr getFailureCaseNames();

		double normalizeScoreByLen() const;

		explicit operator bool() const;

		Json::Value toJson() const;
		static const uint64_t npos = std::numeric_limits<uint64_t>::max();
	};

	struct MidInfo{

		MidInfo(const std::string & midName, const std::string & barcode);

		std::string midName_; /**< the name of the barcode*/

		std::unique_ptr<motif> bar_; /**< the barcode */
		std::unique_ptr<motif> rcompBar_; /**< reverse complement of the barcode */

		std::unique_ptr<motif> shortenFrontBar_; /**< the barcode with the first base missing*/
		std::unique_ptr<motif> shortenFrontRCompBar_; /**< the reverse complement of the barcode with the first base missing */

		std::unique_ptr<motif> shortenBackBar_; /**< the barcode with the last base missing*/
		std::unique_ptr<motif> shortenBackRCompBar_; /**< the reverse complement of the barcode with the last base missing */

		bool rCompSame_;
	};

	struct MID{
		MID(const std::string & name);

		MID(const std::string & name, const std::string & forwardBar);

		MID(const std::string & name,
				const std::string & forwardBar,
				const std::string & reverseBar);

		std::string name_;

		std::unique_ptr<MidInfo> forwardBar_;
		std::unique_ptr<MidInfo> reverseBar_;

		bool forSameAsRev_ { false };
		bool forSameAsRevShorten_ { false };

		bool dualBarcoded() const;

	};

	struct MidDeterminePars{
		uint32_t searchStart_ = 0;
		uint32_t searchStop_ = 0;
		uint32_t allowableErrors_ = 0;
		bool checkComplement_ = false;
		bool checkForShorten_ = false;
	};


	struct MidSearchRes{
		std::vector<MidDeterminator::midPos> forward_;
		std::vector<MidDeterminator::midPos> reverse_;
	};

	struct ProcessedRes{
		enum class PROCESSED_CASE{
			NONE,
			NO_MATCHING,
			PARTIALDUAL,
			TOO_MANY_MATCHING,
			MISMATCHING_MIDS,
			MISMATCHING_DIRECTION,
			MATCH
		};

		static std::string getProcessedCaseName(PROCESSED_CASE pCase);
		static VecStr getProcessedCaseNames();

		bool rcomplement_{false};
		std::string midName_;
		PROCESSED_CASE case_{PROCESSED_CASE::NONE};


	};


  MidDeterminator(const table & mids, const MidDeterminePars & searchPars);
	MidDeterminator(const bfs::path & idFileFnp, const MidDeterminePars & searchPars);

	std::unordered_map<std::string, MID> mids_;
	MidDeterminePars searchPars_;
	MidDeterminePars shortenSearchPars_;

private:
public:

	void containsMidByNameThrow(const std::string & name, const std::string & funcName) const;
	bool containsMidByName(const std::string & name) const;
	bool containsMidByBarcode(const std::string & barcode) const;

	std::string getMidName(const std::string & barcode) const;

	void addForwardReverseBarcode(const std::string & name, const std::string & forward, const std::string & reverse);
	void addForwardBarcode(const std::string & name, const std::string & forward);
	void addReverseBarcode(const std::string & name, const std::string & reverse);



	static std::vector<midPos> frontDeterminePosMIDPos (
			const std::string & seq,
			const motif & bar,
			const std::string & midName,
			const MidDeterminePars & mPars);

	static std::vector<midPos> backDeterminePosMIDPos (
			const std::string & seq,
			const motif & bar,
			const std::string & midName,
			const MidDeterminePars & mPars);



	MidSearchRes searchPairedEndRead(const PairedRead & seq) const;
	ProcessedRes processSearchPairedEndRead(PairedRead & seq, const MidSearchRes & res) const;

	MidSearchRes searchRead(const seqInfo & seq) const;
	ProcessedRes processSearchRead(seqInfo & seq, const MidSearchRes & res) const;

};




} /* namespace bibseq */


