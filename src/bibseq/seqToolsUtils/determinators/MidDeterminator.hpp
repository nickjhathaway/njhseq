#pragma once
/*
 * MidDeterminator.hpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */
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

	struct MidDeterminePars{
		uint32_t variableStop_ = 0;
		bool checkComplement_ = false;
		bool barcodesBothEnds_ = false;
		bool checkForShorten_ = false;
	};

	/**@b create MidDterminator from a table with ids and barcodes
	 *
	 * @param mids a table that has at least two columns titled id(string name of mid) and barcode(nucleoties)
	 * @todo check for whether there are several mids of the same barcode by accident
	 */
	explicit MidDeterminator(const table & mids);

	MidDeterminator();

	std::unordered_map<std::string, MidInfo> mids_;
private:
	uint32_t allowableMismatches_ = 0;
	bool midEndsRevComp_ = false;
public:

	void setAllowableMismatches(uint32_t allowableMismatches);
	void setMidEndsRevComp(bool midEndsRevComp);

	bool containsMidByName(const std::string & name) const;
	bool containsMidByBarcode(const std::string & barcode) const;

	std::string getMidName(const std::string & barcode) const;



	void addBarcode(const std::string & name, const std::string & barcode);



	std::vector<midPos> determinePossibleMidPos(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosBack(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosComp(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosCompBack(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosShorten(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosBackShorten(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosCompShorten(const std::string & seq, uint32_t within);

	std::vector<midPos> determinePossibleMidPosCompBackShorten(const std::string & seq, uint32_t within);

	void processInfoWithMidPos(seqInfo & info,
			const MidDeterminator::midPos & pos,
			MidDeterminator::MidDeterminePars pars);

	void processInfoWithMidPos(seqInfo & info,
			const MidDeterminator::midPos & frontPos,
			const MidDeterminator::midPos & backPos);

	void processInfoWithMidPos(PairedRead & info,
			const MidDeterminator::midPos & pos,
			MidDeterminator::MidDeterminePars pars);

	void processInfoWithMidPos(PairedRead & info,
			const MidDeterminator::midPos & frontPos,
			const MidDeterminator::midPos & backPos);

	std::pair<midPos, midPos> fullDetermine(seqInfo & info, MidDeterminePars pars);
	std::pair<midPos, midPos> fullDetermine(PairedRead & info, MidDeterminePars pars);

	template<typename T>
	std::pair<midPos, midPos> fullDetermine(T & read, MidDeterminePars pars){
		return fullDetermine(getSeqBase(read), pars);
	}

	static void increaseFailedBarcodeCounts(const MidDeterminator::midPos & pos, std::unordered_map<std::string, uint32_t> & counts);

};

} /* namespace bibseq */


