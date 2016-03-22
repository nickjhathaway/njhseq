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
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"


namespace bibseq {


struct midPos {
	midPos() ;
	midPos(const std::string & midName, uint64_t midPos, uint64_t barcodeSize) ;
	std::string midName_;
	uint64_t midPos_;
	uint64_t barcodeSize_;

	explicit operator bool() const ;

	Json::Value toJson()const;
	static const uint64_t npos = std::numeric_limits<uint64_t>::max();
};

class MidDeterminator {
public:
	/**@b create MidDterminator from a table with ids and barcodes
	 *
	 * @param mids a table that has at least two columns titled id(string name of mid) and barcode(nucleoties)
	 * @todo check for whether there are several mids of the same barcode by accident
	 */
	explicit MidDeterminator(const table & mids);

	MidDeterminator(){};

	std::unordered_map<std::string, motif> mids_;
	std::unordered_map<std::string, motif> compMids_;

	bool containsMidByName(const std::string & name)const;
	bool containsMidByBarcode(const std::string & barcode, bool checkComp)const;

	std::string getMidName(const std::string & barcode, bool checkComp)const;

	void addBarcode(const std::string & name, const std::string & barcode,
			bool addToComp) ;




	midPos determineMidSimple(const std::string & seq);
	midPos determineMidSimple(const std::string & seq, uint32_t allowableMismatches);

	midPos determineMidSimpleBack(const std::string & seq);
	midPos determineMidSimpleBack(const std::string & seq,  uint32_t allowableMismatches);

	midPos determineMidSimpleComp(const std::string & seq);
	midPos determineMidSimpleComp(const std::string & seq,  uint32_t allowableMismatches);

	midPos determineMidSimpleCompFront(const std::string & seq);
	midPos determineMidSimpleCompFront(const std::string & seq, uint32_t allowableMismatches);

	midPos determineMidPosVarStart(const std::string & seq, uint32_t varStop);
	midPos determineMidPosVarStart(const std::string & seq, uint32_t varStop, uint32_t allowableMismatches);

	midPos determineMidPosVarStartBack(const std::string & seq, uint32_t varStop);
	midPos determineMidPosVarStartBack(const std::string & seq, uint32_t varStop, uint32_t allowableMismatches);

	midPos determineMidPosVarStartComp(const std::string & seq, uint32_t varStop);
	midPos determineMidPosVarStartComp(const std::string & seq, uint32_t varStop, uint32_t allowableMismatches);

	midPos determineMidPosVarStartCompFront(const std::string & seq, uint32_t varStop);
	midPos determineMidPosVarStartCompFront(const std::string & seq, uint32_t varStop, uint32_t allowableMismatches);


	midPos fullDetermine(seqInfo & info, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds);
	midPos fullDetermine(seqInfo & info, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds, uint32_t allowableMismatches);

	template<typename T>
	midPos fullDetermine(T & read, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds){
		return fullDetermine(read.seqBase_, variableStart, variableStop, checkComplement, barcodesBothEnds);
	}
	template<typename T>
	midPos fullDetermine(T & read, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds, uint32_t allowableMismatches){
		return fullDetermine(read.seqBase_, variableStart, variableStop, checkComplement, barcodesBothEnds, allowableMismatches);
	}
};

} /* namespace bibseq */


