#pragma once
/*
 * MidDeterminator.hpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */

#include "bibseq/objects/helperObjects/motif.hpp"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include "bibseq/objects/dataContainers/table.hpp"


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
	MidDeterminator(const table & mids);

	MidDeterminator(){};

	std::unordered_map<std::string, motif> mids_;
	std::unordered_map<std::string, motif> compMids_;

	bool containsMidByName(const std::string & name)const;
	bool containsMidByBarcode(const std::string & barcode, bool checkComp)const;

	std::string getMidName(const std::string & barcode, bool checkComp)const;

	void addBarcode(const std::string & name, const std::string & barcode,
			bool addToComp) ;




	midPos determineMidSimple(const std::string & seq) ;

	midPos determineMidSimpleBack(const std::string & seq) ;

	midPos determineMidSimpleComp(const std::string & seq) ;

	midPos determineMidSimpleCompFront(const std::string & seq);

	midPos determineMidPosVarStart(const std::string & seq, uint32_t varStop);

	midPos determineMidPosVarStartBack(const std::string & seq, uint32_t varStop);

	midPos determineMidPosVarStartComp(const std::string & seq, uint32_t varStop) ;

	midPos determineMidPosVarStartCompFront(const std::string & seq, uint32_t varStop);

	midPos fullDetermine(seqInfo & info, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds );

	template<typename T>
	midPos fullDetermine(T & read, bool variableStart, uint32_t variableStop, bool checkComplement, bool barcodesBothEnds ){
		return fullDetermine(read.seqBase_, variableStart, variableStop, checkComplement, barcodesBothEnds);
	}
};

} /* namespace bibseq */


