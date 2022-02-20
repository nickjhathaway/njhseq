#pragma once
/*
 * HmmerTableDomainHit.hpp
 *
 *  Created on: Jan 09, 2022
 *      Author: nick
 */



#include "njhseq/common.h"
#include "njhseq/utils.h"


#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"

#include "njhseq/objects/BioDataObject/HmmerObjs/HmmerHitBase.hpp"

namespace njhseq {


/**
 * @class HmmerSeqHitTab
 * @brief class to hold the per line info for hits output by the DNA search functions of hmmer3
 *
 */
class HmmerSeqHitTab : public HmmerHitBase{

public:
	HmmerSeqHitTab();
	HmmerSeqHitTab(const std::string & line);


//	(1) target name: The name of the target sequence or profile.
//	(2) accession: The accession of the target sequence or profile, or ’-’ if none.
//	(3) query name: The name of the query sequence or profile.
//	(4) accession: The accession of the query sequence or profile, or ’-’ if none.
//	(5) hmmfrom: The position in the hmm at which the hit starts.
//	(6) hmm to: The position in the hmm at which the hit ends.
//	(7) alifrom: The position in the target sequence at which the hit starts.
//	(8) ali to: The position in the target sequence at which the hit ends.
//	(9) envfrom: The position in the target sequence at which the surrounding envelope starts.
//	(10) env to: The position in the target sequence at which the surrounding envelope ends.
//	(11) modlen: The length of the model.
//	(12) strand: The strand on which the hit was found (“-" when alifrom>ali to).
//	(13) E-value: The expectation value (statistical significance) of the target, as above.
//	(14) score (full sequence): The score (in bits) for this hit. It includes the biased-composition correction.
//	(15) Bias (full sequence): The biased-composition correction, as above
//	(16) description of target: The remainder of the line is the target’s description line, as free text.



	std::string targetName_{""};
	std::string targetAcc_{""};

	std::string queryName_{""};
	std::string queryAcc_{""};
//	uint32_t hmmFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
//	uint32_t hmmTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec
//
//	uint32_t alignFrom_ {std::numeric_limits<uint32_t>::max()};//1-based as per file spec
//	uint32_t alignTo_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
//
//	uint32_t envFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
//	uint32_t envTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	uint32_t modelLen_{std::numeric_limits<uint32_t>::max()};

	char strand_ {' '};

	double modelEvalue_ {std::numeric_limits<double>::max()};
	double modelScore_ {std::numeric_limits<double>::max()};
	double modelBias_ {std::numeric_limits<double>::max()};

	std::string targetDesc_;

//	bool isReverseStrand() const;

	virtual Json::Value toJson() const;
	virtual std::string toDelimStr() const;

	static VecStr toDelimStrHeader ();

	virtual ~HmmerSeqHitTab(){

	}
	double modelCoverage() const;
	double queryCoverageEnv(uint32_t queryLen) const;
	double queryCoverageAln(uint32_t queryLen) const;


//	uint32_t env0BasedPlusStrandStart() const;
//	uint32_t env0BasedPlusStrandEnd() const;
//	uint32_t align0BasedPlusStrandStart() const;
//	uint32_t align0BasedPlusStrandEnd() const;
//
//	uint32_t zeroBasedHmmFrom() const; //per file specs the positions are 1-based
//
//	uint32_t envLen() const;
//	uint32_t hmmLen() const;
//	uint32_t alignLen() const;



	Bed6RecordCore genBed6_env() const;
	Bed6RecordCore genBed6_aln() const;
private:
	Bed6RecordCore genBed6(uint32_t start, uint32_t end) const;
public:
};


}  // namespace njhseq



