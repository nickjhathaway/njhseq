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


namespace njhseq {


class HmmerTableDomainHit{
public:
	HmmerTableDomainHit();
	HmmerTableDomainHit(const std::string & line);


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
	uint32_t hmmFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t hmmTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	uint32_t alignFrom_ {std::numeric_limits<uint32_t>::max()};//1-based as per file spec
	uint32_t alignTo_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec

	uint32_t envFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t envTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	uint32_t modelLen_{std::numeric_limits<uint32_t>::max()};

	char strand_ {' '};

	double modelEvalue_ {std::numeric_limits<double>::max()};
	double modelScore_ {std::numeric_limits<double>::max()};
	double modelBias_ {std::numeric_limits<double>::max()};

	std::string targetDesc_;

	bool isReverseStrand() const;

	Json::Value toJson() const;
	std::string toDelimStr() const;

	static VecStr toDelimStrHeader ();


	Bed3RecordCore genBed3() const;
	Bed6RecordCore genBed6() const;

};


}  // namespace njhseq



