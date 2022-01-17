#pragma once
/*
 * HmmerDomainHitTab.hpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */



#include "njhseq/common.h"
#include "njhseq/utils.h"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"

namespace njhseq {


/**
 * @class HmmerDomainHitTab
 * @brief class to hold the the per line info for the per domain hit for the protein functions of hmmer3
 *
 */
class HmmerDomainHitTab{
	/*
	(1) target name: The name of the target sequence or profile.
	(2) target accession: Accession of the target sequence or profile, or ’-’ if none is available.
	(3) tlen: Length of the target sequence or profile, in residues. This (together with the query length) is
	useful for interpreting where the domain coordinates (in subsequent columns) lie in the sequence.
	(4) query name: Name of the query sequence or profile.
	(5) accession: Accession of the target sequence or profile, or ’-’ if none is available.
	(6) qlen: Length of the query sequence or profile, in residues.
	(7) E-value: E-value of the overall sequence/profile comparison (including all domains).
	(8) score: Bit score of the overall sequence/profile comparison (including all domains), inclusive of a
	null2 bias composition correction to the score.
	(9) bias: The biased composition score correction that was applied to the bit score.
	(10) #: This domain’s number (1..ndom).
	(11) of: The total number of domains reported in the sequence, ndom.
	47
	(12) c-Evalue: The “conditional E-value”, a permissive measure of how reliable this particular domain
	may be. The conditional E-value is calculated on a smaller search space than the independent Evalue. The conditional E-value uses the number of targets that pass the reporting thresholds. The
	null hypothesis test posed by the conditional E-value is as follows. Suppose that we believe that
	there is already sufficient evidence (from other domains) to identify the set of reported sequences as
	homologs of our query; now, how many additional domains would we expect to find with at least this
	particular domain’s bit score, if the rest of those reported sequences were random nonhomologous
	sequence (i.e. outside the other domain(s) that were sufficient to identified them as homologs in the
	first place)?
	(13) i-Evalue: The “independent E-value”, the E-value that the sequence/profile comparison would have
	received if this were the only domain envelope found in it, excluding any others. This is a stringent
	measure of how reliable this particular domain may be. The independent E-value uses the total
	number of targets in the target database.
	(14) score: The bit score for this domain.
	(15) bias: The biased composition (null2) score correction that was applied to the domain bit score.
	(16) from (hmm coord): The start of the MEA alignment of this domain with respect to the profile, numbered 1..N for a profile of N consensus positions.
	(17) to (hmm coord): The end of the MEA alignment of this domain with respect to the profile, numbered
	1..N for a profile of N consensus positions.
	(18) from (ali coord): The start of the MEA alignment of this domain with respect to the sequence,
	numbered 1..L for a sequence of L residues.
	(19) to (ali coord): The end of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
	(20) from (env coord): The start of the domain envelope on the sequence, numbered 1..L for a sequence of L residues. The envelope defines a subsequence for which their is substantial probability
	mass supporting a homologous domain, whether or not a single discrete alignment can be identified.
	The envelope may extend beyond the endpoints of the MEA alignment, and in fact often does, for
	weakly scoring domains.
	(21) to (env coord): The end of the domain envelope on the sequence, numbered 1..L for a sequence
	of L residues.
	(22) acc: The mean posterior probability of aligned residues in the MEA alignment; a measure of how
	reliable the overall alignment is (from 0 to 1, with 1.00 indicating a completely reliable alignment
	according to the model).
	(23) description of target: The remainder of the line is the target’s description line, as free text.*/


public:

	HmmerDomainHitTab();
	HmmerDomainHitTab(const std::string & line);

	std::string targetName_{""};
	std::string targetAcc_{""};
	uint32_t targetLen_ {std::numeric_limits<uint32_t>::max()};

	std::string queryName_{""};
	std::string queryAcc_{""};
	uint32_t queryLen_ {std::numeric_limits<uint32_t>::max()};

	double seqEvalue_ {std::numeric_limits<double>::max()};
	double seqScore_ {std::numeric_limits<double>::max()};
	double seqBias_ {std::numeric_limits<double>::max()};

	uint32_t domainId_ {std::numeric_limits<uint32_t>::max()};
	uint32_t domainsTotals_ {std::numeric_limits<uint32_t>::max()};
	double domain_c_evalue_ {std::numeric_limits<double>::max()};
	double domain_i_evalue_ {std::numeric_limits<double>::max()};
	double domainScore_ {std::numeric_limits<double>::max()};
	double domainBias_ {std::numeric_limits<double>::max()};

	uint32_t hmmFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t hmmTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	uint32_t alignFrom_ {std::numeric_limits<uint32_t>::max()};//1-based as per file spec
	uint32_t alignTo_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec

	uint32_t envFrom_ {std::numeric_limits<uint32_t>::max()};  //1-based as per file spec
	uint32_t envTo_ {std::numeric_limits<uint32_t>::max()};    //1-based as per file spec

	double acc_{std::numeric_limits<double>::max()};
	std::string targetDesc_{""};

	uint32_t zeroBasedHmmFrom() const; //per file specs the positions are 1-based
	uint32_t zeroBasedAlignFrom() const;//per file specs the positions are 1-based
	uint32_t zeroBasedEnvFrom() const;//per file specs the positions are 1-based

	uint32_t envLen() const;
	uint32_t hmmLen() const;
	uint32_t alignLen() const;

	double modelCoverage() const;

	Bed6RecordCore genBed6_env() const;
	Bed6RecordCore genBed6_aln() const;
private:
	Bed6RecordCore genBed6(uint32_t start, uint32_t end) const;
public:


	Json::Value toJson() const;
	std::string toDelimStr() const;
	static VecStr toDelimStrHeader ();

};


}  // namespace njhseq



