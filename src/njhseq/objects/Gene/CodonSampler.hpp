#pragma once

//
// Created by Nicholas Hathaway on 8/8/25.
//

#include "njhseq/objects/Gene/aminoAcidInfo.hpp"
#include <njhcpp/simulation/randObjGen.hpp>

namespace njhseq {
namespace aminoAcidInfo {


class CodonSampler {
public:

	CodonSampler();
	explicit CodonSampler(const std::unordered_map<std::string, double> &codon_usage);
	CodonSampler(const std::unordered_map<std::string, double> &codon_usage,
		const std::unordered_map<char, aminoAcid> & amino_acid_info);

	std::unordered_map<char, std::shared_ptr<njh::randObjectGen<std::string, double>>> codon_generators_;//!< the generators for all codons per amino acid
	std::unordered_map<std::string, double> codon_usage_;//!< the codon usage, can be unnormalized
	std::unordered_map<char, aminoAcid> amino_acid_info_;//!< the table defining which amino acid each codon_usage goes with

	std::unordered_map<std::string, char> reverse_dna_codon_lookup_;//!< reverse codon look up table

	void set_stop_codon_usage_to_amber_only(bool reset_generators = false);

	void set_seed(uint64_t seed);
	void set_generators(const std::unordered_map<std::string, double> &codon_usage,
		const std::unordered_map<char, aminoAcid> & amino_acid_info);
	void set_generators();

	std::string gen(char aa);
	std::string gen(const std::string & protein);
	std::string gen_no_check(const std::string & protein);

	VecStr get_alt_codons(const std::string & codon) const;
};



} //namespace aminoAcidInfo
} //namespace njhseq

