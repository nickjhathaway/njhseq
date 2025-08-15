//
// Created by Nicholas Hathaway on 8/8/25.
//


#include "CodonSampler.hpp"

namespace njhseq {
namespace aminoAcidInfo {

CodonSampler::CodonSampler(): amino_acid_info_(infos::allInfo) {
  std::unordered_map<std::string, double> codon_usage;
  for (const auto & amino_acid : amino_acid_info_) {
    for (const auto & dna_codon : amino_acid.second.dnaCodons_) {
      codon_usage_[dna_codon] = 1.0;
    }
  }
  set_generators();
}

CodonSampler::CodonSampler(const std::unordered_map<std::string, double> &codon_usage): codon_usage_(codon_usage),
amino_acid_info_(infos::allInfo) {
  set_generators();
}

CodonSampler::CodonSampler(const std::unordered_map<std::string, double> &codon_usage,
const std::unordered_map<char, aminoAcid> & amino_acid_infos): codon_usage_(codon_usage),
amino_acid_info_(amino_acid_infos) {
  set_generators();
}


void CodonSampler::set_generators(
  const std::unordered_map<std::string, double> &codon_usage,
  const std::unordered_map<char, aminoAcid> & amino_acid_info) {
  codon_generators_.clear();
  //first check that all codons in usage are 3 letter codons
  VecStr warnings;
  for (const auto & codon : codon_usage) {
    if (codon.first.size() != 3) {
      warnings.emplace_back(njh::pasteAsStr("codons must be length 3, not ", codon.first.size(), " for ", codon.first));
    }
  }
  //check that all codons in amino_acid_info can be found within codon_usage
  for (const auto & amino_acid : amino_acid_info) {
    for (const auto & dna_codon : amino_acid.second.dnaCodons_) {
      if (njh::notIn(dna_codon, codon_usage)) {
        warnings.emplace_back(njh::pasteAsStr("missing codon ", dna_codon, " for ", amino_acid.second.letCode_));
      }
    }
  }
  if (!warnings.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", errors " << " " << "\n";
    ss << njh::conToStr(warnings, "\n") << "\n";
    throw std::runtime_error{ss.str()};
  }
  for (const auto & amino_acid : amino_acid_info) {
    std::vector<double> dna_codon_weights;
    dna_codon_weights.reserve(amino_acid.second.dnaCodons_.size());
    for (const auto & dna_codon : amino_acid.second.dnaCodons_) {
      dna_codon_weights.emplace_back(codon_usage_[dna_codon]);
    }
    codon_generators_.emplace(amino_acid.second.letCode_, std::make_shared<njh::randObjectGen<std::string, double>>(amino_acid.second.dnaCodons_, dna_codon_weights));
  }
  //create a reverse look up table
  for (const auto & amino_acid : amino_acid_info_) {
    for (const auto & dna_codon : amino_acid.second.dnaCodons_) {
      reverse_dna_codon_lookup_[dna_codon] = amino_acid.first;
    }
  }
}

void CodonSampler::set_generators() {
  set_generators(codon_usage_, amino_acid_info_);
}


void CodonSampler::set_seed(uint64_t seed) {
  for (auto & gen : codon_generators_) {
    gen.second->set_seed(seed);
  }
}


void CodonSampler::set_stop_codon_usage_to_amber_only(bool reset_generators) {
  codon_usage_[OPAL_STOP_DNA_CODON] = 0;
  codon_usage_[OCHRE_STOP_DNA_CODON] = 0;

  if (reset_generators) {
    set_generators();
  }
}


std::string CodonSampler::gen(char aa) {
  return codon_generators_[aa]->genObj();
}

std::string CodonSampler::gen(const std::string & protein) {
  std::string ret;
  ret.reserve(protein.size() * 3);
  for (const auto aa : protein) {
    if (njh::notIn(aa, codon_generators_)) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " aa: " << aa << " not found in codon_generators_" << "\n";
      ss << "option are: " << njh::conToStr(njh::getSetOfMapKeys(codon_generators_), ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
    ret.append(gen(aa));
  }
  return ret;
}

std::string CodonSampler::gen_no_check(const std::string & protein) {
  std::string ret;
  ret.reserve(protein.size() * 3);
  for (const auto aa : protein) {
    ret.append(gen(aa));
  }
  return ret;
}

VecStr CodonSampler::get_alt_codons(const std::string & codon) const {
  VecStr ret;
  if (3 != codon.size()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " codon must be size 3 not " << codon.size() << " for " << codon << "\n";
    throw std::runtime_error{ss.str()};
  }

  if (njh::notIn(codon, reverse_dna_codon_lookup_)) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " can't find " << codon << " in reverse look up table" << "\n";
    ss << "options are: " << njh::conToStr(njh::getSetOfMapKeys(reverse_dna_codon_lookup_), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  ret = amino_acid_info_.at(reverse_dna_codon_lookup_.at(codon)).dnaCodons_;
  //remove current codon
  removeElement(ret, codon);
  //sort by usage
  njh::sort(ret, [this](const std::string & codon1, const std::string & codon2) {
    return codon_usage_.at(codon1) < codon_usage_.at(codon2);
  });
  return ret;
}


} //namespace aminoAcidInfo
} //namespace njhseq



