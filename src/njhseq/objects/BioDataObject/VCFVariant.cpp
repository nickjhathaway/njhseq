//
// Created by Nicholas Hathaway on 5/21/23.
//


#include "VCFVariant.hpp"


namespace njhseq {

std::vector<VCFVariant> VCFVariant::readVCFLine(const std::string line){
	std::vector<VCFVariant> ret;
	auto toks = tokenizeString(line, "\t");
	if(toks.size() < 8){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "line needs to be at least 8 elements" << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::string chrom = toks[0];
	uint32_t location = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);

	std::string ref = toks[3];
	std::string seq = toks[4];
	std::string name = toks[2];


	VecStr infoToks = tokenizeString(toks[7], ";");
	if(std::string::npos == seq.find(",")){
		double freq = 1;
		//single variant
		for(const auto & info : infoToks){
			auto subToks = tokenizeString(info, "=");
			if(2 != subToks.size()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "info toks should be size 2, not " << subToks.size() << " for " << info << "\n";
				throw std::runtime_error{ss.str()};
			}
			if(subToks.front() == "AF"){
				freq = njh::StrToNumConverter::stoToNum<double>(subToks.back());
			}
		}
		ret.emplace_back(VCFVariant(GenomicRegion(name, chrom, location - 1, location - 1 + seq.size(), false	),
																ref, seq, freq
		));
	}else{
		//multiple variants at this location

		auto seqToks = tokenizeString(seq, ",");
		std::vector<double> freqs;
		for(const auto & info : infoToks){
			auto subToks = tokenizeString(info, "=");
			if(2 != subToks.size()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "info toks should be size 2, not " << subToks.size() << " for " << info << "\n";
				throw std::runtime_error{ss.str()};
			}
			if(subToks.front() == "AF"){
				auto freqToks = tokenizeString(subToks.back(), ",");
				if(freqToks.size() != seqToks.size()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "frequency toks, " << freqToks.size() << " do not equal seq toks, " << seqToks.size() << " for " << line << "\n";
					throw std::runtime_error{ss.str()};
				}
				for(const auto & freq : freqToks){
					freqs.emplace_back(njh::StrToNumConverter::stoToNum<double>(freq));
				}
			}
		}
		if(freqs.size() != seqToks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "frequencies, " << freqs.size() << " do not equal seq toks, " << seqToks.size() << " for " << line << "\n";
			throw std::runtime_error{ss.str()};
		}
		for(const auto seqTok : iter::enumerate(seqToks)){
			ret.emplace_back(VCFVariant(GenomicRegion(name, chrom, location - 1, location - 1 + seqTok.element.size(), false	),
																	ref, seqTok.element, freqs[seqTok.index]
			));
		}
	}

	return ret;
}

} //namespace njhseq



