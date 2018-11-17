#pragma once
/*
 * STOCKHOLMFileParser.hpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */


#include "njhseq/IO.h"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

namespace njhseq {

class STOCKHOLMFileParser{
	//parser made based off of https://en.wikipedia.org/wiki/Stockholm_format

	/**
	 * 	#=GF <feature> <Generic per-File annotation, free text>
			#=GC <feature> <Generic per-Column annotation, exactly 1 char per column>
			#=GS <seqname> <feature> <Generic per-Sequence annotation, free text>
			#=GR <seqname> <feature> <Generic per-Residue annotation, exactly 1 char per residue>
	 */
public:

	STOCKHOLMFileParser(const bfs::path & fnp);
	const bfs::path fnp_;

	std::unordered_map<std::string, std::string> perFileAnnotation_;
	std::unordered_map<std::string, std::string> perColumnAnnotation_;

	std::vector<seqInfo> seqs_;
	VecStr readInSeqs_;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> perResidueAnnotation_;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> perSeqAnnotation_;

	uint32_t numColumns_{std::numeric_limits<uint32_t>::max()};


	void parseFile();
	void writeOutSeqFile(const SeqIOOptions & opts) const;
	void writeOutFile(const OutOptions & opts) const;

};



}  // namespace njhseq



