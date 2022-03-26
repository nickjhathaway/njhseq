//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include <njhcpp/bashUtils.h>
#include "njhseq/IO/fileUtils.hpp"
#include "SeqIO.hpp"

namespace njhseq {

SeqIO::SeqIO(const SeqIOOptions & options) :
		ioOptions_(options),in_(options), out_(options) {
}

void SeqIO::openIn(){
	in_.openIn();
}

void SeqIO::closeIn(){
	in_.closeIn();
}

void SeqIO::openOut(){
	out_.openOut();
}

void SeqIO::closeOut(){
	out_.closeOut();
}

std::unordered_map<uint32_t, std::string> SeqIO::rewriteSeqsWithIndexAsName(const SeqIOOptions & inOpts, const SeqIOOptions & outOpts){
	std::unordered_map<uint32_t, std::string> seqKey;
	//rewrite the seqs so name is just number
	seqInfo seq;
	SeqInput reader(inOpts);
	SeqOutput writer(outOpts);
	reader.openIn();
	writer.openOut();
	uint32_t count = 0;
	while(reader.readNextRead(seq)){
		seqKey[count] = seq.name_;
		seq.name_ = njh::pasteAsStr(count);
		writer.write(seq);
		++count;
	}
	return seqKey;
}

std::unordered_map<uint32_t, std::string> SeqIO::rewriteSeqsWithIndexAsName(const SeqIOOptions & inOpts, const SeqIOOptions & outOpts, const OutOptions & seqKeyFileOutOpts) {
	auto seqKey = rewriteSeqsWithIndexAsName(inOpts, outOpts);
	OutputStream seqNameKeyOut(seqKeyFileOutOpts);
	seqNameKeyOut << "oldName\tnewName" << "\n";
	for(const auto key : iter::range(0UL, seqKey.size())){
		seqNameKeyOut << seqKey[key] << "\t" << key << "\n";
	}
	return seqKey;
}

std::unordered_map<uint32_t, std::string> SeqIO::rewriteSeqsWithIndexAsName(const std::vector<seqInfo> & seqs, const SeqIOOptions & outOpts){
	std::unordered_map<uint32_t, std::string> seqKey;
	SeqOutput writer(outOpts);
	writer.openOut();
	uint32_t count = 0;
	for(auto seq : seqs){
		seqKey[count] = seq.name_;
		seq.name_ = njh::pasteAsStr(count);
		writer.write(seq);
		++count;
	}
	return seqKey;
}

std::unordered_map<uint32_t, std::string> SeqIO::rewriteSeqsWithIndexAsName(const std::vector<seqInfo> & seqs, const SeqIOOptions & outOpts, const OutOptions & seqKeyFileOutOpts) {
	auto seqKey = rewriteSeqsWithIndexAsName(seqs, outOpts);
	OutputStream seqNameKeyOut(seqKeyFileOutOpts);
	seqNameKeyOut << "oldName\tnewName" << "\n";
	for(const auto key : iter::range(0UL, seqKey.size())){
		seqNameKeyOut << seqKey[key] << "\t" << key << "\n";
	}
	return seqKey;
}

}  // namespace njhseq
