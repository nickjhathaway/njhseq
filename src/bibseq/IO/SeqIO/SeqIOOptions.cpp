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
#include "SeqIOOptions.hpp"

#include <bibcpp/debug.h>

namespace bibseq {


bool SeqIOOptions::inExists() const{
	return bib::files::bfs::exists(firstName_);
}

bool SeqIOOptions::outExists() const{
	return bib::files::bfs::exists(bib::appendAsNeededRet(out_.outFilename_, out_.outExtention_));
}

SeqIOOptions::inFormats SeqIOOptions::getInFormat(const std::string & format){
	inFormats in = inFormats::NOFORMAT;
	if (format == "fastq" || format == "fq"
			|| format == "fnq") {
		in =  inFormats::FASTQ;
	} else if (format == "fastqPaired") {
		in =  inFormats::FASTQPAIRED;
	}  else if (format == "fastqgz") {
		in =  inFormats::FASTQGZ;
	} else if (format == "fasta") {
		in =  inFormats::FASTA;
	} else if (format == "fastagz") {
		in =  inFormats::FASTAGZ;
	} else if (format == "bam") {
		in =  inFormats::BAM;
	} else if (format == "fastaQual") {
		in =  inFormats::FASTAQUAL;
	} else if (format == "sff") {
		in =  inFormats::SFFTXT;
	} else if (format == "sffbin"){
		in =  inFormats::SFFBIN;
	}else {
		std::stringstream ss;
		ss << "Unrecognized file type : " << format
				<< ", in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "Acceptable types are fasta,fastaQual,fastq,fastqPaired, bam, fastagz, sff, and sffbin" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	return in;
}

SeqIOOptions::outFormats SeqIOOptions::getOutFormat(const std::string & format){
	outFormats out = outFormats::NOFORMAT;
	if (format == "fasta") {
		out = outFormats::FASTA;
	} else if (format == "fastq") {
		out = outFormats::FASTQ;
	} else if (format == "fastaQual") {
		out = outFormats::FASTAQUAL;
	} else if (format == "fastqPaired") {
		out = outFormats::FASTQPAIRED;
	} else if (format == "flow") {
		out = outFormats::FLOW;
	} else if (format == "flowMothur") {
		out = outFormats::FLOWMOTHUR;
	} else {
		std::stringstream ss;
		ss << "Unrecognized out file type : " << format << ", in " << __PRETTY_FUNCTION__
				<< std::endl;
		ss << "Acceptable types are fasta, fastq, fastqPaired, and fastaQual" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	return out;
}

std::string SeqIOOptions::getOutExtension(outFormats format){
	std::string out = "noformat";
	std::stringstream ss;
	switch (format) {
	case outFormats::FASTQ:
		out = ".fastq";
		break;
	case outFormats::FASTA:
		out = ".fasta";
		break;
	case outFormats::FASTAQUAL:
		out = ".fasta";
		break;
	case outFormats::FASTQPAIRED:
		out = "_R1.fastq";
		break;
	case outFormats::FLOW:
		out = ".dat";
		break;
	case outFormats::FLOWMOTHUR:
		out = ".flow";
		break;
	case outFormats::NOFORMAT:
		ss << "Format set NORMAT in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
		break;
	default:
		ss << "Unrecognized out file type : " << ", in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
		break;
	}
	return out;
}

std::string SeqIOOptions::getInFormat(inFormats format){
	std::string in = "noformat";
	std::stringstream ss;
	switch (format) {
	case inFormats::FASTQ:
		in = "fastq";
		break;
	case inFormats::FASTA:
		in = "fasta";
		break;
	case inFormats::BAM:
		in = "bam";
		break;
	case inFormats::FASTQPAIRED:
		in = "fastqPaired";
		break;
	case inFormats::FASTQGZ:
		in = "fastqgz";
		break;
	case inFormats::FASTAGZ:
		in = "fastagz";
		break;
	case inFormats::FASTAQUAL:
		in = "FASTAQUAL";
		break;
	case inFormats::SFFTXT:
		in = "sff";
		break;
	case inFormats::SFFBIN:
		in = "sffbin";
		break;
	case inFormats::NOFORMAT:
		ss << "Format set to NOFORMAT in " << __PRETTY_FUNCTION__ << std::endl;
		throw std::runtime_error { ss.str() };
		break;
	default:
		ss << "Unrecognized file type in " << __PRETTY_FUNCTION__ << std::endl;
		throw std::runtime_error { ss.str() };
		break;
	}
	return in;
}

SeqIOOptions::outFormats SeqIOOptions::getOutFormat(inFormats format){
	SeqIOOptions::outFormats out = outFormats::NOFORMAT;
	switch (format) {
	case inFormats::FASTQ:
		out = outFormats::FASTQ;
		break;
	case inFormats::FASTA:
		out = outFormats::FASTA;
		break;
	case inFormats::BAM:
		out = outFormats::FASTQ;
		break;
	case inFormats::FASTQPAIRED:
		out = outFormats::FASTQPAIRED;
		break;
	case inFormats::FASTQGZ:
		out = outFormats::FASTQ;
		break;
	case inFormats::FASTAGZ:
		out = outFormats::FASTA;
		break;
	case inFormats::FASTAQUAL:
		out = outFormats::FASTAQUAL;
		break;
	case inFormats::SFFTXT:
		out = outFormats::FASTQ;
		break;
	case inFormats::SFFBIN:
		out = outFormats::FASTQ;
		break;
	default:
		out = outFormats::NOFORMAT;
		break;
	}
	return out;
}
std::string SeqIOOptions::getOutFormat(outFormats format){
	std::string out = "noformat";
	std::stringstream ss;
	switch (format) {
	case outFormats::FASTQ:
		out = "fastq";
		break;
	case outFormats::FASTA:
		out = "fasta";
		break;
	case outFormats::FASTAQUAL:
		out = "fastaQual";
		break;
	case outFormats::FASTQPAIRED:
		out = "fastqPaired";
		break;
	case outFormats::FLOW:
		out = "flow";
		break;
	case outFormats::FLOWMOTHUR:
		out = "flowMothur";
		break;
	case outFormats::NOFORMAT:
		ss << "Format set NORMAT in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
		break;
	default:
		ss << "Unrecognized out file type : " << ", in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
		break;
	}
	return out;
}

SeqIOOptions SeqIOOptions::genFastqIn(const std::string & inFilename, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = inFilename;
	ret.inFormat_ = inFormats::FASTQ;
	ret.outFormat_ = outFormats::FASTQ;
	ret.out_.outExtention_ = ".fastq";
	ret.processed_ = processed;
	return ret;
}
SeqIOOptions SeqIOOptions::genFastqInOut(const std::string & inFilename,
		const std::string & outFilename, bool processed){
	SeqIOOptions ret = genFastqIn(inFilename);
	ret.out_.outFilename_ =  outFilename;
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastqOut(const std::string & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTQ;
	ret.outFormat_ = outFormats::FASTQ;
	ret.out_.outExtention_ = ".fastq";
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaIn(const std::string & inFilename, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = inFilename;
	ret.inFormat_ = inFormats::FASTA;
	ret.outFormat_ = outFormats::FASTA;
	ret.out_.outExtention_ = ".fasta";
	ret.processed_ = processed;
	return ret;
}
SeqIOOptions SeqIOOptions::genFastaInOut(const std::string & inFilename,
		const std::string & outFilename, bool processed){
	SeqIOOptions ret = genFastaIn(inFilename);
	ret.out_.outFilename_ = outFilename;
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaOut(const std::string & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTA;
	ret.outFormat_ = outFormats::FASTA;
	ret.out_.outExtention_ = ".fasta";
	return ret;
}

SeqIOOptions::SeqIOOptions() {

}

SeqIOOptions::SeqIOOptions(const std::string & jsonStr) {
	Json::Reader reader;
	Json::Value root;
	auto stats =  reader.parse(jsonStr,root);
	if(!stats){
		std::cerr << "Error in parsing jsonStr for readObjectIOOptions in " << __PRETTY_FUNCTION__ << std::endl;
		std::cerr << jsonStr << std::endl;
		throw bib::err::Exception(bib::err::F() << "Could not construct from jsonStr");
	}
	firstName_ = root.get("firstName_", "").asString();
	secondName_ = root.get("secondName_", "").asString();
	inFormat_ = getInFormat(root.get("inFormat_", "").asString());
	outFormat_ = getOutFormat(root.get("outFormat_", "").asString());
	out_ = OutOptions(root.get("out_",""));
	complementMate_ = root.get("complementMate_", false).asBool();
	processed_ = root.get("processed_", false).asBool();
	lowerCaseBases_ = root.get("lowerCaseBases_", "").asString();
	removeGaps_ = root.get("removeGaps_", false).asBool();
	includeWhiteSpaceInName_ = root.get("includeWhiteSpaceInName_", true).asBool();
	extra_ = root.get("extra_", true).asInt();
}

Json::Value SeqIOOptions::toJson() const {
	Json::Value ret;
	ret["firstName_"] = bib::json::toJson(firstName_);
	ret["secondName_"] = bib::json::toJson(secondName_);
	ret["inFormat_"] = bib::json::toJson(getInFormat(inFormat_));
	ret["outFormat_"] = bib::json::toJson(getOutFormat(outFormat_));
	ret["out_"] = out_.toJson();
	ret["complementMate_"] = bib::json::toJson(complementMate_);
	ret["processed_"] = bib::json::toJson(processed_);
	ret["lowerCaseBases_"] = bib::json::toJson(lowerCaseBases_);
	ret["removeGaps_"] = bib::json::toJson(removeGaps_);
	ret["includeWhiteSpaceInName_"] = bib::json::toJson(includeWhiteSpaceInName_);
	ret["extra_"] = bib::json::toJson(extra_);
	return ret;
}

SeqIOOptions::SeqIOOptions(const OutOptions & out, outFormats outFormat): outFormat_(outFormat), out_(out) {

}

SeqIOOptions::SeqIOOptions(const std::string & outFilename,
		outFormats outFormat, const OutOptions & out):outFormat_(outFormat),out_(out) {
	out_.outFilename_ = outFilename;
}

SeqIOOptions::SeqIOOptions(const std::string & firstName,
		 inFormats inFormat, bool processed) :
		firstName_(firstName), inFormat_(inFormat),outFormat_(SeqIOOptions::getOutFormat(inFormat)), processed_(processed) {

}


}  // namespace bibseq