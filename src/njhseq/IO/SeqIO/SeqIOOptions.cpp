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
#include "SeqIOOptions.hpp"

#include <njhcpp/debug.h>


namespace njhseq {

bool SeqIOOptions::isInGz() const{
	return inFormat_ == inFormats::FASTAGZ
			|| inFormat_ == inFormats::FASTQGZ
			|| inFormat_ == inFormats::FASTQPAIREDGZ;
}
bool SeqIOOptions::isOutGz() const{
	return outFormat_ == outFormats::FASTQPAIRED
			|| outFormat_ == outFormats::FASTQPAIREDGZ
			|| outFormat_ == outFormats::FASTQPAIREDGZ;
}

bool SeqIOOptions::isPairedIn() const {
	return inFormat_ == inFormats::FASTQPAIRED
			|| inFormat_ == inFormats::FASTQPAIREDGZ;
}

bool SeqIOOptions::isPairedOut() const {
	return outFormat_ == outFormats::FASTQPAIRED
			|| outFormat_ == outFormats::FASTQPAIREDGZ;
}

bool SeqIOOptions::inExists() const{
	return njh::files::bfs::exists(firstName_);
}

bool SeqIOOptions::outExists() const{
	return out_.outExists();
}

SeqIOOptions::inFormats SeqIOOptions::getInFormat(const std::string & format){
	inFormats in = inFormats::NOFORMAT;
	if (format == "fastq" || format == "fq" || format == "fnq") {
		in =  inFormats::FASTQ;
	} else if (format == "fastqPaired" || "_R1.fastq"  == format || "_1.fastq" == format) {
		in =  inFormats::FASTQPAIRED;
	} else if (format == "fastqPairedgz" || "_R1.fastq.gz"  == format || "_1.fastq.gz" == format) {
		in =  inFormats::FASTQPAIREDGZ;
	} else if (format == "fastqgz" || ".fastq.gz" == format || "fastq.gz" == format) {
		in =  inFormats::FASTQGZ;
	} else if (format == "fasta" || format == "fa") {
		in =  inFormats::FASTA;
	} else if (format == "fastagz" || ".fasta.gz" == format || "fasta.gz" == format) {
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
		ss << "Acceptable types are fastq,fq,fnq,fastqPaired,_R1.fastq,_1.fastq,fastqPairedgz,_R1.fastq.gz,_1.fastq.gz,fastqgz,fastq.gz,.fastq.gz,.fasta,fa,fasta,fa,fastagz,.fasta.gz,bam,fastaQual,sff, and sffbin" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	return in;
}

SeqIOOptions::outFormats SeqIOOptions::getOutFormat(const std::string & format){
	outFormats out = outFormats::NOFORMAT;
	if (format == "fasta" || format == "fa") {
		out = outFormats::FASTA;
	} else if (format == "fastagz") {
		out = outFormats::FASTAGZ;
	} else if (format == "fastq") {
		out = outFormats::FASTQ;
	} else if (format == "fastqgz") {
		out = outFormats::FASTQGZ;
	} else if (format == "fastaQual") {
		out = outFormats::FASTAQUAL;
	} else if (format == "fastqPaired") {
		out = outFormats::FASTQPAIRED;
	} else if (format == "fastqPairedgz") {
		out = outFormats::FASTQPAIREDGZ;
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

SeqIOOptions::inFormats SeqIOOptions::getInFormatFromFnp(const bfs::path & fnp){
	auto fnpStr = fnp.string();
	inFormats in = inFormats::NOFORMAT;
	if (njh::endsWith(fnpStr, ".fastq") || njh::endsWith(fnpStr, ".fq")
			|| njh::endsWith(fnpStr, ".fnq")) {
		in =  inFormats::FASTQ;
	} else if (njh::endsWith(fnpStr, "_R1.fastq") || njh::endsWith(fnpStr, "_1.fastq")) {
		in =  inFormats::FASTQPAIRED;
	} else if (njh::endsWith(fnpStr, "_R1.fastq.gz") || njh::endsWith(fnpStr, "_1.fastq.gz")) {
		in =  inFormats::FASTQPAIREDGZ;
	} else if (njh::endsWith(fnpStr, ".fastq.gz")) {
		in =  inFormats::FASTQGZ;
	} else if (njh::endsWith(fnpStr, ".fasta") || njh::endsWith(fnpStr, ".fa")) {
		in =  inFormats::FASTA;
	} else if (njh::endsWith(fnpStr, ".fasta.gz")) {
		in =  inFormats::FASTAGZ;
	} else if (njh::endsWith(fnpStr, ".bam")) {
		in =  inFormats::BAM;
	} else if (njh::endsWith(fnpStr, ".sff.txt")) {
		in =  inFormats::SFFTXT;
	} else if (njh::endsWith(fnpStr, ".sff")){
		in =  inFormats::SFFBIN;
	}else {
		std::stringstream ss;
		ss << "Unrecognized file extension for : " << fnp
				<< ", in " << __PRETTY_FUNCTION__ << std::endl;
		throw std::runtime_error { ss.str() };
	}
	return in;
}

SeqIOOptions::outFormats SeqIOOptions::getOutFormatFromFnp(const bfs::path & fnp){
	auto fnpStr = fnp.string();
	outFormats out = outFormats::NOFORMAT;
	if (njh::endsWith(fnpStr, ".fastq") || njh::endsWith(fnpStr, ".fq") || njh::endsWith(fnpStr, ".fnq")) {
		out =  outFormats::FASTQ;
	} else if (njh::endsWith(fnpStr, "_R1.fastq") || njh::endsWith(fnpStr, "_1.fastq")) {
		out =  outFormats::FASTQPAIRED;
	} else if (njh::endsWith(fnpStr, "_R1.fastq.gz") || njh::endsWith(fnpStr, "_1.fastq.gz")) {
		out =  outFormats::FASTQPAIREDGZ;
	} else if (njh::endsWith(fnpStr, ".fastq.gz")) {
		out =  outFormats::FASTQGZ;
	} else if (njh::endsWith(fnpStr, ".fasta") || njh::endsWith(fnpStr, ".fa")) {
		out =  outFormats::FASTA;
	} else if (njh::endsWith(fnpStr, ".fasta.gz")) {
		out =  outFormats::FASTAGZ;
	} else if (njh::endsWith(fnpStr, ".bam")) {
		out =  outFormats::FASTQ;
	} else if (njh::endsWith(fnpStr, ".sff.txt")) {
		out =  outFormats::FASTQ;
	} else if (njh::endsWith(fnpStr, ".sff")){
		out =  outFormats::FASTQ;
	}else {
		std::stringstream ss;
		ss << "Unrecognized file extension for : " << fnp
				<< ", in " << __PRETTY_FUNCTION__ << std::endl;
		throw std::runtime_error { ss.str() };
	}
	return out;
}


std::string SeqIOOptions::getOutExtension() const{
	return getOutExtension(outFormat_);
}



bfs::path SeqIOOptions::getPriamryOutName() const{
	return njh::appendAsNeededRet(out_.outFilename_.string(), getOutExtension(outFormat_));
}

bfs::path SeqIOOptions::getSecondaryOutName() const{
	return njh::appendAsNeededRet(out_.outFilename_.string(), getOutExtensionSecondary(outFormat_));
}


std::string SeqIOOptions::getOutExtensionSecondary(outFormats format){
	std::string out = "noformat";
	std::stringstream ss;
	switch (format) {
	case outFormats::FLOWMOTHUR:
	case outFormats::FLOW:
	case outFormats::FASTAGZ:
	case outFormats::FASTQGZ:
	case outFormats::FASTQ:
	case outFormats::FASTA:
		ss << __PRETTY_FUNCTION__ << ": error, no secondary out for " << getOutFormat(format)
						<< std::endl;
				throw std::runtime_error { ss.str() };
		break;
	case outFormats::FASTAQUAL:
		out = ".qual";
		break;
	case outFormats::FASTQPAIRED:
		out = "_R2.fastq";
		break;
	case outFormats::FASTQPAIREDGZ:
		out = "_R2.fastq.gz";
		break;
	case outFormats::NOFORMAT:
		ss << "Format set NORMAT in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
		break;
	default:
		ss << "Unrecognized out file type or " << ", in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
		break;
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
	case outFormats::FASTQGZ:
		out = ".fastq.gz";
		break;
	case outFormats::FASTA:
		out = ".fasta";
		break;
	case outFormats::FASTAGZ:
		out = ".fasta.gz";
		break;
	case outFormats::FASTAQUAL:
		out = ".fasta";
		break;
	case outFormats::FASTQPAIRED:
		out = "_R1.fastq";
		break;
	case outFormats::FASTQPAIREDGZ:
		out = "_R1.fastq.gz";
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
	case inFormats::FASTQPAIREDGZ:
		in = "fastqPairedgz";
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
	case inFormats::FASTQPAIREDGZ:
		out = outFormats::FASTQPAIREDGZ;
		break;
	case inFormats::FASTQGZ:
		out = outFormats::FASTQGZ;
		break;
	case inFormats::FASTAGZ:
		out = outFormats::FASTAGZ;
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

SeqIOOptions::inFormats SeqIOOptions::getInFormat(outFormats format){
	SeqIOOptions::inFormats out = SeqIOOptions::inFormats::NOFORMAT;
	std::stringstream ss;
	switch (format) {
	case outFormats::FASTQ:
		out = SeqIOOptions::inFormats::FASTQ;
		break;
	case outFormats::FASTA:
		out = SeqIOOptions::inFormats::FASTA;
		break;
	case outFormats::FASTQGZ:
		out = SeqIOOptions::inFormats::FASTQGZ;
		break;
	case outFormats::FASTAGZ:
		out = SeqIOOptions::inFormats::FASTAGZ;
		break;
	case outFormats::FASTAQUAL:
		out = SeqIOOptions::inFormats::FASTAQUAL;
		break;
	case outFormats::FASTQPAIRED:
		out = SeqIOOptions::inFormats::FASTQPAIRED;
		break;
	case outFormats::FASTQPAIREDGZ:
		out = SeqIOOptions::inFormats::FASTQPAIREDGZ;
		break;
	case outFormats::FLOW:
		out = SeqIOOptions::inFormats::SFFBIN;
		break;
	case outFormats::FLOWMOTHUR:
		out = SeqIOOptions::inFormats::SFFBIN;
		break;
	case outFormats::NOFORMAT:
		out = SeqIOOptions::inFormats::NOFORMAT;
		break;
	default:
		ss << "Unrecognized out file type : " << ", in " << __PRETTY_FUNCTION__
				<< std::endl;
		throw std::runtime_error { ss.str() };
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
	case outFormats::FASTQGZ:
		out = "fastqgz";
		break;
	case outFormats::FASTAGZ:
		out = "fastagz";
		break;
	case outFormats::FASTAQUAL:
		out = "fastaQual";
		break;
	case outFormats::FASTQPAIRED:
		out = "fastqPaired";
		break;
	case outFormats::FASTQPAIREDGZ:
		out = "fastqPairedgz";
		break;
	case outFormats::FLOW:
		out = "flow";
		break;
	case outFormats::FLOWMOTHUR:
		out = "flowMothur";
		break;
	case outFormats::NOFORMAT:
		ss << "Format set NOFORMAT in " << __PRETTY_FUNCTION__
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

SeqIOOptions SeqIOOptions::genFastqIn(const bfs::path & inFilename, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = inFilename;
	ret.inFormat_ = inFormats::FASTQ;
	ret.outFormat_ = outFormats::FASTQ;
	ret.out_.outExtention_ = ".fastq";
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastqInOut(const bfs::path & inFilename,
		const bfs::path & outFilename, bool processed){
	SeqIOOptions ret = genFastqIn(inFilename);
	ret.out_.outFilename_ =  outFilename;
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastqOut(const bfs::path & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTQ;
	ret.outFormat_ = outFormats::FASTQ;
	ret.out_.outExtention_ = ".fastq";
	return ret;
}


SeqIOOptions SeqIOOptions::genFastqInGz(const bfs::path & inFilename, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = inFilename;
	ret.inFormat_ = inFormats::FASTQGZ;
	ret.outFormat_ = outFormats::FASTQGZ;
	ret.out_.outExtention_ = ".fastq.gz";
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastqOutGz(const bfs::path & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTQGZ;
	ret.outFormat_ = outFormats::FASTQGZ;
	ret.out_.outExtention_ = ".fastq.gz";
	return ret;
}

SeqIOOptions SeqIOOptions::genFastqInOutGz(const bfs::path & inFilename,
		const bfs::path & outFilename, bool processed){
	SeqIOOptions ret = genFastqInGz(inFilename, processed);
	ret.out_.outFilename_ =  outFilename;
	return ret;
}


SeqIOOptions SeqIOOptions::genPairedOut(const bfs::path & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTQPAIRED;
	ret.outFormat_ = outFormats::FASTQPAIRED;
	ret.out_.outExtention_ = "_R1.fastq";
	return ret;
}

SeqIOOptions SeqIOOptions::genPairedIn(const bfs::path & r1reads, const bfs::path & r2reads, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = r1reads;
	ret.secondName_ = r2reads;
	ret.inFormat_ = inFormats::FASTQPAIRED;
	ret.outFormat_ = outFormats::FASTQPAIRED;
	ret.out_.outExtention_ = "_R1.fastq";
	ret.processed_ = processed;
	return ret;
}
SeqIOOptions SeqIOOptions::genPairedInOut(const bfs::path & r1reads, const bfs::path & r2reads,const bfs::path & outFilename, bool processed){
	SeqIOOptions ret = genPairedIn(r1reads, r2reads, processed );
	ret.out_.outFilename_ =  outFilename;
	return ret;
}


SeqIOOptions SeqIOOptions::genPairedOutGz(const bfs::path & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTQPAIREDGZ;
	ret.outFormat_ = outFormats::FASTQPAIREDGZ;
	ret.out_.outExtention_ = "_R1.fastq.gz";
	return ret;
}

SeqIOOptions SeqIOOptions::genPairedInGz(const bfs::path & r1reads, const bfs::path & r2reads, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = r1reads;
	ret.secondName_ = r2reads;
	ret.inFormat_ = inFormats::FASTQPAIREDGZ;
	ret.outFormat_ = outFormats::FASTQPAIREDGZ;
	ret.out_.outExtention_ = "_R1.fastq.gz";
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genPairedInOutGz(const bfs::path & r1reads, const bfs::path & r2reads,const bfs::path & outFilename, bool processed){
	SeqIOOptions ret = genPairedInGz(r1reads, r2reads, processed );
	ret.out_.outFilename_ =  outFilename;
	return ret;
}



SeqIOOptions SeqIOOptions::genFastaIn(const bfs::path & inFilename, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = inFilename;
	ret.inFormat_ = inFormats::FASTA;
	ret.outFormat_ = outFormats::FASTA;
	ret.out_.outExtention_ = ".fasta";
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaOut(const bfs::path & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTA;
	ret.outFormat_ = outFormats::FASTA;
	ret.out_.outExtention_ = ".fasta";
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaInOut(const bfs::path & inFilename,
		const bfs::path & outFilename, bool processed){
	SeqIOOptions ret = genFastaIn(inFilename);
	ret.out_.outFilename_ = outFilename;
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaInGz(const bfs::path & inFilename, bool processed){
	SeqIOOptions ret;
	ret.firstName_ = inFilename;
	ret.inFormat_ = inFormats::FASTAGZ;
	ret.outFormat_ = outFormats::FASTAGZ;
	ret.out_.outExtention_ = ".fasta.gz";
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaOutGz(const bfs::path & outFilename){
	SeqIOOptions ret;
	ret.out_.outFilename_ = outFilename;
	ret.inFormat_ = inFormats::FASTAGZ;
	ret.outFormat_ = outFormats::FASTAGZ;
	ret.out_.outExtention_ = ".fasta.gz";
	return ret;
}

SeqIOOptions SeqIOOptions::genFastaInOutGz(const bfs::path & inFilename,
		const bfs::path & outFilename, bool processed){
	SeqIOOptions ret = genFastaInGz(inFilename);
	ret.out_.outFilename_ = outFilename;
	ret.processed_ = processed;
	return ret;
}

SeqIOOptions::SeqIOOptions() {

}

SeqIOOptions::SeqIOOptions(const std::string & jsonStr) {
	auto root = njh::json::parse(jsonStr);
	firstName_ = root.get("firstName_", "").asString();
	secondName_ = root.get("secondName_", "").asString();
	inFormat_ = getInFormat(root.get("inFormat_", "").asString());
	outFormat_ = getOutFormat(root.get("outFormat_", "").asString());
	out_ = OutOptions(root.get("out_",""));
	//revComplMate_ = root.get("complementMate_", false).asBool();
	processed_ = root.get("processed_", false).asBool();
	lowerCaseBases_ = root.get("lowerCaseBases_", "").asString();
	removeGaps_ = root.get("removeGaps_", false).asBool();
	includeWhiteSpaceInName_ = root.get("includeWhiteSpaceInName_", true).asBool();
	extra_ = root.get("extra_", true).asInt();
}

Json::Value SeqIOOptions::toJson() const {
	Json::Value ret;
	ret["firstName_"] = njh::json::toJson(firstName_);
	ret["secondName_"] = njh::json::toJson(secondName_);
	ret["inFormat_"] = njh::json::toJson(inFormat_ == inFormats::NOFORMAT ? std::string("NOFORMAT") : getInFormat(inFormat_));
	ret["outFormat_"] = njh::json::toJson(outFormat_ == outFormats::NOFORMAT ? std::string("NOFORMAT") :getOutFormat(outFormat_));
	ret["out_"] = out_.toJson();
	ret["revComplMate_"] = njh::json::toJson(revComplMate_);
	ret["processed_"] = njh::json::toJson(processed_);
	ret["lowerCaseBases_"] = njh::json::toJson(lowerCaseBases_);
	ret["removeGaps_"] = njh::json::toJson(removeGaps_);
	ret["includeWhiteSpaceInName_"] = njh::json::toJson(includeWhiteSpaceInName_);
	ret["extra_"] = njh::json::toJson(extra_);
	return ret;
}

SeqIOOptions::SeqIOOptions(const OutOptions & out, outFormats outFormat) :
		outFormat_(outFormat), out_(out) {
	inFormat_ = getInFormat(outFormat);
}

SeqIOOptions::SeqIOOptions(const bfs::path & outFilename, outFormats outFormat) :
		outFormat_(outFormat), out_(outFilename) {
	inFormat_ = getInFormat(outFormat);
	if("" == out_.outExtention_&&
			outFormat != SeqIOOptions::outFormats::NOFORMAT){
		out_.outExtention_ = getOutExtension(outFormat_);
	}else if("" != out_.outExtention_ && !njh::beginsWith(out_.outExtention_,".f") &&
			(outFormat == SeqIOOptions::outFormats::FASTA || outFormat == SeqIOOptions::outFormats::FASTQ)){
		out_.outExtention_ = getOutExtension(outFormat_);
	}
}

SeqIOOptions::SeqIOOptions(const bfs::path & outFilename, outFormats outFormat,
		const OutOptions & out) :
		outFormat_(outFormat), out_(out) {
	out_.outFilename_ = outFilename;
	if("" == out_.outExtention_&&
			outFormat != SeqIOOptions::outFormats::NOFORMAT){
		out_.outExtention_ = getOutExtension(outFormat_);
	}else if("" != out_.outExtention_ && !njh::beginsWith(out_.outExtention_,".f") &&
			(outFormat == SeqIOOptions::outFormats::FASTA || outFormat == SeqIOOptions::outFormats::FASTQ)){
		out_.outExtention_ = getOutExtension(outFormat_);
	}

	inFormat_ = getInFormat(outFormat);
}

SeqIOOptions::SeqIOOptions(const bfs::path & firstName, inFormats inFormat,
		bool processed) :
		firstName_(firstName), inFormat_(inFormat), outFormat_(
				SeqIOOptions::getOutFormat(inFormat)), processed_(processed) {


}


}  // namespace njhseq
