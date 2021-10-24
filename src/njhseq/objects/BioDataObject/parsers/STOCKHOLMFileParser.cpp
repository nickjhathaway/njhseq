/*
 * STOCKHOLMFileParser.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */

#include "STOCKHOLMFileParser.hpp"


namespace njhseq {

STOCKHOLMFileParser::STOCKHOLMFileParser(const bfs::path & fnp):fnp_(fnp){
	if("STDIN" != fnp.string()){
		njh::files::checkExistenceThrow(fnp, __PRETTY_FUNCTION__);
	}
}

void STOCKHOLMFileParser::parseFile(){
	perFileAnnotation_.clear();
	perColumnAnnotation_.clear();
	seqs_.clear();
	perResidueAnnotation_.clear();
	InputStream in{InOptions(fnp_)};
	std::string line  = "";
	//get header
	uint32_t lineNumber = 1;
	njh::files::crossPlatGetline(in, line);
	if("# STOCKHOLM 1.0" != line ){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", first line in " << fnp_ << " was not \"# STOCKHOLM 1.0\" reading not supported" << "\n";
		throw std::runtime_error{ss.str()};
	}

	while(njh::files::crossPlatGetline(in, line)){
		++lineNumber;
		//end of reading
		if(njh::beginsWith(line, "//")){
			break;
		}
		//skip over empty lines
		if("" == line || allWhiteSpaceStr(line)){
			continue;
		}
		if(line.size() >4 && '#' == line.front()){
			std::string front = line.substr(0, 4);
			if("#=GF" == front){
				auto toks = tokenizeString(line, "whitespace");
				if(toks.size() < 3){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", lines should have at least 3 items separated by whitespace" << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto feature = toks[1];
				if(njh::in(feature, perFileAnnotation_)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", already have feature: " << feature << "\n";
					throw std::runtime_error{ss.str()};
				}
				perFileAnnotation_[feature]= njh::conToStr(getSubVector(toks, 2), " ");
			}else if("#=GC" == front){
				auto toks = tokenizeString(line, "whitespace");
				if(toks.size() != 3){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", lines should have 3 items separated by whitespace" << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(std::numeric_limits<uint32_t>::max() == numColumns_){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", " << front << " lines should come after alignment" << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(toks[2].size() != numColumns_){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line "
							<< lineNumber << ", for " << front << ", " << front
							<< " 3rd item should be the same size as the number of columns, #of cols: "
							<< numColumns_ << ", item size: " << toks[2].size() << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto feature = toks[1];
				if(njh::in(feature, perColumnAnnotation_)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", already have column feature: " << feature << "\n";
					throw std::runtime_error{ss.str()};
				}
				perColumnAnnotation_[feature] = toks[2];
			}else if("#=GS" == front){
				auto toks = tokenizeString(line, "whitespace");
				if(toks.size() < 4){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", lines should have at least 4 items separated by whitespace" << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto seqname = toks[1];
				auto feature = toks[2];
				if(njh::in(feature, perSeqAnnotation_[seqname])){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", already have feature: " << feature << " for seq " << seqname << "\n";
					throw std::runtime_error{ss.str()};
				}
				perSeqAnnotation_[seqname][feature]= njh::conToStr(getSubVector(toks, 3), " ");
			}else if("#=GR" == front){
				auto toks = tokenizeString(line, "whitespace");
				if(toks.size() != 4){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", lines should have 4 items separated by whitespace" << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto seqname = toks[1];
				auto feature = toks[2];
				if(!njh::in(seqname, readInSeqs_)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << " have not read in seq: " << seqname  << ", yet"<< "\n";
					throw std::runtime_error{ss.str()};
				}
				if(njh::in(feature, perResidueAnnotation_[seqname])){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", for " << front << ", already have feature: " << feature << " for seq " << seqname << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(toks[3].size() != numColumns_){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line "
							<< lineNumber << ", for " << front << ", " << front
							<< " 4th item should be the same size as the number of columns, #of cols: "
							<< numColumns_ << ", item size: " << toks[3].size() << "\n";
					throw std::runtime_error{ss.str()};
				}
				perResidueAnnotation_[seqname][feature] = toks[3];
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", started with " << front << " which is not a recognized pattern" << "\n";
				throw std::runtime_error{ss.str()};
			}
		}else{
			auto toks = tokenizeString(line, "whitespace");
			if(toks.size() != 2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", lines should have 2 items separated by whitespace" << "\n";
				throw std::runtime_error{ss.str()};
			}
			if(njh::in(toks[0], readInSeqs_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error in processing line " << lineNumber << ", already read in seq " << toks[0] << "\n";
				throw std::runtime_error{ss.str()};
			}
			seqs_.emplace_back(toks[0], toks[1]);
			//first seq
			if (readInSeqs_.empty()) {
				numColumns_ = len(seqs_.back());
			} else {
				if (numColumns_ != len(seqs_.back())) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line "
							<< lineNumber
							<< ", sequence is the same number of residudes as the first the seq, size: "
							<< len(seqs_.back()) << ", num of cols: " << numColumns_
							<< "\n";
					throw std::runtime_error { ss.str() };
				}
			}
			readInSeqs_.emplace_back(seqs_.back().name_);
		}
	}
}

void STOCKHOLMFileParser::writeOutSeqFile(const SeqIOOptions & opts) const{
	SeqOutput::write(seqs_, opts);
}

void STOCKHOLMFileParser::writeOutFile(const OutOptions & opts) const{
	OutputStream out(opts);
	out << "# STOCKHOLM 1.0" << std::endl;
	for(const auto & gf : perFileAnnotation_){
		out << "#=GF " << gf.first << " " << gf.second << std::endl;
	}
	out << "" << std::endl;
	uint32_t maxnameSize = 0;
	for (const auto & name : readInSeqs_) {
		if (maxnameSize < name.size()) {
			maxnameSize = name.size();
		}
	}
	for(const auto & seq : seqs_){
		if(njh::in(seq.name_, perSeqAnnotation_)){
			for (const auto & feature : perSeqAnnotation_.at(seq.name_)) {
				out << "#=GS " << seq.name_
						<< std::string(maxnameSize - seq.name_.size(), ' ') << " "
						<< feature.first << " " << feature.second << std::endl;
			}
		}
	}

	out << "" << std::endl;
	for(const auto & seq : seqs_){
		out << seq.name_ << std::string(maxnameSize + 5 - seq.name_.size(), ' ')
		<< " " << "  " << " "
		<< seq.seq_ << std::endl;
		if(njh::in(seq.name_, perResidueAnnotation_)){
			for(const auto & feature : perResidueAnnotation_.at(seq.name_)){
				out << "#=GR "  << seq.name_ << std::string(maxnameSize - seq.name_.size(), ' ')
				 << " " << feature.first << " " << feature.second << std::endl;
			}
		}
	}
	for(const auto & gc : perColumnAnnotation_){
		out << "#=GC " << gc.first << std::string(maxnameSize - gc.first.size(), ' ')
		<< " " << "  " << " " << gc.second << std::endl;
	}
	out << "//" << std::endl;;
}

}  // namespace njhseq

