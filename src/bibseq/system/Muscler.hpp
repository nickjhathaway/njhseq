#pragma once
/*
 * Muscler.hpp
 *
 *  Created on: Jan 30, 2017
 *      Author: nick
 */

#include <bibcpp/system.h>
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/alignment/alnCache/alnInfoGlobal.hpp"
#include "bibseq/IO/SeqIO/SeqInput.hpp"
#include "bibseq/alignment/aligner/alignCalc.hpp"


namespace bibseq {


class Muscler{

	alnInfoGlobal genGlobalAlnInfo(const std::string & seq);

	bfs::path musclePath_;
public:



	Muscler(const bfs::path & musclePath);

	Muscler();

	void setMusclePath(const bfs::path & musclePath);

	/**@brief run muscle on this file
	 *
	 * @param filename name of the fasta file
	 * @return a vector of the algined seqs
	 */
	std::vector<seqInfo> muscleSeqs(const SeqIOOptions & opts);

	/**@brief muscle the sequences in seqs, leave the original seqs alone
	 *
	 * @param seqs the sequence to align
	 * @return a vector of the alinged seqs
	 */
	template<typename T>
	std::vector<T> muscleSeqsRet(std::vector<T> seqs){
		muscleSeqs(seqs);
		return seqs;
	}

	/**@brief muscle the sequences in seqs
	 *
	 * @param seqs the sequence to align
	 * @param selected only muscle these selected sequences at these positions, will throw if out of range
	 */
	template<typename T, typename POSTYPE>
	void muscleSeqs(std::vector<T> & seqs,
			const std::vector<POSTYPE> & selected){
		//create temporary file, the last 6 xs will be randomized characters
		char *tmpname = strdup("/tmp/tmpfileXXXXXX");
		mkstemp(tmpname);
		uint32_t seqsWritten = 0;
		{
			//in it's own scope so that that tFile gets flushed at termination
			std::ofstream tFile(tmpname);
			if(!tFile){
				throw std::runtime_error{bib::bashCT::boldRed("Error in opening " + std::string(tmpname))};
			}
			//make name the read position as muscle will reorganize the seqs afterwards

			for (const auto & pos : selected) {
				if (pos >= seqs.size()) {
					throw std::out_of_range {
							"Error in bibseq::sys::muscleSeqs, position out of range, pos: "
									+ estd::to_string(pos) + ", size: "
									+ estd::to_string(seqs.size()) };
				}
				//hack because muscle doesn't like stop codons
				if(!bib::containsSubString(getSeqBase(seqs[pos]).seq_, "*")){
					++seqsWritten;
					tFile << ">" << pos << "\n";
					tFile << getSeqBase(seqs[pos]).seq_ << "\n";
				}
			}
		}
		//this is for when there are no sequences written due to all having stop codons which muscle won't align;
		if(seqsWritten > 0){
			std::vector<readObject> tempObjs;
			try {
				std::vector<std::string> cmds { musclePath_.string(), "-quiet", "-in", tmpname };
				auto rOut = bib::sys::run(cmds);
				if(!rOut.success_){
					std::stringstream sErr;
					sErr << bib::bashCT::red << "failure:" << std::endl;
					sErr << rOut.stdErr_ << std::endl;
					sErr << bib::bashCT::reset << std::endl;
					throw std::runtime_error{sErr.str()};
				}
				std::stringstream ss(rOut.stdOut_);
				SeqIOOptions opts;
				SeqInput reader(opts);
				seqInfo seq;
				while(reader.readNextFastaStream(ss, seq,false)){
					auto & currentRead = getSeqBase(seqs[std::stoul(seq.name_)]);
					auto gAlnInfo = genGlobalAlnInfo(seq.seq_);
					alignCalc::rearrangeGlobalQueryOnly(currentRead.seq_, '-', gAlnInfo );
					alignCalc::rearrangeGlobalQueryOnly(currentRead.qual_, 0, gAlnInfo );
				}
			} catch (std::exception & e) {
				bib::files::bfs::remove(tmpname);
				throw e;
			}
			bib::files::bfs::remove(tmpname);
		}
	}


	struct MusPosSize {

		MusPosSize();
		MusPosSize(size_t pos);
		MusPosSize(size_t pos, uint32_t size);

		size_t pos_;
		uint32_t size_ = std::numeric_limits<uint32_t>::max();

	};

	template<typename T>
	void muscleSeqs(std::vector<T> & seqs,
			const std::unordered_map<uint32_t, MusPosSize> & posSizes){
		//create temporary file, the last 6 xs will be randomized characters
		char *tmpname = strdup("/tmp/tmpfileXXXXXX");
		mkstemp(tmpname);
		uint32_t seqsWritten = 0;
		std::unordered_map<uint32_t, std::shared_ptr<seqInfo>> subInfos;
		{
			//in it's own scope so that that tFile gets flushed at termination
			std::ofstream tFile(tmpname);
			if(!tFile){
				throw std::runtime_error{bib::bashCT::boldRed("Error in opening " + std::string(tmpname))};
			}
			//make name the read position as muscle will reorganize the seqs afterwards

			for (const auto & pos : posSizes) {
				if (pos.first >= seqs.size()) {
					throw std::out_of_range {
							"Error in " + std::string(__PRETTY_FUNCTION__) + ", position out of range, pos: "
									+ estd::to_string(pos.first) + ", size: "
									+ estd::to_string(seqs.size()) };
				}
				if(pos.second.pos_ > len(getSeqBase(seqs[pos.first]))){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error pos.second.pos_ " << pos.second.pos_ << " is greater than the length of "
							<< getSeqBase(seqs[pos.first]).name_ << ", " << len(getSeqBase(seqs[pos.first])) << "\n";
					throw std::out_of_range {ss.str()};
				}
				std::shared_ptr<seqInfo> subSeq;
				if(std::numeric_limits<uint32_t>::max() == pos.second.size_ ){
					subSeq = std::make_shared<seqInfo>(getSeqBase(seqs[pos.first]).getSubRead(pos.second.pos_));
				}else{
					subSeq = std::make_shared<seqInfo>(getSeqBase(seqs[pos.first]).getSubRead(pos.second.pos_, pos.second.size_));
				}
				//hack because muscle doesn't like stop codons
				if(!bib::containsSubString(subSeq->seq_, "*")){
					++seqsWritten;
					tFile << ">" << pos.first << "\n";
					tFile << subSeq->seq_ << "\n";
					subInfos[pos.first] = subSeq;
				}
			}
		}
		//this is for when there are no sequences written due to all having stop codons which muscle won't align;
		if(seqsWritten > 0){
			std::vector<readObject> tempObjs;
			try {
				std::vector<std::string> cmds { musclePath_.string(), "-quiet", "-in", tmpname };
				auto rOut = bib::sys::run(cmds);
				if(!rOut.success_){
					std::stringstream sErr;
					sErr << bib::bashCT::red << "failure:" << std::endl;
					sErr << rOut.stdErr_ << std::endl;
					sErr << bib::bashCT::reset << std::endl;
					throw std::runtime_error{sErr.str()};
				}
				std::stringstream ss(rOut.stdOut_);
				SeqIOOptions opts;
				SeqInput reader(opts);
				seqInfo seq;
				while(reader.readNextFastaStream(ss, seq,false)){
					uint32_t pos = bib::lexical_cast<uint32_t>(seq.name_);
					auto gAlnInfo = genGlobalAlnInfo(seq.seq_);
					alignCalc::rearrangeGlobalQueryOnly(subInfos[pos]->seq_, '-', gAlnInfo );
					alignCalc::rearrangeGlobalQueryOnly(subInfos[pos]->qual_, 0, gAlnInfo );
					if(std::numeric_limits<uint32_t>::max() == posSizes.at(pos).size_){
						getSeqBase(seqs[pos]).trimBack(posSizes.at(pos).pos_);
					}else{
						getSeqBase(seqs[pos]).clipOut(posSizes.at(pos).pos_, posSizes.at(pos).size_);
					}
					getSeqBase(seqs[pos]).insert(posSizes.at(pos).pos_, *(subInfos[pos]));
				}
			} catch (std::exception & e) {
				bib::files::bfs::remove(tmpname);
				throw e;
			}
			bib::files::bfs::remove(tmpname);
		}
	}

	/**@brief muscle the sequences in seqs
	 *
	 * @param seqs the sequence to align
	 */
	template<typename T>
	void muscleSeqs(std::vector<T> & seqs){
		std::vector<uint64_t> allSelected(seqs.size());
		bib::iota<uint64_t>(allSelected, 0);
		muscleSeqs(seqs, allSelected);
	}

};

} /* namespace bibseq */

