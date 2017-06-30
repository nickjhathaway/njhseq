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
#include "bibseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp"

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
//		auto mkTempRet = mkstemp(tmpname);
//		if(-1 == mkTempRet){
//			std::stringstream sErr;
//			sErr << __PRETTY_FUNCTION__ << ", error in creating file name from template " << tmpname << "\n";
//			throw std::runtime_error{sErr.str()};
//		}
//		close(mkTempRet);
		auto mkTempRet = mktemp(tmpname);
		if(nullptr == mkTempRet){
			std::stringstream sErr;
			sErr << __PRETTY_FUNCTION__ << ", error in creating file name from template " << tmpname << "\n";
			throw std::runtime_error{sErr.str()};
		}
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
					sErr << __PRETTY_FUNCTION__ << ", error " << "\n";
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
				std::cerr << e.what() << std::endl;
				bib::files::bfs::remove(tmpname);
				throw std::runtime_error{e.what()};
			}
		}
		bib::files::bfs::remove(tmpname);
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


	struct StartStopMALNPos{
		uint32_t start_;
		uint32_t stop_;
	};

	/**@brief Find the first and last non gap base in a seqs that have been multiple algined
	 *
	 * @param alnSeqs the vector of aligned seqs
	 * @return the start and stops
	 */
	template<typename T>
	static std::vector<StartStopMALNPos> getMAlnStartsAndStops(const std::vector<T> & alnSeqs){
		std::vector<StartStopMALNPos> ret;
		for(const auto & seq : alnSeqs){
			StartStopMALNPos pos;
			pos.start_ = getSeqBase(seq).seq_.find_first_not_of('-');
			pos.stop_ = getSeqBase(seq).seq_.find_last_not_of('-');
			ret.emplace_back(pos);
		}
		return ret;
	}

	/**@brief The pileup at a multiple alignment position
	 *
	 */
	struct AlnPosScore {

		AlnPosScore(seqInfo::size_type pos);

		seqInfo::size_type pos_;/**< The position in the multiple alignment */


		charCounter counter_ { std::vector<char> { 'A', 'C', 'G', 'T', '-' } };


		uint32_t baseCount_ = 0;/**< The number of bases here*/
		uint32_t gapCount_ = 0;/**< The number of gaps here*/

		/**@brief Set the number of bases and gaps at this position from the counter_ object
		 *
		 */
		void setCounts();

		/**@brief The number of bases and gaps here
		 *
		 * @return the count of bases and gaps,
		 */
		uint32_t getSpanningCount() const;


	};

	template<typename T>
	static void checkAlignSeqsLensThrow(const std::vector<T> & alnSeqs, const std::string & funcName){
		bool failed = false;
		std::stringstream ss;
		ss << funcName << ", errors found, all input seqs must have the same length" << "\n";
		for(const auto seqPos : iter::range(alnSeqs.size())){
			const auto & seq = alnSeqs[seqPos];
			if(len(getSeqBase(seq)) != len(getSeqBase(alnSeqs.front()))){
				failed = true;
				ss << "Pos: " << seqPos << ", Name " << getSeqBase(alnSeqs[seqPos]).name_ << "\n";
			}
		}
		if(failed){
			throw std::runtime_error{ss.str()};
		}
	}

	/**@brief Score each position in the multiple alignment, don't count gaps at the end and beginging of seqs so spanning counts can be calculated
	 *
	 * @param alnSeqs the aligned seqs
	 * @return a vector of counts per position
	 */
	template<typename T>
	static std::vector<std::shared_ptr<AlnPosScore>> getPileupCounts(
			const std::vector<T> & alnSeqs) {

		checkAlignSeqsLensThrow(alnSeqs, __PRETTY_FUNCTION__);

		auto startsAndStops = getMAlnStartsAndStops(alnSeqs);
		std::vector<std::shared_ptr<AlnPosScore>> ret;
		for (const auto basePos : iter::range(len(alnSeqs.front()))) {
			std::shared_ptr<AlnPosScore> score = std::make_shared<AlnPosScore>(
					basePos);
			for (const auto seqPos : iter::range(alnSeqs.size())) {
				const auto & seq = alnSeqs[seqPos];
				if (basePos >= startsAndStops[seqPos].start_
						&& basePos <= startsAndStops[seqPos].stop_) {
					score->counter_.increaseCountOfBase(seq.seq_[basePos]);
				}
			}
			score->setCounts();
			ret.emplace_back(score);
		}
		return ret;
	}

	/**@brief A streak of multiple alignment positions that pass a certain predicate
	 *
	 */
	struct AlnPosScoreStreak {
		AlnPosScoreStreak();
		AlnPosScoreStreak(const std::shared_ptr<AlnPosScore> & firstScore);
		uint32_t start_;
		uint32_t end_;
		uint32_t getLne() const;
		std::vector<std::shared_ptr<AlnPosScore>> scores_;

	};

	/**@brief Get streaks of positions that pass a certain predicate
	 *
	 * @param scores the scores to get streaks from
	 * @param scorePred the predicate, a function object takes a std::shared_ptr<AlnPosScore> and returns a bool
	 * @param streakLenCutOff the length a streak must reach to be included in the streaks
	 * @return a vector of streaks that pass the predicate and length cut off
	 */
	static std::vector<AlnPosScoreStreak> getAlignmentStreaksPositions(
			const std::vector<std::shared_ptr<AlnPosScore>> & scores,
			const std::function<bool(const std::shared_ptr<AlnPosScore>&)> & scorePred,
			uint32_t streakLenCutOff);

	/**@brief Get streaks of positions that pass a certain predicate
	 *
	 * @param alnSeqs the aligned seqs that scores will then be calculated from and then evauluated for streaks
	 * @param scorePred the predicate, a function object takes a std::shared_ptr<AlnPosScore> and returns a bool
	 * @param streakLenCutOff the length a streak must reach to be included in the streaks
	 * @return a vector of streaks that pass the predicate and length cut off
	 */
	template<typename T>
	static std::vector<AlnPosScoreStreak> getAlignmentStreaksPositions(
			const std::vector<T> & alnSeqs,
			const std::function<bool(const std::shared_ptr<AlnPosScore>&)> & scorePred,
			uint32_t streakLenCutOff){
		auto scores = getPileupCounts(alnSeqs);
		return getAlignmentStreaksPositions(scores, scorePred, streakLenCutOff);
	}


	/**@brief Trim multiple aligned sequences (must be same length)
	 *
	 * @param alnSeqs the alignment seqs to trim
	 * @param streaks the streaks calculated
	 */
	template<typename T>
	static void trimAlnSeqsToFirstAndLastStreak(std::vector<T> & alnSeqs,
			const std::vector<AlnPosScoreStreak> & streaks
			){
		if(alnSeqs.empty()){
			return;
		}
		if(streaks.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, streaks shouldn't be empty" << "\n";
			throw std::runtime_error{ss.str()};
		}
		checkAlignSeqsLensThrow(alnSeqs, __PRETTY_FUNCTION__);
		auto startsAndStops = getMAlnStartsAndStops(alnSeqs);

		uint32_t trimEnd = len(getSeqBase(alnSeqs.front())) - streaks.back().end_;
		uint32_t trimFront = streaks.front().start_;
		readVecTrimmer::trimEnds(alnSeqs, trimFront, trimEnd);
		for(const auto seqPos : iter::range(len(alnSeqs))){
			if(startsAndStops[seqPos].start_ >  streaks.front().start_ | startsAndStops[seqPos].stop_ + 1 < streaks.back().end_){
				alnSeqs[seqPos].on_ = false;
			}else{
				alnSeqs[seqPos].on_ = true;
			}
		}
		readVec::removeGapsFromReads(alnSeqs);
	}


	/**@brief Trim seqs to streaks calculated by a streak length cutoff and score predicate
	 *
	 * @param alnSeqs the alignment streaks
	 * @param scorePred the predicate to consider an alignment position
	 * @param streakLenCutOff the number of position in a row that must pass the predicate to be considered a streak
	 */
	template<typename T>
	static void trimAlnSeqsToFirstAndLastStreak(std::vector<T> & alnSeqs,
			const std::function<bool(const std::shared_ptr<AlnPosScore>&)> & scorePred,
						uint32_t streakLenCutOff
			){
		if(alnSeqs.empty()){
			return;
		}
		auto streaks = getAlignmentStreaksPositions(alnSeqs, scorePred, streakLenCutOff);
		trimAlnSeqsToFirstAndLastStreak(alnSeqs, streaks);
	}


	/**@brief Write out scores of multiple alignment
	 *
	 * @param out the out file
	 * @param alnSeqs the aligned seqs
	 * @param scores the scores from the alignemnt
	 */
	template<typename T>
	static void writeScores(std::ostream & out, const std::vector<T> & alnSeqs, const std::vector<std::shared_ptr<AlnPosScore>> & scores ){
		auto totalInputSeqs = len(alnSeqs);
		out << "position\tentropy\tbaseCount\tbaseFrac\tgapCount\tspanningCount\tpercentSpanning\ttotalReadCount" << std::endl;
		for(const auto scorePos : iter::range(scores.size())){
			auto & score = scores[scorePos];
			out << score->pos_
					<< "\t" << score->counter_.computeEntrophy()
					<< "\t" << score->baseCount_
					<< "\t" << static_cast<double>(score->baseCount_)/score->getSpanningCount()
					<< "\t" << score->gapCount_
					<< "\t" << score->getSpanningCount()
					<< "\t" << static_cast<double>(score->getSpanningCount())/totalInputSeqs
					<< "\t" << totalInputSeqs << std::endl;
		}
	}

	/**@brief Trim Input sequences to the approaximate end and starts of a multiple alignment to reference sequences
	 *
	 * @param inputSeqs the input sequences to trim
	 * @param refSeqs the ref sequence to trim to
	 */
	template<typename INPUTSEQ, typename REF>
	void trimSeqsToMultiAlnRef(std::vector<INPUTSEQ> & inputSeqs,
			const std::vector<REF> & refSeqs) {
		bool fail = false;
		std::stringstream ss;
		if (refSeqs.empty()) {
			fail = true;
			ss << __PRETTY_FUNCTION__ << ", error refSeqs is empty " << "\n";
		}

		if (inputSeqs.empty()) {
			fail = true;
			ss << __PRETTY_FUNCTION__ << ", error inputSeqs is emtpy " << "\n";
		}
		if (fail) {
			throw std::runtime_error { ss.str() };
		}
		std::vector<seqInfo> allSeqs;
		for(const auto & refSeq : refSeqs){
			allSeqs.emplace_back(getSeqBase(refSeq));
		}
		for(const auto & inputSeq : inputSeqs){
			allSeqs.emplace_back(getSeqBase(inputSeq));
		}
		muscleSeqs(allSeqs);

		std::vector<seqInfo> alignedRefs = std::vector<seqInfo>(allSeqs.begin(), allSeqs.begin() + refSeqs.size());


		uint32_t streakLenCutOff = 3; // at least 3 positions in a row must pass the threshold below
		double spanningCutOff = .50; // at least 50% of the ref seqs must start here
		double baseCutOff = .50; // at least 50% of the bases at this location must have a base
		uint32_t hardGapCutOff = 2; //there must at least be 2 refs here

		auto totalInputSeqs = len(alignedRefs);
		std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred = [&baseCutOff,&hardGapCutOff,&totalInputSeqs,&spanningCutOff](const std::shared_ptr<Muscler::AlnPosScore> & score){
			if((static_cast<double>(score->baseCount_)/score->getSpanningCount() >= baseCutOff | score->gapCount_ <= hardGapCutOff) &&
					static_cast<double>(score->getSpanningCount())/totalInputSeqs > spanningCutOff){
				return true;
			}else{
				return false;
			}
		};

		auto refStartsStop = getMAlnStartsAndStops(alignedRefs);
		//count up each location
		auto scores = Muscler::getPileupCounts(alignedRefs);
		auto streaks = Muscler::getAlignmentStreaksPositions(alignedRefs, scorePred, streakLenCutOff);

		if(streaks.empty()){
			std::cerr << "No streaks passed current filters, try being less stringent" << std::endl;
		}else{
			Muscler::trimAlnSeqsToFirstAndLastStreak(allSeqs, streaks);
			readVec::removeGapsFromReads(allSeqs);
			for(const auto pos : iter::range(refSeqs.size(), allSeqs.size())){
				const auto inputSeqPos =  pos - refSeqs.size();
				getSeqBase(inputSeqs[inputSeqPos] ) = allSeqs[pos];
			}
		}
	}
};

} /* namespace bibseq */

