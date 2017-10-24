/*
 * Muscler.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: nick
 */

#include "Muscler.hpp"

namespace bibseq {



alnInfoGlobal Muscler::genGlobalAlnInfo(const std::string & seq){
	uint64_t pos = 0;
	uint32_t gapOffSet = 0;
	uint32_t gapSize = 0;
	uint32_t gapStartSiteInSeq = 0;
	//first is pos, second is gap size
	std::vector<std::pair<uint32_t, uint32_t>> gapInfos;
	while(pos != seq.size()){
		if(seq[pos] == '-'){
			if(gapSize == 0){
				gapStartSiteInSeq = pos - gapOffSet;
			}
			++gapSize;
			++gapOffSet;
		} else {
			//log gap if currently building one
			if(gapSize > 0){
				gapInfos.emplace_back(gapStartSiteInSeq, gapSize);
			}
			//reset
			gapSize = 0;
			gapStartSiteInSeq = 0;
		}
		++pos;
	}
	if(gapSize > 0){
		gapInfos.emplace_back(gapStartSiteInSeq, gapSize);
	}
	//sort backwards and add appropriate number of zeros, starting from back so pos is not invalidated
	std::sort(gapInfos.begin(),
			gapInfos.end(), [](const std::pair<uint32_t, uint32_t> & p1,
					const std::pair<uint32_t, uint32_t> & p2){ return p1.first > p2.first;});
	std::vector<gapInfo> gInfos;
	for(const auto & g : gapInfos){
		gInfos.emplace_back(gapInfo(g.first, g.second, false));
	}
	return alnInfoGlobal(gInfos, 1, false);
}


Muscler::Muscler(const bfs::path & musclePath) {
	setMusclePath(musclePath);
}

Muscler::Muscler() :
		Muscler("muscle") {
}

void Muscler::setMusclePath(const bfs::path & musclePath) {
	musclePath_ = musclePath;
	auto hasProgram = bib::sys::hasSysCommand(musclePath_.string());
	if (!hasProgram) {
		std::stringstream ss;
		ss << bib::bashCT::boldBlack(musclePath_.string())
				<< bib::bashCT::boldRed(
						" is not in path or may not be executable, cannot be used")
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
}

std::vector<seqInfo> Muscler::muscleSeqs(const SeqIOOptions & opts){
	SeqInput reader(opts);
	auto seqs = reader.readAllReads<seqInfo>();
	muscleSeqs(seqs);
	return seqs;
}



Muscler::MusPosSize::MusPosSize() :
		MusPosSize(std::numeric_limits < size_t > ::max()) {
}

Muscler::MusPosSize::MusPosSize(size_t pos) :
		pos_(pos) {
}

Muscler::MusPosSize::MusPosSize(size_t pos, uint32_t size) :
		pos_(pos), size_(size) {
}



Muscler::AlnPosScore::AlnPosScore(size_t pos) :
	pos_(pos) {
}


void Muscler::AlnPosScore::setCounts() {
	counter_.resetAlphabet(true);
	counter_.setFractions();
	for (const auto & let : counter_.alphabet_) {
		if ('-' == let) {
			gapCount_ += counter_.chars_[let];
		} else {
			baseCount_ += counter_.chars_[let];
		}
	}
}


double Muscler::AlnPosScore::getBaseSpannedPerc() const{
	return static_cast<double>(baseCount_)/getSpanningCount();
}

double Muscler::AlnPosScore::getPercentOfSequencesSpanningPosition(uint32_t totalInputSeqs) const{
	return static_cast<double>(getSpanningCount())/totalInputSeqs;
}

uint32_t Muscler::AlnPosScore::getSpanningCount() const{
	return baseCount_ + gapCount_;
}

Muscler::AlnPosScoreStreak::AlnPosScoreStreak() :
		start_(std::numeric_limits<uint32_t>::max()), end_(
				std::numeric_limits<uint32_t>::max()) {

}
Muscler::AlnPosScoreStreak::AlnPosScoreStreak(const std::shared_ptr<AlnPosScore> & firstScore) :
		start_(firstScore->pos_), end_(firstScore->pos_ + 1), scores_{firstScore} {

}

uint32_t Muscler::AlnPosScoreStreak::getLen() const{
	return end_ - start_;
}


std::vector<Muscler::AlnPosScoreStreak> Muscler::getAlignmentStreaksPositions(const std::vector<std::shared_ptr<AlnPosScore>> & scores,
		const std::function<bool(const std::shared_ptr<AlnPosScore>&)> & scorePred,
		uint32_t streakLenCutOff){
	std::vector<AlnPosScoreStreak> streaks;
	if(scores.empty()){
		return streaks;
	}

	auto startCut = std::find_if(scores.begin(), scores.end(), scorePred);
	auto endCut = std::find_if(scores.rbegin(), scores.rend(), scorePred);
	uint32_t startCutSeqPos =  startCut - scores.begin();
	uint32_t stopCutSeqPos  = len(scores) - (endCut - scores.rbegin()) - 1;

	AlnPosScoreStreak currentStreak(scores[startCutSeqPos]);
	for(const auto & streakSearchPos : iter::range(startCutSeqPos + 1, stopCutSeqPos + 1)){
		if(scorePred(scores[streakSearchPos])){
			if(currentStreak.end_ == scores[streakSearchPos]->pos_){
				++currentStreak.end_;
				currentStreak.scores_.emplace_back(scores[streakSearchPos]);
			}else{
				if(currentStreak.getLen() >= streakLenCutOff){
					streaks.emplace_back(currentStreak);
				}
				currentStreak = AlnPosScoreStreak(scores[streakSearchPos]);
			}
		}
	}
	if(currentStreak.getLen() >= streakLenCutOff){
		streaks.emplace_back(currentStreak);
	}
	return streaks;
}

} /* namespace bibseq */
