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

} /* namespace bibseq */
