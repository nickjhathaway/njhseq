/*
 * tandemRepeatUtils.cpp
 *
 *  Created on: May 31, 2020
 *      Author: nick
 */



#include "tandemRepeatUtils.hpp"

namespace njhseq {

SimpleTandemRepeatFinder::RepeatStartLenRepeats SimpleTandemRepeatFinder::getRepeatInfo(
	const std::string& seq, const motif& m) {
	auto locs = m.findPositionsFull(seq, 0);
	njh::sort(locs);
	uint32_t maxRepeatNumber = 0;
	uint32_t start = std::numeric_limits<uint32_t>::max();
	if (!locs.empty()) {
		uint32_t length = 1;
		start = locs.front();
		for (const auto pos: iter::range<uint32_t>(1, locs.size())) {
			if (locs[pos] == locs[pos - 1] + m.size()) {
				++length;
			} else {
				if (length > maxRepeatNumber) {
					maxRepeatNumber = length;
				}
				length = 1;
				start = locs[pos];
			}
		}
		if (length > maxRepeatNumber) {
			maxRepeatNumber = length;
		}
	}
	return RepeatStartLenRepeats{start, static_cast<uint32_t>(maxRepeatNumber * m.size()), maxRepeatNumber};
}

void SimpleTandemRepeatFinder::SimpleTRFinderLocsPars::setDefaultOpts(seqSetUp & setUp){
	setUp.setOption(numThreads, "--numThreads", "Number of CPUs to use", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.setOption(lengthCutOff, "--lengthCutOff", "The minimum length to report(inclusive)");
	setUp.setOption(minNumRepeats, "--minNumRepeats", "The minimum number of repeated units(inclusive)");
	setUp.setOption(minRepeatUnitSize, "--minRepeatUnitSize", "Min Repeat Unit Size to search for(inclusive)");
	setUp.setOption(maxRepeatUnitSize, "--maxRepeatUnitSize", "Max Repeat Unit Size to search for(inclusive)");
	setUp.setOption(searchAllUnits, "--searchAllUnits", "Search All Units, by default similar tandems are collapsed e.g. searching only for AT rather than AT,TA,ATAT etc");
	setUp.setOption(doNotAddFlankingSeq, "--doNotAddFlankingSeq", "Don't Add Flanking Seq, by default partially matching before and after sequence is added e.g. adding the precedding GCT if the unit search is AGCT for ATGCT-AGCT-AGCT-AGCT");
	setUp.setOption(addFullSeqToOuput, "--addFullSeqToOuput", "Add Full Seq To Ouput");


	if(minRepeatUnitSize > maxRepeatUnitSize){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("maxRepeatUnitSize: ", maxRepeatUnitSize, " can't be greater than minRepeatUnitSize: ", minRepeatUnitSize));
	}
	setUp.setOption(alphabet, "--alphabet", "alphabet");

	setUp.processWritingOptions(outOpts);
}


SimpleTandemRepeatFinder::SimpleTandemRepeatFinder(SimpleTRFinderLocsPars pars): pars_(std::move(pars)){

}

SimpleTandemRepeatFinder::MinimalUnitsAndAltMotifs SimpleTandemRepeatFinder::genMinimalUnitsNeededForSearch() const {
	MinimalUnitsAndAltMotifs ret;
	auto allUnits = genAllUnitsPossible();
		for(const auto & unit : allUnits){
			bool add = true;
			if(unit.size() > 1){
				for(const auto & otherUnit : *ret.allUnits){
					//check rotating if the same size
					if(otherUnit.size() == unit.size()){
						if(!checkTwoRotatingStrings(unit, otherUnit, 0).empty()){
							ret.altMots[otherUnit].emplace_back(unit);
							add = false;
							break;
						}
					}else if(unit.size() > otherUnit.size() &&  unit.size() % otherUnit.size() == 0){
						motif otherUnitMot(otherUnit);
						auto locs = otherUnitMot.findPositionsFull(unit, 0);
						if(!locs.empty()){
							njh::sort(locs);
							uint32_t length = 1;
							size_t start = locs.front();
							for(const auto pos : iter::range<uint32_t>(1, locs.size())){
								if(locs[pos] == locs[pos - 1] + otherUnitMot.size() ){
									++length;
								} else {
									length = 1;
									start = locs[pos];
								}
							}
							uint32_t outputStart = start;
							uint32_t outputEnd = start + otherUnitMot.size() * length;

							if(start > 0){
								uint32_t walkbackPos = 0;
								while(outputStart > 0 && walkbackPos + 1 < otherUnitMot.size()){
									if(unit[outputStart - 1] == otherUnit[otherUnitMot.size() - walkbackPos - 1] ){
										++walkbackPos;
										--outputStart;
									}else{
										break;
									}
								}
							}
							if(outputEnd != unit.size()){
								uint32_t walkforwardPos = 0;
								while(outputEnd < unit.size() && walkforwardPos < otherUnit.size()){
									if(unit[outputEnd] == otherUnit[walkforwardPos]){
										++outputEnd;
										++walkforwardPos;
									}else{
										break;
									}
								}
							}
//							if("AT" == otherUnit && "TATATA" == unit){
//								std::cout << __FILE__ << " " << __LINE__ << std::endl;
//								std::cout << "outputStart: " << outputStart << std::endl;
//								std::cout << "outputEnd:   " << outputEnd << std::endl;
//								std::cout << "unit.size(): " << unit.size() << std::endl;
//								exit(1);
//							}
							if(0 == outputStart && unit.size() == outputEnd){
								add = false;
								break;
							}
						}
					}
				}
			}
			if(add){
				ret.allUnits->emplace_back(unit);
			}
		}
	return ret;
}
VecStr SimpleTandemRepeatFinder::genAllUnitsPossible() const {
	VecStr allUnits;
	for(uint32_t repeatUnitSize = pars_.minRepeatUnitSize; repeatUnitSize < pars_.maxRepeatUnitSize + 1; ++repeatUnitSize){
		auto units = permuteVector(pars_.alphabet, repeatUnitSize);
		for(const auto & unit : units){
			allUnits.emplace_back(njh::pasteAsStr(unit));
		}
	}
	return allUnits;
}

void SimpleTandemRepeatFinder::runSimpleTRFinderLocs(const SeqIOOptions & seqInput) const{

	//first create tandems that will be searched for
	MinimalUnitsAndAltMotifs minimalUnits;
	if (pars_.searchAllUnits) {
		minimalUnits.allUnits = std::make_shared<VecStr>(genAllUnitsPossible());
	} else {
		minimalUnits = genMinimalUnitsNeededForSearch();
	}

	if(pars_.debug){
		std::cout << "Searching for: " << std::endl;
		printVector(*minimalUnits.allUnits,"\n");
	}


	seqInfo seq;
	SeqInput reader(seqInput);
	reader.openIn();
	OutputStream out(pars_.outOpts);
	std::mutex outMut;
	while(reader.readNextRead(seq)){
		njh::concurrent::LockableQueue<std::string> unitQueue(*minimalUnits.allUnits);
		std::function<void()> findTandems = [&unitQueue,&seq,this,&out,&outMut,
																				 &minimalUnits](){
			std::string motifstr;
			std::vector<Bed6RecordCore> repeatUnitLocs;
			std::stringstream currentOut;
			while(unitQueue.getVal(motifstr)){
				motif mot(motifstr);
				uint32_t currentAllowableError = 0;
				auto locs = mot.findPositionsFull(seq.seq_, currentAllowableError);
				njh::sort(locs);
				if(!locs.empty()){
					uint32_t length = 1;
					size_t start = locs.front();
					for(const auto pos : iter::range<uint32_t>(1, locs.size())){
						if(locs[pos] == locs[pos - 1] + mot.size() ){
							++length;
						} else {
							uint32_t outputStart = start;
							uint32_t outputEnd = start + mot.size() * length;
							if(!pars_.doNotAddFlankingSeq){
								if(start > 0){
									uint32_t walkbackPos = 0;
									while(outputStart > 0 && walkbackPos + 1 < mot.size()){
										if(seq.seq_[outputStart - 1] == motifstr[mot.size() - walkbackPos - 1] ){
											++walkbackPos;
											--outputStart;
										}else{
											break;
										}
									}
								}
								if(outputEnd != seq.seq_.size()){
									uint32_t walkforwardPos = 0;
									while(outputEnd < seq.seq_.size() && walkforwardPos < motifstr.size()){
										if(seq.seq_[outputEnd] == motifstr[walkforwardPos]){
											++outputEnd;
											++walkforwardPos;
										}else{
											break;
										}
									}
								}
							}
							//check to see that there is at least one alt tandem that equals the min number of required repeats
							uint32_t repeatNumber = length;
							std::string outMotifstr = motifstr;
							//					if(repeatNumber < minNumRepeats && (outputEnd - outputStart)/mot.size() >= minNumRepeats && njh::in(motifstr, altMots)){
							//						auto subSeq = seq.seq_.substr(outputStart, outputEnd - outputStart);
							//						for(const auto & altMot : altMots.at(motifstr)){
							//							auto altRepeatInfo = getRepeatInfo(subSeq, altMot);
							//							if (altRepeatInfo.repeats_ > repeatNumber) {
							//								repeatNumber = altRepeatInfo.repeats_;
							//								outMotifstr = altMot.motifOriginal_;
							//							} else if (altRepeatInfo.repeats_ == repeatNumber && outputStart + altRepeatInfo.start_ < start) {
							//								outMotifstr = altMot.motifOriginal_;
							//							}
							//						}
							//					}
							if(outputEnd - outputStart >=pars_.lengthCutOff && njh::in(motifstr, minimalUnits.altMots) ){
								auto subSeq = seq.seq_.substr(outputStart, outputEnd - outputStart);
								auto minStart = start;
								for(const auto & altMot : minimalUnits.altMots.at(motifstr)){
									auto altRepeatInfo = getRepeatInfo(subSeq, altMot);
									if (altRepeatInfo.repeats_ > repeatNumber) {
										repeatNumber = altRepeatInfo.repeats_;
										outMotifstr = altMot.motifOriginal_;
									} else if (altRepeatInfo.repeats_ == repeatNumber && outputStart + altRepeatInfo.start_ < minStart) {
										outMotifstr = altMot.motifOriginal_;
										minStart = outputStart + altRepeatInfo.start_;
									}
								}
							}

							if(repeatNumber >= pars_.minNumRepeats && outputEnd - outputStart >=pars_.lengthCutOff){
								currentOut << seq.name_
										<< "\t" << outputStart
										<< "\t" << outputEnd
										<< "\t" << njh::pasteAsStr(seq.name_,"-", outputStart, "-", outputEnd) << "__" << outMotifstr << "_x" << static_cast<double>(outputEnd - outputStart)/static_cast<double>(mot.size())
										<< "\t" << outputEnd - outputStart
										<< "\t" << "+";
								if(pars_.addFullSeqToOuput){
									currentOut << "\t" << seq.seq_.substr(outputStart, outputEnd - outputStart);
								}
								currentOut << '\n';
							}
							length = 1;
							start = locs[pos];
						}
					}
					uint32_t outputStart = start;
					uint32_t outputEnd = start + mot.size() * length;
					if(!pars_.doNotAddFlankingSeq){
						if(start > 0){
							uint32_t walkbackPos = 0;
							while(outputStart > 0 && walkbackPos + 1 < mot.size()){
								if(seq.seq_[outputStart - 1] == motifstr[mot.size() - walkbackPos - 1] ){
									++walkbackPos;
									--outputStart;
								}else{
									break;
								}
							}
						}
						if(outputEnd != seq.seq_.size()){
							uint32_t walkforwardPos = 0;
							while(outputEnd < seq.seq_.size() && walkforwardPos < motifstr.size()){
								if(seq.seq_[outputEnd] == motifstr[walkforwardPos]){
									++outputEnd;
									++walkforwardPos;
								}else{
									break;
								}
							}
						}
					}
					//check to see that there is at least one alt tandem that equals the min number of required repeats
					uint32_t repeatNumber = length;
					std::string outMotifstr = motifstr;
//					if(repeatNumber < minNumRepeats && (outputEnd - outputStart)/mot.size() >= minNumRepeats && njh::in(motifstr, altMots)){
//						auto subSeq = seq.seq_.substr(outputStart, outputEnd - outputStart);
//						for(const auto & altMot : altMots.at(motifstr)){
//							auto altRepeatInfo = getRepeatInfo(subSeq, altMot);
//							if (altRepeatInfo.repeats_ > repeatNumber) {
//								repeatNumber = altRepeatInfo.repeats_;
//								outMotifstr = altMot.motifOriginal_;
//							} else if (altRepeatInfo.repeats_ == repeatNumber && outputStart + altRepeatInfo.start_ < start) {
//								outMotifstr = altMot.motifOriginal_;
//							}
//						}
//					}
					if(outputEnd - outputStart >=pars_.lengthCutOff && njh::in(motifstr, minimalUnits.altMots)){
						auto subSeq = seq.seq_.substr(outputStart, outputEnd - outputStart);
						auto minStart = start;
						for(const auto & altMot : minimalUnits.altMots.at(motifstr)){
							auto altRepeatInfo = getRepeatInfo(subSeq, altMot);
							if (altRepeatInfo.repeats_ > repeatNumber) {
								repeatNumber = altRepeatInfo.repeats_;
								outMotifstr = altMot.motifOriginal_;
							} else if (altRepeatInfo.repeats_ == repeatNumber && outputStart + altRepeatInfo.start_ < minStart) {
								outMotifstr = altMot.motifOriginal_;
								minStart = outputStart + altRepeatInfo.start_;
							}
						}
					}
					if(repeatNumber >= pars_.minNumRepeats && outputEnd - outputStart >=pars_.lengthCutOff){
						currentOut << seq.name_
								<< "\t" << outputStart
								<< "\t" << outputEnd
								<< "\t" << njh::pasteAsStr(seq.name_,"-", outputStart, "-", outputEnd) << "__" << outMotifstr << "_x" << static_cast<double>(outputEnd - outputStart)/static_cast<double>(mot.size())
								<< "\t" << outputEnd - outputStart
								<< "\t" << "+";
						if(pars_.addFullSeqToOuput){
							currentOut << "\t" << seq.seq_.substr(outputStart, outputEnd - outputStart);
						}
						currentOut << '\n';
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(outMut);
				out << currentOut.str();
			}
		};
		njh::concurrent::runVoidFunctionThreaded(findTandems, pars_.numThreads);
	}
}


}  // namespace njhseq
