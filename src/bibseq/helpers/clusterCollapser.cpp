#include "clusterCollapser.hpp"

namespace bibseq {

std::vector<identicalCluster> clusterCollapser::collapseIdenticalReads(
    const std::vector<readObject> &inputReads, const std::string &repQual,
    const std::string &lower) {
  std::vector<identicalCluster> finalClusterVec;
  // handelLowerCaseBases(inputReads, lower);
  int count = 0;
  for (const auto &seqIter : inputReads) {
    ++count;
    if (count == 1) {
      finalClusterVec.emplace_back(identicalCluster(seqIter));
      continue;
    }
    bool foundMatch = false;
    for (auto &clusterIter : finalClusterVec) {
      if (seqIter.seqBase_.seq_ == clusterIter.seqBase_.seq_) {
        clusterIter.addRead(seqIter);
        foundMatch = true;
        break;
      }
    }
    if (!foundMatch) {
      finalClusterVec.emplace_back(identicalCluster(seqIter));
    }
  }
  identicalCluster::setIdneticalClusterQual(finalClusterVec, repQual);
  readVec::allSetLetterCount(finalClusterVec);
  readVec::allUpdateName(finalClusterVec);
  return finalClusterVec;
}

std::vector<cluster> clusterCollapser::markChimeras(
    std::vector<cluster> &processedReads, aligner &alignerObj,
    double parentFreqs, bool kmersByPos, int kLength, int runCutOff, bool local,
    bool weighHomopolyer) {
  std::vector<cluster> output;
  kmerCalculator kCalc = kmerCalculator();
  kmerMaps chiMaps =
      kCalc.indexKmerMpas(processedReads, kLength, runCutOff, runCutOff);
  auto oldMaps = alignerObj.getKmerMaps();
  alignerObj.setKmerMpas(chiMaps);

  int counter = 0;
  for (std::vector<cluster>::iterator clusIter = processedReads.begin();
       clusIter != processedReads.end(); ++clusIter) {
    if (clusIter == processedReads.begin()) {
      output.push_back(*clusIter);
      continue;
    }
    ++counter;
    // size_t furestMismatch = clusIter->seqBase_.seq_.length() / 4;
    size_t furestMismatch = 0;
    bool foundMatch = false;
    bool foundChimera = false;
    for (std::vector<cluster>::iterator secondClusIter = output.begin();
         secondClusIter != output.end(); ++secondClusIter) {
      if (secondClusIter == clusIter) {
        continue;
      }
      if (parentFreqs * clusIter->seqBase_.frac_ >
          secondClusIter->seqBase_.frac_) {
        continue;
      }
      // alignerObj.align(<#const simpleReadObject &firstRead#>, <#const
      // simpleReadObject &secondRead#>, <#bool usingQuality#>, <#bool local#>,
      // <#bool simple#>)
      alignerObj.alignVec(*clusIter, *secondClusIter, local);
      alignerObj.profileAlignment(*clusIter, *secondClusIter, kLength,
                                  kmersByPos, true, true, false,
                                  weighHomopolyer);
      /*
       if (kmersByPos) {
       alignerObj.profileAlignment(clusIter->seqBase_.seq_,
       secondClusIter->seqBase_.seq_,
       kmerMapByPos, kLength, runCutOff, true, false);
       }else{
       alignerObj.profileAlignment(clusIter->seqBase_.seq_,
       secondClusIter->seqBase_.seq_,
       kmerMap, kLength, runCutOff, true, false);
       }*/

      std::string tailEnd = "";
      if (alignerObj.comp_.hqMismatches_ > 1 &&
          alignerObj.comp_.largeBaseIndel_ < 2 &&
          alignerObj.alignmentGaps_.begin()->second.size_ < 50) {
        // if (alignerObj.comparison.highQualityMismatch>1 &&
        // clusIter->seqBase_.seq_.length()==secondClusIter->seqBase_.seq_.length()
        // &&
        // alignerObj.comparison.numberOfLargeGaps<1) {
        // if (alignerObj.errors_.hqMismatches>1) {
        size_t firstMismatchPos =
            alignerObj.mismatches_.begin()->second.refBasePos;
        if (firstMismatchPos > furestMismatch) {
          furestMismatch = firstMismatchPos;
          clusIter->chimeras.clear();
          clusIter->endChimeras.clear();
          /*
          std::cout << std::endl << ">"<< clusIter->name << std::endl;
          std::cout << clusIter->seqBase_.seq_ << std::endl;
          std::cout << std::endl << ">"<< secondClusIter->name << std::endl;
          std::cout << secondClusIter->seqBase_.seq_ << std::endl;
          std::cout << firstMismatchPos << std::endl;
          std::cout << alignerObj.mismatches.begin()->second.refBasePos <<
          std::endl;
          std::cout <<
          alignerObj.mismatches.begin()->second.seqBase_.seq_BasePos <<
          std::endl;*/
          tailEnd = clusIter->seqBase_.seq_.substr(firstMismatchPos);
          // std::cout<<"tailEnd"<<tailEnd<<std::endl;
          // tailEnd=alignerObj.alignObjectA.seqBase_.seq_.substr(firstMismatchPos);
          for (auto &third : processedReads) {
            if (third.seqBase_.name_ == secondClusIter->seqBase_.name_ ||
                third.seqBase_.name_ == clusIter->seqBase_.name_) {
              continue;
            }
            if (third.seqBase_.seq_.rfind(tailEnd) == std::string::npos) {
              continue;
            }
            if (parentFreqs * clusIter->seqBase_.frac_ > third.seqBase_.frac_) {
              continue;
            }
            if ((third.seqBase_.seq_.rfind(tailEnd) + tailEnd.size()) ==
                third.seqBase_.seq_.size()) {

              // clusIter->name=combineStrings("CHI_",clusIter->name);
              // secondClusIter->addCluster(*clusIter);
              clusIter->seqBase_.markAsChimeric();
              clusIter->endChimeras.push_back(third);
              foundMatch = false;
              foundChimera = true;
              break;
            }
          }
          clusIter->chimeras.push_back(*secondClusIter);
        }
      }
      if (foundChimera) {
        break;
      }
    }
    if (!foundMatch) {
      output.push_back(*clusIter);
    }
  }

  alignerObj.setKmerMpas(oldMaps);

  return output;
}
void markChimerasAdvancedNew(
    std::vector<cluster> &processedReads, aligner &alignerObj,
    double parentFreqs, double runCutOff,
    const comparison &chiOverlap, uint32_t overLapSizeCutoff,
    bool weighHomopolyer, uint32_t &chimeraCount, uint32_t allowableError, bool verbose) {

  for (const auto &pos : iter::range(processedReads.size())) {
  	if(verbose){
    	std::cout << pos << ":" << processedReads.size()<< "\r";
    	std::cout.flush();
  	}
    for (const auto &subPos : iter::range(pos + 1, processedReads.size())) {
    	if(pos == subPos){
    		continue;
    	}
    	//std::cout << "part 1" << std::endl;
      if ((processedReads[pos].seqBase_.frac_ /
           processedReads[subPos].seqBase_.frac_) < parentFreqs ||
          processedReads[subPos].seqBase_.cnt_ <= runCutOff) {
        // std::cout << processedReads[pos].seqBase_.frac_ /
        // processedReads[subPos].seqBase_.frac_ << std::endl;
      	//std::cout << "part 1.5" << std::endl;
        continue;
      }
      //std::cout << "part 2" << std::endl;
      alignerObj.alignVec(processedReads[pos], processedReads[subPos], false);
      // alignerObj.profilePrimerAlignment(processedReads[pos],
      // processedReads[subPos], setUp.weightHomopolymers_);
      alignerObj.profileAlignment(processedReads[pos], processedReads[subPos],
                                  11, true, false, true, false,
                                  weighHomopolyer);
      /*if (alignerObj.errors_.largeBaseIndel_ > 0) {
        continue;
      }*/
      // std::cout << alignerObj.mismatches_.size() << std::endl;
      if (alignerObj.mismatches_.size() > allowableError) {
      	//std::cout << "part 3" << std::endl;
      	//std::cout << processedReads[subPos].seqBase_.name_ << std::endl;
      	//std::cout <<  processedReads[pos].seqBase_.name_ << std::endl;
      	//std::cout << "first mismatch " <<  alignerObj.mismatches_.begin()->first << std::endl;
      	//std::cout << "last mismatch " << alignerObj.mismatches_.rbegin()->first << std::endl;
      	//std::cout << "front overlap size " << alignerObj.mismatches_.begin()->first + 1 << std::endl;
      	//std::cout << "end overlap size " << alignerObj.alignObjectA_.seqBase_.seq_.size() -
            //alignerObj.mismatches_.rbegin()->first << std::endl;
      	//std::cout << "overlap cut off : " << overLapSizeCutoff << std::endl;
      	bool failedEnd = true;
      	bool failedFront = true;
      	std::map<uint32_t, mismatch> savedMismatches = alignerObj.mismatches_;
        if (savedMismatches.begin()->first + 1 <= overLapSizeCutoff){
        	failedFront = true;
        	//std::cout << "failed front due to overlap" << std::endl;
        }else{
        	//std::cout << "this part 1" << std::endl;
          alignerObj.profileAlignment(
              processedReads[pos], processedReads[subPos], 11, true, false, true,
              false, weighHomopolyer, 0, savedMismatches.begin()->first);
          //std::cout << "that part 1" << std::endl;
          failedFront =
              chiOverlap.passErrorProfile(alignerObj.comp_);
        }
        //std::cout << "here?" << std::endl;
        if(overLapSizeCutoff >=
            (alignerObj.alignObjectA_.seqBase_.seq_.size() -
            		savedMismatches.rbegin()->first)){
        	failedEnd = true;
        	//std::cout << "failed end due to overlap" << std::endl;
        }else{
         //std::cout << "this part 2" << std::endl;
          alignerObj.profileAlignment(
              processedReads[pos], processedReads[subPos], 11, true, false, true,
              false, weighHomopolyer, savedMismatches.rbegin()->first + 1);
          //std::cout << "that part 2" << std::endl;
          // std::cout << "after these? " << std::endl;
          failedEnd = chiOverlap.passErrorProfile(alignerObj.comp_);

        }
        //std::cout << "part 4" << std::endl;
        //std::cout << "this part" << std::endl;
        if (!failedEnd && !failedFront) {
          continue;
        }
        // processedReads[subPos].outPutSeq(std::cout);
        // processedReads[pos].outPutSeq(std::cout);
        if (failedFront) {
          processedReads[subPos].frontChiPos.insert(
              {savedMismatches.begin()->second.seqBasePos, pos});
          // baseReadObject subReadAFront =
          // baseReadObject(alignerObj.alignObjectA_.seqBase_.getSubRead(0,
          // savedMismatches.begin()->first));
          // baseReadObject subReadBFront =
          // baseReadObject(alignerObj.alignObjectB_.seqBase_.getSubRead(0,
          // savedMismatches.begin()->first));
          // subReadAFront.outPutSeq(std::cout);
          // subReadBFront.outPutSeq(std::cout);
        }
        if (failedEnd) {
          processedReads[subPos].endChiPos.insert(
              {savedMismatches.rbegin()->second.seqBasePos, pos});
          // baseReadObject subReadABack =
          // baseReadObject(alignerObj.alignObjectA_.seqBase_.getSubRead(savedMismatches.rbegin()->first
          // + 1));
          // baseReadObject subReadBBack =
          // baseReadObject(alignerObj.alignObjectB_.seqBase_.getSubRead(savedMismatches.rbegin()->first
          // + 1));
          // subReadBBack.outPutSeq(std::cout);
          // subReadABack.outPutSeq(std::cout);
        }
      }
    }
  }
  if(verbose){
  	std::cout << std::endl;
  }
  bool print = false;
  if (print) {
    for (auto &clus : processedReads) {
      std::cout << clus.seqBase_.name_;
      for (const auto &front : iter::reverse(clus.frontChiPos)) {
        std::cout << "\t" << front.first << "\t"
                  << processedReads[front.second].seqBase_.name_;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for (auto &clus : processedReads) {
      std::cout << clus.seqBase_.name_;
      for (const auto &end : clus.endChiPos) {
        std::cout << "\t" << end.first << "\t"
                  << processedReads[end.second].seqBase_.name_;
      }
      std::cout << std::endl;
    }
  }
  for (auto &clus : processedReads) {
    bool foundChi = false;
    for (const auto &front : iter::reverse(clus.frontChiPos)) {
      for (const auto &pos : iter::range(processedReads.size())) {
        if ((processedReads[pos].seqBase_.frac_ / clus.seqBase_.frac_) <
                parentFreqs ||
            pos == front.second) {
          continue;
        }
        alignerObj.alignVec(processedReads[pos], clus, false);
        alignerObj.profileAlignment(processedReads[pos], clus, 11, true, false,
                                    true, false, weighHomopolyer, front.first);
        if (chiOverlap.passErrorProfile(alignerObj.comp_)) {
          foundChi = true;
          clus.seqBase_.markAsChimeric();
          ++chimeraCount;
          break;
        }
      }
      if (foundChi) {
        break;
      }
    }
    if (foundChi) {
      continue;
    }
    for (const auto &end : clus.endChiPos) {
      for (const auto &pos : iter::range(processedReads.size())) {
        if ((processedReads[pos].seqBase_.frac_ / clus.seqBase_.frac_) <
                parentFreqs ||
            pos == end.second) {
          continue;
        }
        alignerObj.alignVec(processedReads[pos], clus, false);
        alignerObj.profileAlignment(processedReads[pos], clus, 11, true, false,
                                    true, false, weighHomopolyer, 0,
                                    end.first + 1);
        if (chiOverlap.passErrorProfile(alignerObj.comp_)) {
          foundChi = true;
          clus.seqBase_.markAsChimeric();
          ++chimeraCount;
          break;
        }
      }
      if (foundChi) {
        break;
      }
    }
  }
}


void clusterCollapser::markChimerasAdvanced(
    std::vector<cluster> &processedReads, aligner &alignerObj,
    double parentFreqs, int runCutOff, bool local,
    const comparison &chiOverlap, uint32_t overLapSizeCutoff,
    bool weighHomopolyer, uint32_t &chimeraCount, uint32_t allowableError) {

  for (const auto &pos : iter::range(processedReads.size())) {
  	std::cout << pos << ":" << processedReads.size()<< "\r";
  	std::cout.flush();
    for (const auto &subPos : iter::range(pos + 1, processedReads.size())) {
    	if(pos == subPos){
    		continue;
    	}
    	//std::cout << "part 1" << std::endl;
      if ((processedReads[pos].seqBase_.frac_ /
           processedReads[subPos].seqBase_.frac_) < parentFreqs ||
          processedReads[subPos].seqBase_.cnt_ <= runCutOff) {
        // std::cout << processedReads[pos].seqBase_.frac_ /
        // processedReads[subPos].seqBase_.frac_ << std::endl;
      	//std::cout << "part 1.5" << std::endl;
        continue;
      }
      //std::cout << "part 2" << std::endl;
      alignerObj.alignVec(processedReads[pos], processedReads[subPos], local);
      // alignerObj.profilePrimerAlignment(processedReads[pos],
      // processedReads[subPos], setUp.weightHomopolymers_);
      alignerObj.profileAlignment(processedReads[pos], processedReads[subPos],
                                  11, true, false, true, false,
                                  weighHomopolyer);
      /*if (alignerObj.errors_.largeBaseIndel_ > 0) {
        continue;
      }*/
      // std::cout << alignerObj.mismatches_.size() << std::endl;
      if (alignerObj.mismatches_.size() > allowableError) {
      	//std::cout << "part 3" << std::endl;
      	//std::cout << processedReads[subPos].seqBase_.name_ << std::endl;
      	//std::cout <<  processedReads[pos].seqBase_.name_ << std::endl;
      	//std::cout << "first mismatch " <<  alignerObj.mismatches_.begin()->first << std::endl;
      	//std::cout << "last mismatch " << alignerObj.mismatches_.rbegin()->first << std::endl;
      	//std::cout << "front overlap size " << alignerObj.mismatches_.begin()->first + 1 << std::endl;
      	//std::cout << "end overlap size " << alignerObj.alignObjectA_.seqBase_.seq_.size() -
            //alignerObj.mismatches_.rbegin()->first << std::endl;
      	//std::cout << "overlap cut off : " << overLapSizeCutoff << std::endl;
      	bool failedEnd = true;
      	bool failedFront = true;
      	std::map<uint32_t, mismatch> savedMismatches = alignerObj.mismatches_;
        if (savedMismatches.begin()->first + 1 <= overLapSizeCutoff){
        	failedFront = true;
        	//std::cout << "failed front due to overlap" << std::endl;
        }else{
        	//std::cout << "this part 1" << std::endl;
          alignerObj.profileAlignment(
              processedReads[pos], processedReads[subPos], 11, true, false, true,
              false, weighHomopolyer, 0, savedMismatches.begin()->first);
          //std::cout << "that part 1" << std::endl;
          failedFront =
              chiOverlap.passErrorProfile(alignerObj.comp_);
        }
        //std::cout << "here?" << std::endl;
        if(overLapSizeCutoff >=
            (alignerObj.alignObjectA_.seqBase_.seq_.size() -
            		savedMismatches.rbegin()->first)){
        	failedEnd = true;
        	//std::cout << "failed end due to overlap" << std::endl;
        }else{
         //std::cout << "this part 2" << std::endl;
          alignerObj.profileAlignment(
              processedReads[pos], processedReads[subPos], 11, true, false, true,
              false, weighHomopolyer, savedMismatches.rbegin()->first + 1);
          //std::cout << "that part 2" << std::endl;
          // std::cout << "after these? " << std::endl;
          failedEnd = chiOverlap.passErrorProfile(alignerObj.comp_);

        }
        //std::cout << "part 4" << std::endl;
        //std::cout << "this part" << std::endl;
        if (!failedEnd && !failedFront) {
          continue;
        }
        // processedReads[subPos].outPutSeq(std::cout);
        // processedReads[pos].outPutSeq(std::cout);
        if (failedFront) {
          processedReads[subPos].frontChiPos.insert(
              {savedMismatches.begin()->second.seqBasePos, pos});
          // baseReadObject subReadAFront =
          // baseReadObject(alignerObj.alignObjectA_.seqBase_.getSubRead(0,
          // savedMismatches.begin()->first));
          // baseReadObject subReadBFront =
          // baseReadObject(alignerObj.alignObjectB_.seqBase_.getSubRead(0,
          // savedMismatches.begin()->first));
          // subReadAFront.outPutSeq(std::cout);
          // subReadBFront.outPutSeq(std::cout);
        }
        if (failedEnd) {
          processedReads[subPos].endChiPos.insert(
              {savedMismatches.rbegin()->second.seqBasePos, pos});
          // baseReadObject subReadABack =
          // baseReadObject(alignerObj.alignObjectA_.seqBase_.getSubRead(savedMismatches.rbegin()->first
          // + 1));
          // baseReadObject subReadBBack =
          // baseReadObject(alignerObj.alignObjectB_.seqBase_.getSubRead(savedMismatches.rbegin()->first
          // + 1));
          // subReadBBack.outPutSeq(std::cout);
          // subReadABack.outPutSeq(std::cout);
        }
      }
    }
  }
  std::cout << std::endl;
  bool print = false;
  if (print) {
    for (auto &clus : processedReads) {
      std::cout << clus.seqBase_.name_;
      for (const auto &front : iter::reverse(clus.frontChiPos)) {
        std::cout << "\t" << front.first << "\t"
                  << processedReads[front.second].seqBase_.name_;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for (auto &clus : processedReads) {
      std::cout << clus.seqBase_.name_;
      for (const auto &end : clus.endChiPos) {
        std::cout << "\t" << end.first << "\t"
                  << processedReads[end.second].seqBase_.name_;
      }
      std::cout << std::endl;
    }
  }
  for (auto &clus : processedReads) {
    bool foundChi = false;
    for (const auto &front : iter::reverse(clus.frontChiPos)) {
      for (const auto &pos : iter::range(processedReads.size())) {
        if ((processedReads[pos].seqBase_.frac_ / clus.seqBase_.frac_) <
                parentFreqs ||
            pos == front.second) {
          continue;
        }
        alignerObj.alignVec(processedReads[pos], clus, local);
        alignerObj.profileAlignment(processedReads[pos], clus, 11, true, false,
                                    true, false, weighHomopolyer, front.first);
        if (chiOverlap.passErrorProfile(alignerObj.comp_)) {
          foundChi = true;
          clus.seqBase_.markAsChimeric();
          ++chimeraCount;
          break;
        }
      }
      if (foundChi) {
        break;
      }
    }
    if (foundChi) {
      continue;
    }
    for (const auto &end : clus.endChiPos) {
      for (const auto &pos : iter::range(processedReads.size())) {
        if ((processedReads[pos].seqBase_.frac_ / clus.seqBase_.frac_) <
                parentFreqs ||
            pos == end.second) {
          continue;
        }
        alignerObj.alignVec(processedReads[pos], clus, local);
        alignerObj.profileAlignment(processedReads[pos], clus, 11, true, false,
                                    true, false, weighHomopolyer, 0,
                                    end.first + 1);
        if (chiOverlap.passErrorProfile(alignerObj.comp_)) {
          foundChi = true;
          clus.seqBase_.markAsChimeric();
          ++chimeraCount;
          break;
        }
      }
      if (foundChi) {
        break;
      }
    }
  }
}

void clusterCollapser::collapseTandems(std::vector<cluster> &processedReads,
                                       aligner &alignerObj, int runCutOff,
                                       int kLength, bool kMersByPosition,
                                       double freqCutoff, bool local,
                                       bool weighHomopolyer) {

  int counter = 0;
  auto kMaps = kmerCalculator::indexKmerMpas(processedReads, kLength, runCutOff,
                                             runCutOff);
  auto currentKMaps = alignerObj.getKmerMaps();
  alignerObj.setKmerMpas(kMaps);
  for (auto &clusIter : iter::reverse(processedReads)) {
    ++counter;
    if (counter % 20 == 0) {
      std::cout << counter << " out of " << processedReads.size() << std::endl;
    }
    for (auto &clusIterSecond : processedReads) {
      if (clusIter == clusIterSecond) {
        continue;
      }
      alignerObj.alignVec(clusIterSecond, clusIter, local);
      alignerObj.profileAlignment(clusIterSecond, clusIter, kLength,
                                  kMersByPosition, true, true, false,
                                  weighHomopolyer);
      // alignerObj.outPutParameterInfo(std::cout);
      if (alignerObj.checkForTandemRepeatGap() &&
          alignerObj.comp_.hqMismatches_ < 1 &&
          alignerObj.comp_.lqMismatches_ < 1 &&
          clusIterSecond.seqBase_.cnt_ >
              (freqCutoff * clusIter.seqBase_.cnt_)) {
        clusIter.remove = true;
        clusIterSecond.addRead(clusIter);
        break;
      }
    }
  }
  alignerObj.setKmerMpas(currentKMaps);
  return;
}
}  // namespace bib
