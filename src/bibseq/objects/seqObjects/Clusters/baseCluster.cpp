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

#include "baseCluster.hpp"
#include "bibseq/helpers/consensusHelper.hpp"
#include "bibseq/readVectorManipulation/readVectorHelpers.h"
#include "bibseq/IO/SeqIO/SeqOutput.hpp"

namespace bibseq {


baseCluster::baseCluster() : readObject() {
  firstReadCount_ = 0;
  needToCalculateConsensus_ = true;
}

baseCluster::baseCluster(const seqInfo& firstRead) : readObject(firstRead) {
  
	firstReadName_ = firstRead.name_;
  firstReadCount_ = firstRead.cnt_;
  reads_.emplace_back(std::make_shared<readObject>(firstRead));
  needToCalculateConsensus_ = true;
  remove = false;
  updateName();
}


void baseCluster::addRead(const baseCluster& otherCluster) {
  seqBase_.cnt_ += otherCluster.seqBase_.cnt_;
  seqBase_.frac_ = ((seqBase_.frac_ * reads_.size()) + otherCluster.seqBase_.frac_) /
                   (reads_.size() + 1);
  if (seqBase_.cnt_ / 2.0 >= firstReadCount_) {
    needToCalculateConsensus_ = true;
  }
  // needToCalculateConsensus = true;
  reads_.insert(reads_.end(), otherCluster.reads_.begin(), otherCluster.reads_.end());
}
void baseCluster::calculateConsensusToCurrent(aligner& alignerObj, bool setToConsensus){
	// if the cluster is only one read, no need to create consensus
	if (reads_.size() <= 1) {
		return;
	}
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus_) {
    return;
  }
  calculateConsensusTo(seqBase_, alignerObj, setToConsensus);

}




void baseCluster::calculateConsensusTo(const seqInfo & seqBase,
		aligner& alignerObj, bool setToConsensus) {
	// create the map for letter counters for each position
	std::map<uint32_t, charCounter> counters;
	// create a map in case of insertions
	std::map<uint32_t, std::map<uint32_t, charCounter>> insertions;
	std::map<int32_t, charCounter> beginningGap;
	auto getSeqBase =
			[](const std::shared_ptr<readObject> & read) ->const seqInfo& {return read->seqBase_;};

	consensusHelper::increaseCounters(seqBase_, reads_, getSeqBase, alignerObj, counters,
			insertions, beginningGap);
	calcConsensusInfo_ = seqBase_;
	for( auto & counter :counters){
		counter.second.resetAlphabet(true);
	}
	for( auto & counter :beginningGap){
		counter.second.resetAlphabet(true);
	}
	for( auto & counter :insertions){
		for( auto & subCounter : counter.second){
			subCounter.second.resetAlphabet(true);
		}
	}
	consensusHelper::genConsensusFromCounters(calcConsensusInfo_, counters, insertions,
			beginningGap);
	if (setToConsensus) {
		if (seqBase_.seq_ != calcConsensusInfo_.seq_) {
			seqBase_.seq_ = calcConsensusInfo_.seq_;
			setLetterCount();
		} else {
			seqBase_.seq_ = calcConsensusInfo_.seq_;
		}
		seqBase_.qual_ = calcConsensusInfo_.qual_;
	}
	previousErrorChecks_.clear();
	needToCalculateConsensus_ = false;
}

void baseCluster::calculateConsensus(aligner& alignerObj, bool setToConsensus) {
	// if the cluster is only one read, no need to create consensus
	if (reads_.size() <= 1) {
		return;
	}
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus_) {
    return;
  }
  auto averageSize = getAverageReadLength();
  seqInfo longestcluster;
  if (uAbsdiff(this->seqBase_.seq_.size(), averageSize) >
      0.1 * averageSize) {
    uint64_t biggest = 0;
    for (const auto& read : reads_) {
      if (read->seqBase_.seq_.length() > biggest) {
        longestcluster = read->seqBase_;
        biggest = read->seqBase_.seq_.length();
      }
    }
  } else {
    longestcluster = seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_);
  }

  calculateConsensusTo(longestcluster, alignerObj, setToConsensus);

}


//////align the current clusters to the curent consensus
std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
baseCluster::calculateAlignmentsToConsensus(aligner& alignObj) {
  std::vector<baseReadObject> withConAlignments;
  std::vector<baseReadObject> withOutConAlignments;

  for (const auto & read : reads_) {
    alignObj.alignCacheGlobal(*this, read);

    baseReadObject tempConsensus = alignObj.alignObjectA_;
    tempConsensus.seqBase_.name_ = seqBase_.name_;
    tempConsensus.seqBase_.name_.append("_consensus");
    baseReadObject tempRead = alignObj.alignObjectB_;
    tempRead.seqBase_.name_ = read->seqBase_.name_;

    withConAlignments.push_back(tempRead);
    withConAlignments.push_back(tempConsensus);
    withOutConAlignments.push_back(tempRead);
  }
  return {withConAlignments, withOutConAlignments};
}
// write out the alignments to the consensus
void baseCluster::writeOutAlignments(const std::string& directoryName,
                                     aligner& alignObj) {
  std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
      alignments = calculateAlignmentsToConsensus(alignObj);
  std::stringstream alignMentsFileName;
  std::stringstream alignmentsOnlyFileName;
  alignMentsFileName << directoryName << "align_" << seqBase_.name_ << ".fasta";
  alignmentsOnlyFileName << directoryName << "noConAlign_" << seqBase_.name_
                         << ".fasta";

  SeqIOOptions alignMentsFileNameOpts = SeqIOOptions::genFastaOut(alignMentsFileName.str());
  SeqOutput firstWriter(alignMentsFileNameOpts);
  firstWriter.openWrite(alignments.first);
  SeqIOOptions alignmentsOnlyFileNameOpts = SeqIOOptions::genFastaOut(alignmentsOnlyFileName.str());;
  SeqOutput secondWriter(alignmentsOnlyFileNameOpts);
  secondWriter.openWrite(alignments.second);
}

void baseCluster::writeClustersInDir(const std::string& workingDir,
		const SeqIOOptions & ioOptions) const {
	SeqIOOptions options(workingDir + seqBase_.name_,ioOptions.outFormat_, ioOptions.out_);
	writeClusters(options);
}

/// output the clusters currently clustered to the this cluster
void baseCluster::writeClusters(const SeqIOOptions & ioOptions) const {
	SeqOutput writer(ioOptions);
	writer.openWrite(reads_);
}

VecStr baseCluster::getReadNames() const { return readVec::getNames(reads_); }

double baseCluster::getAverageReadLength() const {
  if (reads_.size() == 1) {
    return (double)seqBase_.seq_.length();
  } else {

    int sumOfLength = 0;
    int numberOfReads = 0;
    double averageReadLength = 0.0;
    for (const auto& read : reads_) {
      sumOfLength += read->seqBase_.seq_.length() * read->seqBase_.cnt_;
      numberOfReads += read->seqBase_.cnt_;
    }
    averageReadLength = (double)sumOfLength / numberOfReads;
    return averageReadLength;
  }
}



readObject baseCluster::createRead() const {
	return readObject(seqBase_);
}

std::string toSlimJsonErrors(const comparison & comp){
	Json::Value ret;
	ret["highQualityMatches_"] = comp.highQualityMatches_;
	ret["hqMismatches_"] = comp.hqMismatches_;
	ret["largeBaseIndel_"] = comp.largeBaseIndel_;
	ret["lowKmerMismatches_"] = comp.lowKmerMismatches_;
	ret["lowQualityMatches_"] = comp.lowQualityMatches_;
	ret["lqMismatches_"] = comp.lqMismatches_;
	ret["oneBaseIndel_"] = comp.oneBaseIndel_;
	ret["twoBaseIndel_"] = comp.twoBaseIndel_;
	Json::FastWriter jWriter;
	return jWriter.write(ret);
}


bool baseCluster::compare(baseCluster & read, aligner & alignerObj,
		const IterPar & runParams, const CollapserOpts & collapserOptsObj) {
	bool ret = false;
	if (previousErrorChecks_.find(read.firstReadName_) != previousErrorChecks_.end()
			&& read.previousErrorChecks_.find(firstReadName_)
					!= read.previousErrorChecks_.end()) {
		ret = runParams.passErrorCheck(
				read.previousErrorChecks_.at(firstReadName_));
		//if(read.previousErrorChecks_.at(firstReadName_).distances_.eventBasedIdentity_ > .95 and read.seqBase_.cnt_ == 1 and !ret){
		//	std::cout << toSlimJsonErrors(read.previousErrorChecks_.at(firstReadName_)) << std::endl;
		//}
	} else {
    if(collapserOptsObj.alignOpts_.noAlign_){
    	alignerObj.noAlignSetAndScore(seqBase_, read.seqBase_);
    }else{
    	alignerObj.alignCacheGlobal(seqBase_, read.seqBase_);
    }
		comparison currentProfile = alignerObj.compareAlignment(seqBase_, read.seqBase_,
				 collapserOptsObj.kmerOpts_.checkKmers_);
		if(bib::containsSubString(read.seqBase_.name_, "lib9_Major_seq.002020_1_")){
			std::cout << currentProfile.toJson() << std::endl;
		}
		if (currentProfile.distances_.query_.coverage_ < 0.50
				|| currentProfile.distances_.ref_.coverage_< 0.50) {
			ret = false;
		} else {
			read.previousErrorChecks_[firstReadName_] = currentProfile;
			previousErrorChecks_[read.firstReadName_] = currentProfile;
			//std::stringstream ss;
			//ss << bib::json::toJson(currentProfile) << std::endl;;
			//std::cout << bib::removeAllWhitespace(ss.str()) << std::endl;;
			//ss.str("");
			//ss << bib::json::toJson(runParams.errors_) << std::endl;
			//std::cout << bib::removeAllWhitespace(ss.str()) << std::endl;;
			ret = runParams.passErrorCheck(currentProfile);
			//if (currentProfile.distances_.eventBasedIdentity_ > .95
			//		&& read.seqBase_.cnt_ == 1 and !ret) {
			//	std::cout << toSlimJsonErrors(currentProfile) << std::endl;
			//}
			//std::cout << ret << std::endl;
		}
	}
	return ret;
}

bool baseCluster::isClusterCompletelyChimeric() {
  for (const auto &read : reads_) {
    if (read->seqBase_.name_.find("CHI") == std::string::npos) {
      return false;
    }
  }
  return true;
}

bool baseCluster::isClusterAtLeastHalfChimeric() {
	size_t chiCount = 0;
	uint32_t chiReadCnt = 0;
	double total = readVec::getTotalReadCount(reads_);
	for (const auto &read : reads_) {
		if (read->seqBase_.name_.find("CHI") != std::string::npos) {
			++chiCount;
			chiReadCnt += read->seqBase_.cnt_;
		}
	}
	if (chiReadCnt >= total / 2.0) {
		return true;
	}
	return false;
}

bool baseCluster::isClusterAtLeastChimericCutOff(double cutOff) {
	size_t chiCount = 0;
	uint32_t chiReadCnt = 0;
	double total = readVec::getTotalReadCount(reads_);
	for (const auto &read : reads_) {
		if (read->seqBase_.name_.find("CHI") != std::string::npos) {
			++chiCount;
			chiReadCnt += read->seqBase_.cnt_;
		}
	}
	if (chiReadCnt/total >= cutOff) {
		return true;
	}
	return false;
}

}  // namespace bib
