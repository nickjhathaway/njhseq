

#include "baseCluster.hpp"
#include "bibseq/helpers/consensusHelper.hpp"

namespace bibseq {


baseCluster::baseCluster() : readObject() {
  firstReadCount_ = 0;
  needToCalculateConsensus_ = true;
}

baseCluster::baseCluster(const readObject& firstRead) : readObject(firstRead.seqBase_) {
  
	firstReadName_ = firstRead.seqBase_.name_;
  firstReadCount_ = firstRead.seqBase_.cnt_;
  reads_.push_back(firstRead);
  needToCalculateConsensus_ = true;
  remove = false;
  updateName();
  //std::cout << "baseCluster constructor: " << std::endl;
  //std::cout << seqBase_.name_ << std::endl;
  //std::cout << seqBase_.cnt_ << std::endl;
  //std::cout << seqBase_.frac_ << std::endl;
}


void baseCluster::addRead(const readObject& newRead) {
  seqBase_.cnt_ += newRead.seqBase_.cnt_;
  seqBase_.frac_ = ((seqBase_.frac_ * reads_.size()) + newRead.seqBase_.frac_) /
                   (reads_.size() + 1);
  if (seqBase_.cnt_ / 2 > firstReadCount_) {
    needToCalculateConsensus_ = true;
  }
  // needToCalculateConsensus = true;
  reads_.push_back(newRead);
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
	std::map<uint32_t, charCounterArray> counters;
	// create a map in case of insertions
	std::map<uint32_t, std::map<uint32_t, charCounterArray>> insertions;
	std::map<int32_t, charCounterArray> beginningGap;
	auto getSeqBase =
			[](const readObject & read) ->const seqInfo& {return read.seqBase_;};

	consensusHelper::increaseCounters(seqBase_, reads_, getSeqBase, alignerObj, counters,
			insertions, beginningGap);
	calcConsensusInfo_ = seqBase_;
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
      if (read.seqBase_.seq_.length() > biggest) {
        longestcluster = read.seqBase_;
        biggest = read.seqBase_.seq_.length();
      }
    }
  } else {
    longestcluster =seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_);
  }

  calculateConsensusTo(longestcluster, alignerObj, setToConsensus);

}


//////align the current clusters to the curent consensus
std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
baseCluster::calculateAlignmentsToConsensus(aligner& alignObj) {
  std::vector<baseReadObject> withConAlignments;
  std::vector<baseReadObject> withOutConAlignments;

  for (std::vector<readObject>::iterator clusterIter = reads_.begin();
       clusterIter != reads_.end(); clusterIter++) {
    alignObj.alignVec(*this, *clusterIter, false);

    baseReadObject tempConsensus = alignObj.alignObjectA_;
    tempConsensus.seqBase_.name_ = seqBase_.name_;
    tempConsensus.seqBase_.name_.append("_consensus");
    baseReadObject tempRead = alignObj.alignObjectB_;
    tempRead.seqBase_.name_ = clusterIter->seqBase_.name_;

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

  readObjectIO reader;
  reader.writeSimpleVector(alignments.first, alignMentsFileName.str(),
                           "fastaQual", false,false, false);
  reader.writeSimpleVector(alignments.second, alignmentsOnlyFileName.str(),
                           "fastaQual", false, false, false);
}

/// output the clusters currently clustered to the this cluster
void baseCluster::writeOutClusters(const std::string& directoryName, const readObjectIOOptions & ioOptions) const {
  readObjectIO::write(reads_, directoryName + seqBase_.name_, ioOptions);
}

VecStr baseCluster::getReadNames() const { return readVec::getNames(reads_); }

double baseCluster::getAverageReadLength() const {
  if (reads_.size() == 1) {
    return (double)seqBase_.seq_.length();
  } else {

    int sumOfLength = 0;
    int numberOfReads = 0;
    double averageReadLength = 0.0;
    for (const auto& readsIter : reads_) {
      sumOfLength += readsIter.seqBase_.seq_.length() * readsIter.seqBase_.cnt_;
      numberOfReads += readsIter.seqBase_.cnt_;
    }
    averageReadLength = (double)sumOfLength / numberOfReads;
    return averageReadLength;
  }
}

void baseCluster::alignmentProfile(const std::string& workingDirectory,
                                   aligner& alignerObj, bool local,
                                   bool kmerChecking, int kLength,
                                   bool kmersByPosition,
                                   bool weighHomopolyers) const {
  alignmentProfiler::getInfoSingleComparison(
      reads_, *this, alignerObj, local,
      workingDirectory + seqBase_.name_ + "_profile.tab.txt", kLength,
      weighHomopolyers);
}

readObject baseCluster::createRead() const {
	return readObject(seqBase_);
}

void baseCluster::printDescription(std::ostream& out, bool deep) const {
  readObject::printDescription(out, deep);
  out << "baseCluster{" << std::endl << "firstReadName:" <<  firstReadName_
      << std::endl << "firstReadCount:" << firstReadCount_ << std::endl;
  if (deep) {
    out << "reads:" << std::endl;
    out << "std::vector<bib::readObject>" << std::endl;
    for (const auto& read : reads_) {
      read.printDescription(out, deep);
    }
    out << "previousErrorChecks" << std::endl;
    out << "std::map<std::string, errorProfile>" << std::endl;
    for (const auto& profile : previousErrorChecks_) {
      out << "\t" << profile.first << ":";
      profile.second.printDescription(out, deep);
    }
  }
  out << "longestAlignmentsQualities:" << std::endl;
  out << "std::vector<std::vector<uint32_t>>" << std::endl;

  out << "needToCalculateConsensus:" << needToCalculateConsensus_ << std::endl;
}

bool baseCluster::compare(baseCluster & read, aligner & alignerObj,
		const comparison & errorThreshold, const collapserOpts & collapserOptsObj) {
	bool ret = false;
	if (previousErrorChecks_.find(read.firstReadName_) != previousErrorChecks_.end()
			&& read.previousErrorChecks_.find(firstReadName_)
					!= read.previousErrorChecks_.end()) {
		ret = errorThreshold.passErrorProfile(
				read.previousErrorChecks_.at(firstReadName_));
	} else {
    if(collapserOptsObj.noAlign_){
    	alignerObj.noAlignSetAndScore(*this, read);
    }else{
    	alignerObj.alignVec(*this, read, collapserOptsObj.local_);
    }
		//alignerObj.alignVec(*this, read, collapserOptsObj.local_);
		alignerObj.profilePrimerAlignment(*this, read,
				collapserOptsObj.weighHomopolyer_);
		if (alignerObj.comp_.distances_.queryCoverage_ < 0.50) {
			ret = false;
		} else {
			comparison currentProfile = alignerObj.compareAlignment(*this, read,
					runningParameters(), collapserOptsObj.checkKmers_,
					collapserOptsObj.kmersByPosition_, collapserOptsObj.weighHomopolyer_);
			read.previousErrorChecks_[firstReadName_] = currentProfile;
			previousErrorChecks_[read.firstReadName_] = currentProfile;
			ret = errorThreshold.passErrorProfile(currentProfile);
		}
	}
	return ret;
}

bool baseCluster::compareId(baseCluster & read, aligner & alignerObj,
		const comparison & errorThreshold, const collapserOpts & collapserOptsObj) {
	bool ret = false;
	if (previousErrorChecks_.find(read.firstReadName_) != previousErrorChecks_.end()
			&& read.previousErrorChecks_.find(firstReadName_)
					!= read.previousErrorChecks_.end()) {
		ret = errorThreshold.passErrorProfile(
				read.previousErrorChecks_.at(firstReadName_));
	} else {
    if(collapserOptsObj.noAlign_){
    	alignerObj.noAlignSetAndScore(*this, read);
    }else{
    	alignerObj.alignVec(*this, read, collapserOptsObj.local_);
    }
		//alignerObj.alignVec(*this, read, collapserOptsObj.local_);
		alignerObj.profilePrimerAlignment(*this, read,
				collapserOptsObj.weighHomopolyer_);
		if (alignerObj.comp_.distances_.queryCoverage_ < 0.50) {
			ret = false;
		} else {
			comparison currentProfile = alignerObj.compareAlignment(*this, read,
					runningParameters(), collapserOptsObj.checkKmers_,
					collapserOptsObj.kmersByPosition_, collapserOptsObj.weighHomopolyer_);
			read.previousErrorChecks_[firstReadName_] = currentProfile;
			previousErrorChecks_[read.firstReadName_] = currentProfile;
			ret = errorThreshold.passIdThreshold(currentProfile);
		}
	}
	return ret;
}

}  // namespace bib
