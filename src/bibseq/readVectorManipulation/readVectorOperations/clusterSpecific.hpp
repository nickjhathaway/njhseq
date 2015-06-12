#pragma once
//
//  clusterSpecific.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 2/3/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//
//////////********************/////////////////////
// consider changing this to an object to alter vector of clusters or something,
// nick 7.17.13
//////////********************/////////////////////


#include "bibseq/objects/helperObjects/kmer.hpp"
#include "bibseq/utils/vectorUtils.hpp"
#include "bibseq/alignment/aligner.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations/massGetters.hpp"
//#include <bibcpp/files/fileUtils.hpp>

/// setting quality representation
namespace bibseq {
namespace clusterVec {
/// writters
template <typename T>
static void allWriteAllInputClusters(const std::vector<T> &reads,
                                     const std::string &workingDir) {
  for_each(reads,
           [&](const T &read) { read.writeOutInputClusters(workingDir); });
}
template <typename T>
static void allWriteOutClusters(const std::vector<T> &reads,
                                const std::string &workingDir, const readObjectIOOptions & ioOptions) {
  for_each(reads, [&](const T &read) { read.writeOutClusters(workingDir, ioOptions); });
}
/*
template <typename T>
static void allWriteOutAlignments(std::vector<T> &reads,
                                  const std::string &workingDirectory,
                                  aligner &alignObj) {
  std::string dir = bib::files::makeDir(workingDirectory, "alignments");
  for_each(reads, [&](T &read) { read.writeOutAlignments(dir, alignObj); });
}
template <typename T>
static void outPutAlignments(std::vector<T> &reads, aligner &alignerObj,
                             const std::string &directory) {
  allWriteOutAlignments(reads, directory, alignerObj);
}
template <typename T>
static void allWriteOutLongestAlignments(std::vector<T> &reads,
                                         const std::string &workingDirectory) {
  std::string dir = bib::files::makeDir(workingDirectory, "longestAlignments");
  for (auto &read : reads) {
    if (read.seqBase_.cnt_ > 1) {
      read.writeOutLongestAlignments(dir);
    }
  }
}*/
// manipulators
template <typename T>
static void allCalculateConsensus(std::vector<T> &reads, aligner &alignerObj,
                                  bool setToConsensus) {
  for_each(reads, [&](T &read) {
    read.calculateConsensus(alignerObj, setToConsensus);
  });
}

template <typename CLUSTER>
static void allSetFractionClusters(std::vector<CLUSTER> &reads) {
  size_t sizeOfRead = readVec::getTotalReadCount(reads);
  int count = 0;
  for (auto &read : reads) {
    read.setFractionByCount(sizeOfRead);
    ++count;
    for (auto &subSub : read.reads_) {
      subSub.setFractionByCount(sizeOfRead);
    }
  }
}

template <typename T>
static std::string returnSuperFromSub(const std::vector<T> &reads,
                                      const std::string &subName) {
  std::string notfound = "";
  for (const auto &read : reads) {
    for (const auto &sread : read.reads_) {
      if (subName == sread.name) {
        return read.name;
      }
    }
  }
  return notfound;
}

template <typename CLUSTER>
bool isClusterCompletelyChimeric(const std::vector<CLUSTER> &clusters) {
  for (const auto &cIter : clusters) {
    if (cIter.name.find("CHI") == std::string::npos) {
      return false;
    }
  }
  return true;
}

template<typename CLUSTER>
bool isClusterAtLeastHalfChimeric(const std::vector<CLUSTER> &clusters) {
	size_t chiCount = 0;
	uint32_t chiReadCnt = 0;
	double total = readVec::getTotalReadCount(clusters);
	for (const auto &clus : clusters) {
		if (clus.seqBase_.name_.find("CHI") != std::string::npos) {
			++chiCount;
			chiReadCnt += clus.seqBase_.cnt_;
		}
	}
	if (chiReadCnt >= total / 2.0) {
		return true;
	}
	return false;
}

template<typename CLUSTER>
bool isClusterAtLeastChimericCutOff(const std::vector<CLUSTER> &clusters, double cutOff) {
	size_t chiCount = 0;
	uint32_t chiReadCnt = 0;
	double total = readVec::getTotalReadCount(clusters);
	for (const auto &clus : clusters) {
		if (clus.seqBase_.name_.find("CHI") != std::string::npos) {
			++chiCount;
			chiReadCnt += clus.seqBase_.cnt_;
		}
	}
	if (chiReadCnt/total >= cutOff) {
		return true;
	}
	return false;
}

template <typename CLUSTER>
void markCompletelyChimericCluster(std::vector<CLUSTER> &clusters) {
  for (auto &clus : clusters) {
    if (isClusterCompletelyChimeric(clus.reads_)) {
      clus.markAsChimeric();
    }
  }
}
template <typename CLUSTER>
void markHalfChimericCluster(std::vector<CLUSTER> &clusters) {
  for (auto &clus : clusters) {
    if (isClusterCompletelyChimeric(clus.reads_)) {
      clus.markAsChimeric();
    }
  }
}

}  // namespace clusterVec
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "clusterSpecific.cpp"
#endif
