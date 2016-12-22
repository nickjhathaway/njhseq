#pragma once
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
//
//  profiler.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/29/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils/utils.hpp"
#include "bibseq/readVectorManipulation.h"
#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/alignment/aligner.h"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/IO/SeqIO/SeqInput.hpp"
#include "bibseq/objects/seqObjects/Clusters/sampleCluster.hpp"
#include "bibseq/objects/helperObjects/hpRun.hpp"

namespace bibseq {

class profiler {

 public:


  // compare to reference sequences
	template<typename READ, typename REF>
	static VecStr compareToRefSingle(const std::vector<REF>& refSeqs,
			const READ& read, bool local, bool eventBased);

	template<typename READ, typename REF>
	static VecStr compareToRefSingle(const std::vector<REF>& refSeqs,
			const READ& read, aligner& alignerObj, bool local, bool eventBased);

	template<typename READ, typename REF>
	static void compareToRef(std::vector<REF> inputRefs,
			std::vector<READ> inputClusters, aligner& alignerObj,
			const std::string& workingDir, bool local, bool verbose, bool eventBased);


  // get info about read count breakdown.
  template <class CLUSTER>
  static void getMapInfo(const std::vector<CLUSTER>& reads,
                         const std::string& seqName,
                         const std::string& directory,
                         const std::string& filename, aligner& alignerObj);
  template <class CLUSTER>
  static void getMapInfo(const std::vector<CLUSTER>& reads,
                         const std::string& seqName,
                         const std::string& directory,
                         const std::string& filename,
                         const std::string& expected,
                         const std::string& barcode, aligner& alignerObj);
  template <class CLUSTER>
  static void getFractionInfo(const std::vector<CLUSTER>& reads,
                              const std::string& directory,
                              const std::string& filename);
  template <class CLUSTER>
  static void getFractionInfo(const std::vector<CLUSTER>& reads,
                              std::ostream & out);
  template <class CLUSTER>
  static void getFractionInfo(const std::vector<CLUSTER>& reads,
                              const std::string& directory,
                              const std::string& filename,
                              const std::string& refFilename, aligner& alignObj,
                              bool local);
  template <class CLUSTER>
  static void getFractionInfo(const std::vector<CLUSTER>& reads,
                              std::ostream & out,
                              std::vector<readObject> & refSeqs, aligner& alignObj,
                              bool local);
  template <class CLUSTER>
  static void getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                              const std::string& directory,
                              const std::string& filename);
  template <class CLUSTER>
  static void getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                              std::ostream & out);
  template <class CLUSTER>
  static void getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                              const std::string& directory,
                              const std::string& filename,
                              const std::string& refFilename, aligner& alignObj,
                              bool local);
  template <class CLUSTER>
  static void getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                              std::ostream & out,
                              std::vector<readObject> & refSeqs, aligner& alignObj,
                              bool local);

  template <class T>
  static void addHPRuns(std::vector<T>& reads,
                        std::multimap<uint32_t, hpRun>& hpRuns,
                        int sizeMinimum = 3);
  template <class T>
  static std::multimap<uint32_t, hpRun> findHPRuns(std::vector<T>& reads,
                                                   int sizeMinimum = 3);

  template <class T>
  static void addHPRunsSimple(std::vector<T>& reads,
                              std::map<std::string, std::vector<hpRun>>& hpRuns,
                              int sizeMinimum = 3);
  template <class T>
  static std::map<std::string, std::vector<hpRun>> findHPRunsSimple(
      std::vector<T>& reads, int sizeMinimum = 3);

  static void outputHpRunsSimpleMap(
      const std::map<std::string, std::vector<hpRun>>& runs, std::ostream& out);
  template <class T>
  static void addHPRunsWQ(std::vector<T>& reads,
                          std::multimap<uint32_t, hpRun>& hpRuns,
                          int sizeMinimum = 3);

  template <class T>
  static std::multimap<uint32_t, hpRun> findHPRunsWQ(std::vector<T>& reads,
                                                     int sizeMinimum = 3);

  template <class T>
  static std::map<uint32_t, hpRun> searchForHpRuns(const T& read);
  template <class T>
  static std::map<std::string, std::map<uint32_t, hpRun>>
      searchForHpRunsMultiple(const std::vector<T>& reads);
  template <class T>
  static void addBleed(std::vector<T>& reads,
                       std::multimap<uint32_t, hpRun>& hpRuns,
                       int sizeMinimum = 3);

  template <class T>
  static void outputHighestGCContent(std::vector<T>& reads, std::ostream& out);

  template <class T>
  static void outputFlowDataForGraphing(const std::vector<T>& reads,
                                        std::ostream& out);

  template <class T>
  static void outputQualityForGraphing(const std::vector<T>& reads,
                                       std::ostream& out);
  template <class T>
  static void allOutputIdenticalReadComp(const std::vector<T>& reads,
                                         const std::string& workingDir);

  template <class T>
  static void outputGraphingQuality(std::vector<T>& reads);
};

template <typename READ, typename REF>
VecStr profiler::compareToRefSingle(const std::vector<REF>& refSeqs,
		const READ& read, aligner& alignerObj, bool local, bool eventBased) {

	double bestScore = std::numeric_limits<double>::min();
	std::vector<uint32_t> bestReadPos;

	for (const auto inputPos : iter::range(refSeqs.size())) {
		const auto & input = refSeqs[inputPos];
		if (getSeqBase(read).name_ == input.seqBase_.name_) {
			continue;
		}
		alignerObj.alignCache(input, read, local);
		double currentScore = 0;
		if (eventBased) {
			currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
		} else {
			currentScore = alignerObj.parts_.score_;
		}
		if (currentScore == bestScore) {
			bestReadPos.emplace_back(inputPos);
		} else if (alignerObj.parts_.score_ > bestScore) {
			bestReadPos.clear();
			bestReadPos.emplace_back(inputPos);
			bestScore = alignerObj.parts_.score_;
		}
	}
  VecStr ret;
  for (const auto& bestPos : bestReadPos) {
  	const auto & best = refSeqs[bestPos];
		alignerObj.alignCache(best, read, local);
		alignerObj.profilePrimerAlignment(best, read);
    std::stringstream profile;
    double score = 0;
    if(eventBased){
    	score = alignerObj.comp_.distances_.eventBasedIdentity_;
    } else {
    	score = alignerObj.parts_.score_;
    }
    profile << best.seqBase_.name_
    		<< "," << score
				<< "," << alignerObj.comp_.oneBaseIndel_
				<< "," << alignerObj.comp_.twoBaseIndel_
				<< "," << alignerObj.comp_.largeBaseIndel_
				<< "," << alignerObj.comp_.lqMismatches_
				<< "," << alignerObj.comp_.hqMismatches_;
    ret.emplace_back(profile.str());
  }
  return ret;
}

template <typename READ, typename REF>
VecStr profiler::compareToRefSingle(const std::vector<REF>& refSeqs,
		const READ& read, bool local, bool eventBased) {
  // set the gap scoring for comparisons so alignment is semi-global
	uint64_t maxReadLen = 0;
	readVec::getMaxLength(read, maxReadLen);
	readVec::getMaxLength(refSeqs, maxReadLen);
	aligner alignerObj(maxReadLen,
			gapScoringParameters(5,1,0,0,0,0),
			substituteMatrix(2,-2));
	alignerObj.countEndGaps_ = false;
	return compareToRefSingle(refSeqs, read, alignerObj, local, eventBased);
}

template <typename READ, typename REF>
void profiler::compareToRef(std::vector<REF> inputRefs,
		std::vector<READ> inputClusters, aligner& alignerObj,
		const std::string& workingDir, bool local, bool verbose, bool eventBased) {
  std::ofstream profileInfoFile;
  openTextFile(profileInfoFile, workingDir + "refComparisonInfo.tab.txt",
               ".txt", true, false);
  std::ofstream tempFile;
  openTextFile(tempFile, workingDir + "tempFilealns.fasta", ".fasta", true,
               false);
  uint32_t counter = 0;
  profileInfoFile
      << "ReadNumber\tReadId\tReadFraction\tBestRef\tscore\thqScore\t1bIndel\t2bI"
         "ndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;

  readVec::lowerCaseBasesToUpperCase(inputClusters);
  readVec::lowerCaseBasesToUpperCase(inputRefs);
  for (const auto& input : inputClusters) {
    if ((counter + 1) % 5 == 0 && verbose) {
      std::cout << "Currently on read " << counter + 1 << " out of "
                << inputClusters.size() << std::endl;
    }
    double bestScore = std::numeric_limits<double>::min();
    std::vector<uint64_t> bestRead;
    for (const auto& refPos : iter::range(inputRefs.size())) {
    	const auto & ref = inputRefs[refPos];
      if (input.seqBase_.name_ == ref.seqBase_.name_) {
        continue;
      }
      alignerObj.alignCache(ref, input, local);
      double currentScore = 0;
      if(eventBased){
        alignerObj.profileAlignment(ref, input, false, true, false);
      	currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      }else{
      	currentScore = alignerObj.parts_.score_;
      }

      if (currentScore == bestScore) {
        bestRead.push_back(refPos);
      }
      if (currentScore > bestScore) {
        bestRead.clear();
        bestRead.push_back(refPos);
        bestScore = currentScore;
      }
    }
    for (const auto& bestPos : bestRead) {
    	const auto & best = inputRefs[bestPos];
      alignerObj.alignCache(best, input, local);
      tempFile << ">" << best.seqBase_.name_ << std::endl;
      tempFile << alignerObj.alignObjectA_.seqBase_.seq_ << std::endl;
      tempFile << ">" << input.seqBase_.name_ << std::endl;
      tempFile << alignerObj.alignObjectB_.seqBase_.seq_ << std::endl;
      double score = 0;
      alignerObj.profileAlignment(best, input, false, true, false);
      if(eventBased){
      	score = alignerObj.comp_.distances_.eventBasedIdentity_;
      } else {
      	score = alignerObj.parts_.score_;
      }


      profileInfoFile << counter << "\t" << input.seqBase_.name_ << "\t"
					<< input.seqBase_.frac_ << "\t" << best.seqBase_.name_
					<< "\t" << score
					<< "\t" << alignerObj.comp_.distances_.eventBasedIdentity_
					<< "\t" << alignerObj.comp_.oneBaseIndel_ << "\t"
					<< alignerObj.comp_.twoBaseIndel_ << "\t"
					<< alignerObj.comp_.largeBaseIndel_ << "\t"
					<< alignerObj.comp_.lqMismatches_ << "\t"
					<< alignerObj.comp_.hqMismatches_ << std::endl;
    }
    counter++;
  }
}



template <class T>
void profiler::addHPRuns(std::vector<T>& reads,
                         std::multimap<uint32_t, hpRun>& hpRuns,
                         int sizeMinimum) {
  std::multimap<uint32_t, hpRun>::iterator hpRunsIter;
  for (typename std::vector<T>::const_iterator iter = reads.begin();
       iter != reads.end(); ++iter) {
    int currentCount = 1;
    size_t pos = 0;
    for (uint32_t i = 1; i < iter->seqBase_.seq_.length(); i++) {
      if (iter->seqBase_.seq_[i] == iter->seqBase_.seq_[i - 1]) {
        currentCount++;
      } else {
        if (currentCount < sizeMinimum) {

        } else {
          bool foundMatch = false;
          for (hpRunsIter = hpRuns.begin(); hpRunsIter != hpRuns.end();
               ++hpRunsIter) {
            if (hpRunsIter->second.pos == (int)pos &&
                hpRunsIter->second.base == iter->seqBase_.seq_[i - 1] &&
                hpRunsIter->second.runSize == currentCount) {
              hpRunsIter->second += iter->seqBase_.cnt_;
              foundMatch = true;
              break;
            }
          }
          if (!foundMatch) {
            hpRuns.insert(
                std::make_pair(pos, hpRun(iter->seqBase_.seq_[i - 1], (int)pos,
                                          currentCount, iter->seqBase_.cnt_)));
          }
        }
        currentCount = 1;
        pos = i;
      }
    }
  }
}
template <class T>
std::multimap<uint32_t, hpRun> profiler::findHPRuns(std::vector<T>& reads,
                                                    int sizeMinimum) {
  std::multimap<uint32_t, hpRun> hpRuns;
  addHPRuns(reads, hpRuns);
  return hpRuns;
};

template <class T>
void profiler::addHPRunsSimple(
    std::vector<T>& reads, std::map<std::string, std::vector<hpRun>>& hpRuns,
    int sizeMinimum) {
  std::map<std::string, std::vector<hpRun>>::iterator hpRunsIter;
  std::vector<hpRun>::iterator hpVecIter;
  for (typename std::vector<T>::const_iterator iter = reads.begin();
       iter != reads.end(); ++iter) {
    int currentCount = 1;
    size_t pos = 0;
    for (uint32_t i = 1; i < iter->seqBase_.seq_.length(); i++) {
      if (iter->seqBase_.seq_[i] == iter->seqBase_.seq_[i - 1]) {
        currentCount++;
      } else {
        if (currentCount < sizeMinimum) {

        } else {
          bool foundMatch = false;
          if (hpRuns.find(iter->seqBase_.seq_.substr(i - 1, 1)) ==
              hpRuns.end()) {
            std::vector<hpRun> tempVec;
            tempVec.push_back(hpRun(iter->seqBase_.seq_[i - 1], (int)pos,
                                    currentCount, iter->seqBase_.cnt_));
            hpRuns.insert(
                std::make_pair(iter->seqBase_.seq_.substr(i - 1, 1), tempVec));
          } else {
            for (hpVecIter =
                     hpRuns[iter->seqBase_.seq_.substr(i - 1, 1)].begin();
                 hpVecIter !=
                     hpRuns[iter->seqBase_.seq_.substr(i - 1, 1)].end();
                 ++hpVecIter) {
              if (hpVecIter->pos == (int)pos &&
                  hpVecIter->base == iter->seqBase_.seq_[i - 1] &&
                  hpVecIter->runSize == currentCount) {
                (*hpVecIter) += iter->seqBase_.cnt_;
                foundMatch = true;
                break;
              }
            }
            if (!foundMatch) {
              hpRuns[iter->seqBase_.seq_.substr(i - 1, 1)]
                  .push_back(hpRun(iter->seqBase_.seq_[i - 1], (int)pos,
                                   currentCount, iter->seqBase_.cnt_));
            }
          }
        }
        currentCount = 1;
        pos = i;
      }
    }
  }
}
template <class T>
std::map<std::string, std::vector<hpRun>> profiler::findHPRunsSimple(
    std::vector<T>& reads, int sizeMinimum) {
  std::map<std::string, std::vector<hpRun>> hpRuns;
  addHPRunsSimple(reads, hpRuns);
  return hpRuns;
};

template <class T>
void profiler::addHPRunsWQ(std::vector<T>& reads,
                           std::multimap<uint32_t, hpRun>& hpRuns,
                           int sizeMinimum) {
  std::multimap<uint32_t, hpRun>::iterator hpRunsIter;
  for (typename std::vector<T>::const_iterator iter = reads.begin();
       iter != reads.end(); ++iter) {
    int currentCount = 1;
    size_t pos = 0;
    std::vector<int> currentQuals;
    for (uint32_t i = 1; i < iter->seqBase_.seq_.length(); i++) {
      if (iter->seqBase_.seq_[i] == iter->seqBase_.seq_[i - 1]) {
        if (currentCount == 1) {
          currentQuals.push_back(iter->qual[i - 1]);
          currentQuals.push_back(iter->qual[i]);
        } else {
          currentQuals.push_back(iter->qual[i]);
        }
        currentCount++;

      } else {
        if (currentCount < sizeMinimum) {

        } else {
          bool foundMatch = false;
          for (hpRunsIter = hpRuns.begin(); hpRunsIter != hpRuns.end();
               ++hpRunsIter) {
            if (hpRunsIter->second.pos == (int)pos &&
                hpRunsIter->second.base == iter->seqBase_.seq_[i - 1] &&
                hpRunsIter->second.runSize == currentCount) {
              hpRunsIter->second += iter->seqBase_.cnt_;
              hpRunsIter->second.increaseCount(0, vectorMean(currentQuals));
              foundMatch = true;
              break;
            }
          }
          if (!foundMatch) {
            std::vector<int> tempQual;
            tempQual.push_back(vectorMean(currentQuals));
            hpRuns.insert(std::make_pair(
                pos, hpRun(iter->seqBase_.seq_[i - 1], (int)pos, currentCount,
                           iter->seqBase_.cnt_, tempQual)));
          }
        }
        currentQuals.clear();
        currentCount = 1;
        pos = i;
      }
    }
  }
}
template <class T>
std::multimap<uint32_t, hpRun> profiler::findHPRunsWQ(std::vector<T>& reads,
                                                      int sizeMinimum) {
  std::multimap<uint32_t, hpRun> hpRuns;
  addHPRunsWQ(reads, hpRuns);
  return hpRuns;
};
template <class T>
void profiler::addBleed(std::vector<T>& reads,
                        std::multimap<uint32_t, hpRun>& hpRuns,
                        int sizeMinimum) {
  std::multimap<uint32_t, hpRun>::iterator hpRunsIter;
  for (typename std::vector<T>::const_iterator iter = reads.begin();
       iter != reads.end(); ++iter) {
    int currentCount = 1;
    size_t pos = 0;
    std::vector<int> currentQuals;
    for (uint32_t i = 1; i < iter->seqBase_.seq_.length(); i++) {
      if (iter->seqBase_.seq_[i] == iter->seqBase_.seq_[i - 1]) {
        if (currentCount == 1) {
          currentQuals.push_back(iter->qual[i - 1]);
          currentQuals.push_back(iter->qual[i]);
        } else {
          currentQuals.push_back(iter->qual[i]);
        }
        currentCount++;

      } else {
        if (currentCount < sizeMinimum) {

        } else {
          bool foundMatch = false;
          for (hpRunsIter = hpRuns.begin(); hpRunsIter != hpRuns.end();
               ++hpRunsIter) {
            if (hpRunsIter->second.pos == (int)pos &&
                hpRunsIter->second.base == iter->seqBase_.seq_[i - 1] &&
                hpRunsIter->second.runSize == currentCount) {
              hpRunsIter->second += iter->seqBase_.cnt_;
              hpRunsIter->second.increaseCount(0, vectorMean(currentQuals));
              foundMatch = true;
              break;
            }
          }
          if (!foundMatch) {
            std::vector<int> tempQual;
            tempQual.push_back(vectorMean(currentQuals));
            hpRuns.insert(std::make_pair(
                pos, hpRun(iter->seqBase_.seq_[i - 1], (int)pos, currentCount,
                           iter->seqBase_.cnt_, tempQual)));
          }
        }
        currentQuals.clear();
        currentCount = 1;
        pos = i;
      }
    }
  }
}
template <class T>
std::map<uint32_t, hpRun> profiler::searchForHpRuns(const T& read) {
  std::map<uint32_t, hpRun> ans;
  int currentCount = 1;
  size_t pos = 0;
  std::vector<int> currentQuals;
  int numOfHPRuns = 1;
  for (uint32_t i = 1; i < read.seqBase_.seq_.length(); i++) {
    if (read.seqBase_.seq_[i] == read.seqBase_.seq_[i - 1]) {
      if (currentCount == 1) {
        currentQuals.push_back(read.seqBase_.qual_[i - 1]);
        currentQuals.push_back(read.seqBase_.qual_[i]);
      } else {
        currentQuals.push_back(read.seqBase_.qual_[i]);
      }
      currentCount++;

    } else {
      ++numOfHPRuns;
      if (currentCount > 1) {
        ans.insert(std::make_pair(
            pos, hpRun(read.seqBase_.seq_[i - 1], (int)pos, numOfHPRuns,
                       currentCount, read.seqBase_.cnt_, currentQuals)));
      }
      currentQuals.clear();
      currentCount = 1;
      pos = i;
    }
  }
  return ans;
}

template <class T>
std::map<std::string, std::map<uint32_t, hpRun>>
profiler::searchForHpRunsMultiple(const std::vector<T>& reads) {
  std::map<std::string, std::map<uint32_t, hpRun>> ans;
  for (typename std::vector<T>::const_iterator rIter = reads.begin();
       rIter != reads.end(); ++rIter) {
    ans.insert(std::make_pair(rIter->seqBase_.name_, searchForHpRuns(*rIter)));
  }
  return ans;
}

template <class T>
void profiler::outputHighestGCContent(std::vector<T>& reads,
                                      std::ostream& out) {
  readObject highestGCContent;
  for (typename std::vector<T>::iterator iter = reads.begin();
       iter != reads.end(); ++iter) {
    iter->setLetterCount();
    iter->counter.setGcContent();
    if (iter == reads.begin()) {
      highestGCContent = *iter;
    } else {
      if (iter->counter.gcContent > highestGCContent.counter_.gcContent_) {
        highestGCContent = *iter;
      }
    }
  }
  // out<<">"<<highestGCContent.seqBase_.name_<<std::endl;
  out << highestGCContent.seqBase_.name_ << std::endl;
  out << highestGCContent.counter_.gcContent_ * 100.00 << "%" << std::endl;
}

template <class T>
void profiler::outputFlowDataForGraphing(const std::vector<T>& reads,
                                         std::ostream& out) {
  typename std::vector<T>::const_iterator rIter;
  size_t max = 0;
  for (rIter = reads.begin(); rIter != reads.end(); ++rIter) {
    if (rIter == reads.begin()) {
      out << rIter->name;
    } else {
      out << "\t" << rIter->name;
    }
    if (rIter->flowValues.size() > max) {
      max = rIter->flowValues.size();
    }
  }
  out << std::endl;
  for (size_t i = 0; i < max; ++i) {
    for (rIter = reads.begin(); rIter != reads.end(); ++rIter) {
      if (rIter == reads.begin()) {
        if (i >= rIter->flowValues.size()) {
          out << "";
        } else {
          out << rIter->flowValues[i];
        }
      } else {
        if (i >= rIter->flowValues.size()) {
          out << "";
        } else {
          out << "\t" << rIter->flowValues[i];
        }
      }
    }
    out << std::endl;
  }
  return;
}
template <class T>
void profiler::outputQualityForGraphing(const std::vector<T>& reads,
                                        std::ostream& out) {
  typename std::vector<T>::const_iterator rIter;
  size_t max = 0;
  for (rIter = reads.begin(); rIter != reads.end(); ++rIter) {
    if (rIter == reads.begin()) {
      out << rIter->seqBase_.name_;
    } else {
      out << "\t" << rIter->seqBase_.name_;
    }
    if (rIter->seqBase_.qual_.size() > max) {
      max = rIter->seqBase_.qual_.size();
    }
  }
  out << std::endl;
  for (size_t i = 0; i < max; ++i) {
    for (rIter = reads.begin(); rIter != reads.end(); ++rIter) {
      if (rIter == reads.begin()) {
        if (i >= rIter->seqBase_.qual_.size()) {
          out << "";
        } else {
          out << rIter->seqBase_.qual_[i];
        }
      } else {
        if (i >= rIter->seqBase_.qual_.size()) {
          out << "";
        } else {
          out << "\t" << rIter->seqBase_.qual_[i];
        }
      }
    }
    out << std::endl;
  }
  return;
}

template <class T>
void profiler::allOutputIdenticalReadComp(const std::vector<T>& reads,
                                          const std::string& workingDir) {
  typename std::vector<T>::const_iterator rIter;
  for (rIter = reads.begin(); rIter != reads.end(); ++rIter) {
    rIter->outputInfoComp(workingDir);
  }

  return;
}
// this doesn't look right
template <class T>
void profiler::outputGraphingQuality(std::vector<T>& reads) {
  typename std::vector<T>::iterator seqIt;
  for (seqIt = reads.begin(); seqIt != reads.end(); seqIt++) {
    StringToUpper(seqIt->seqBase_.seq_);
  }
}
template <class CLUSTER>
void profiler::getFractionInfo(const std::vector<CLUSTER>& reads,
                               const std::string& directory,
                               const std::string& filename,
                               const std::string& refFilename,
                               aligner& alignObj, bool local) {
  std::ofstream outFile;
  openTextFile(outFile, directory + filename, ".tab.txt", false, false);
  SeqIOOptions opts = SeqIOOptions::genFastaIn(refFilename);
  SeqInput reader(opts);
  reader.openIn();
  auto refReads = reader.readAllReads<readObject>();

  getFractionInfo(reads, outFile, refReads, alignObj, local);
}

template <class CLUSTER>
void profiler::getFractionInfo(const std::vector<CLUSTER>& reads,
                               const std::string& directory,
                               const std::string& filename) {
  std::ofstream outFile;
  openTextFile(outFile, directory + filename, ".tab.txt", false, false);
  getFractionInfo(reads, outFile);
}
template <class CLUSTER>
void profiler::getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                               const std::string& directory,
                               const std::string& filename,
                               const std::string& refFilename,
                               aligner& alignObj, bool local) {
  std::ofstream outFile;
  openTextFile(outFile, directory + filename, ".tab.txt", false, false);
  SeqIOOptions opts = SeqIOOptions::genFastaIn(refFilename);
  SeqInput reader(opts);
  reader.openIn();
  auto refReads = reader.readAllReads<readObject>();

  getFractionInfoCluster(reads, outFile, refReads, alignObj, local);
}

template <class CLUSTER>
void profiler::getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                               const std::string& directory,
                               const std::string& filename) {
  std::ofstream outFile;
  openTextFile(outFile, directory + filename, ".tab.txt", false, false);
  getFractionInfoCluster(reads, outFile);
}

template <class CLUSTER>
void profiler::getFractionInfo(const std::vector<CLUSTER>& reads,
                               std::ostream & out) {
  out << "ClusterNumber\tClusterId\tClusterSize\tclusterFraction" << std::endl;
  int currentCluster = 0;
  int totalReads = readVec::getTotalReadCount(reads);

  for (const auto& clus : reads) {
    out << currentCluster << "\t" << clus.seqBase_.name_ << "\t"
            << clus.seqBase_.cnt_ << "\t"
            << (double)clus.seqBase_.cnt_ / totalReads<< std::endl;
    currentCluster++;
  }
}
template <class CLUSTER>
void profiler::getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                               	 	 	 	 std::ostream & out) {
  out << "ClusterNumber\tClusterId\tClusterSize\tNumOfClusters\tclusterFraction" << std::endl;
  int currentCluster = 0;
  int totalReads = readVec::getTotalReadCount(reads);

  for (const auto& clus : reads) {
    out << currentCluster << "\t" << clus.seqBase_.name_ << "\t"
            << clus.seqBase_.cnt_ << "\t" << clus.reads_.size() << "\t"
            << (double)clus.seqBase_.cnt_ / totalReads<< std::endl;
    currentCluster++;
  }
}
template <class CLUSTER>
void profiler::getFractionInfo(const std::vector<CLUSTER>& reads,
                            std::ostream & out,
                            std::vector<readObject> & refSeqs, aligner& alignObj,
                            bool local){
  out << "ClusterNumber\tClusterId\tClusterSize\tclusterFraction\tB"
             "estRef\t1bIndel\t2bIndel\t>2bIndel\tmismatch" << std::endl;
  int currentCluster = 0;
  int totalReads = readVec::getTotalReadCount(reads);

  for (const auto& clusIter : reads) {
    VecStr expectsComp = compareToRefSingle(refSeqs, clusIter, alignObj,
                                            local, false);
    out << currentCluster << "\t" << clusIter.seqBase_.name_ << "\t"
            << clusIter.seqBase_.cnt_ << "\t"
            << (double)clusIter.seqBase_.cnt_ / totalReads;
    for (const auto& expect : expectsComp) {
      VecStr toks = tokenizeString(expect, ",");
      out << "\t" << vectorToString(toks, "\t");
    }
    out << std::endl;
    currentCluster++;
  }
}
template <class CLUSTER>
void profiler::getFractionInfoCluster(const std::vector<CLUSTER>& reads,
                            std::ostream & out,
                            std::vector<readObject> & refSeqs, aligner& alignObj,
                            bool local){
  out << "ClusterNumber\tClusterId\tClusterSize\tNumOfClusters\tclusterFraction\tB"
             "estRef\t1bIndel\t2bIndel\t>2bIndel\tmismatch" << std::endl;
  int currentCluster = 0;
  int totalReads = readVec::getTotalReadCount(reads);

  for (const auto& clusIter : reads) {
    VecStr expectsComp = compareToRefSingle(refSeqs, clusIter, alignObj,
                                            local, false);
    out << currentCluster << "\t" << clusIter.seqBase_.name_ << "\t"
            << clusIter.seqBase_.cnt_ << "\t" << clusIter.reads_.size() << "\t"
            << (double)clusIter.seqBase_.cnt_ / totalReads;
    for (const auto& expect : expectsComp) {
      VecStr toks = tokenizeString(expect, ",");
      out << "\t" << vectorToString(toks, "\t");
    }
    out << std::endl;
    currentCluster++;
  }
}
template <class CLUSTER>
void profiler::getMapInfo(const std::vector<CLUSTER>& reads,
                          const std::string& seqName,
                          const std::string& directory,
                          const std::string& filename, aligner& alignerObj) {
  std::ofstream outFile;
  openTextFile(outFile, directory + filename, ".tab.txt", false, false);
  outFile
      << "seqName\trefId\treadsMapped\tmappedFraction\tmeanErrors\tmedianErros"
      << std::endl;
  int currentCluster = 0;
  //int totalReads = readVec::getTotalReadCount(reads);
  double totalReads = std::accumulate(reads.begin(), reads.end(),
  		0, [&](double res, const CLUSTER & read){return res+= read.seqBase_.cnt_;});
  for (const auto& clusIter : reads) {
    std::vector<double> mismatches;
    for (const auto& read : clusIter.reads_) {
      if (read.seqBase_.name_ == clusIter.seqBase_.name_) {
        continue;
      }
      alignerObj.alignCache(clusIter.seqBase_, read.seqBase_, false);
      alignerObj.profileAlignment(clusIter.seqBase_, read.seqBase_, false, true, false);
      for (int i = 0; i < std::ceil(read.seqBase_.cnt_); ++i) {
        mismatches.emplace_back(alignerObj.comp_.hqMismatches_);
      }
      /*
      if ((int) read.seqBase_.cnt_ < read.seqBase_.cnt_){
        mismatches.emplace_back(alignerObj.errors_.hqMismatches_ *
      (read.seqBase_.cnt_ - (int) read.seqBase_.cnt_ ) );
      }*/
    }

    outFile << seqName << "\t" << clusIter.seqBase_.name_ << "\t"
            << clusIter.seqBase_.cnt_ << "\t"
            << (double)clusIter.seqBase_.cnt_ / totalReads << "\t"
            << vectorMean(mismatches) << "\t" << vectorMedianRef(mismatches)
            << std::endl;
    currentCluster++;
  }
}
template <class CLUSTER>
void profiler::getMapInfo(const std::vector<CLUSTER>& reads,
                          const std::string& seqName,
                          const std::string& directory,
                          const std::string& filename,
                          const std::string& expected,
                          const std::string& barcode, aligner& alignerObj) {
  std::ofstream outFile;
  openTextFile(outFile, directory + filename, ".tab.txt", false, false);
  outFile << "seqName\tbarcode\texpected\trefId\treadsMapped\tmappedFraction\tm"
             "eanErrors\tmedianErros" << std::endl;
  int currentCluster = 0;
  //int totalReads = readVec::getTotalReadCount(reads);
  double totalReads = std::accumulate(reads.begin(), reads.end(),
    		0, [&](double res, const CLUSTER & read){return res+= read.seqBase_.cnt_;});
  for (const auto& clusIter : reads) {
    std::vector<double> mismatches;
    for (const auto& read : clusIter.reads_) {
      if (read.seqBase_.name_ == clusIter.seqBase_.name_) {
        continue;
      }
      alignerObj.alignCache(clusIter.seqBase_, read.seqBase_, false);
      alignerObj.profileAlignment(clusIter.seqBase_, read.seqBase_, false, true, false);
      for (int i = 0; i < std::ceil(read.seqBase_.cnt_); ++i) {
        mismatches.emplace_back(alignerObj.comp_.hqMismatches_);
      }
    }
    outFile << seqName << "\t" << barcode << "\t" << expected << "\t"
            << clusIter.seqBase_.name_ << "\t" << clusIter.seqBase_.cnt_ << "\t"
            << (double)clusIter.seqBase_.cnt_ / totalReads << "\t"
            << vectorMedianRef(mismatches) << "\t" << vectorMedianRef(mismatches)
            << std::endl;
    currentCluster++;
  }
}

}  // namespace bibseq


