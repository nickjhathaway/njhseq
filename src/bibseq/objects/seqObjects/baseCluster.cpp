//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

namespace bibseq {

void baseCluster::addRead(const readObject& newRead) {
  seqBase_.cnt_ += newRead.seqBase_.cnt_;
  seqBase_.frac_ = ((seqBase_.frac_ * reads_.size()) + newRead.seqBase_.frac_) /
                   (reads_.size() + 1);
  cumulativeFraction += newRead.cumulativeFraction;
  if (seqBase_.cnt_ / 2 > firstReadCount) {
    needToCalculateConsensus = true;
  }
  // needToCalculateConsensus = true;
  reads_.push_back(newRead);
}
void baseCluster::calculateConsensusNew(aligner& alignerObj, bool setToConsensus) {
  // if the cluster is only one read, no need to create consensus
  if (seqBase_.cnt_ == 1 || reads_.size() <= 1) {
    return;
  }
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus) {
    return;
  }
  // find the longest cluster
  // if the the current consensus is much shorter than the average amount of reads than choose the longest
  // read, this is to prevent consensus being built with only the a small portion of the sequences in the
  // the cluster
  auto averageSize = getAverageReadLength();
  readObject longestcluster;
  if (true) {
  //if (uAbsdiff(this->seqBase_.seq_.size(),averageSize) >
    //  0.1 * averageSize) {
    uint64_t biggest = 0;
    for (const auto& read : reads_) {
      if (read.seqBase_.seq_.length() > biggest) {
        longestcluster = read;
        biggest = read.seqBase_.seq_.length();
      }
    }
  } else {
    longestcluster =
        readObject(seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_));
  }


  // clear the vectors in case of previous consensus calculations
  longestAlignments.clear();
  longestAlignmentsQualities.clear();
  longestAlingmentsRef.clear();

  calculatedConsensus.clear();
  calculatedConsensusQuality.clear();

  for (const auto & read : reads_) {
    alignerObj.alignVec(longestcluster, read, false);
    longestAlignments.push_back(alignerObj.alignObjectB_.seqBase_.seq_);
    longestAlignmentsQualities.push_back(
        alignerObj.alignObjectB_.seqBase_.qual_);
    longestAlingmentsRef.push_back(alignerObj.alignObjectA_.seqBase_.seq_);
  }

  // create the map for letter counters for each position
  std::map<uint32_t, charCounterArray> counter;
  // create a map in case of insertions
  std::map<uint32_t, std::map<uint32_t, charCounterArray>> insertions;
  // counter of the reads
  //int numberOfReads = 0;
  // now iterator over the various clusters and alignments to count the bases
  // and their qualities
  // VecStr
  // stringsWithoutInsertions(longestAlingmentsRef.size(),"");
  int count = -1;
  for (const auto & readPos : iter::range(len(reads_)) ) {
    ++count;
    // std::cout<<count+1<<"/"<<longestAlingmentsRef.size()<<std::endl;
    // the offset for the insertions
    uint32_t offSet = 0;
    // get the size of the cluster being processed
    //int sizeAdjustment = (int)clusterIter->seqBase_.cnt_;
    //numberOfReads += (int)clusterIter->seqBase_.cnt_;
    int currentOffset = 1;
    for (uint32_t i = 0; i < len(longestAlingmentsRef[readPos]); i++) {
      // if the longest reference has an insertion in it put it in the
      // insertions letter counter map
      if (longestAlingmentsRef[readPos][i] == '-') {
      	insertions[i - offSet][currentOffset].increaseCountOfBaseQual(longestAlignments[readPos][i],
      			longestAlignmentsQualities[readPos][i], reads_[readPos].seqBase_.cnt_ );

        currentOffset++;
        offSet++;
        continue;
      }
      currentOffset = 1;
      counter[i - offSet].increaseCountOfBaseQual(longestAlignments[readPos][i],
    			longestAlignmentsQualities[readPos][i],reads_[readPos].seqBase_.cnt_ );
    }
  }
  // the iterators to over the letter counter maps
  for (auto & count : counter) {
    uint32_t bestQuality = 0;
    char bestBase = ' ';
    // if there is an insertion look at those if there is a majority of reads
    // with that insertion
    if (insertions.find(count.first) != insertions.end()) {
      for (auto & counterInsert : insertions[count.first]) {
      	counterInsert.second.getBest(bestBase, bestQuality, std::round(seqBase_.cnt_));
        if (bestBase == ' ') {
          continue;
        } else {
          calculatedConsensus.push_back(bestBase);
          calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
        }
      }
    }
    count.second.getBest(bestBase, bestQuality);
    if (bestBase == '-') {
      continue;
    }
    calculatedConsensus.push_back(bestBase);
    calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
  }
  bool debug = false;
  if(debug){
  	std::ofstream outFile;
  	openTextFile(outFile, "outFile.txt", ".txt", true, false);
  	outFile << "pos\tchar\tcount" << std::endl;
  	uint32_t posing = 0;
    for (auto & count : counter) {
      uint32_t bestQuality = 0;
      char bestBase = ' ';
      // if there is an insertion look at those if there is a majority of reads
      // with that insertion

      if (insertions.find(count.first) != insertions.end()) {
        for (auto & counterInsert : insertions[count.first]) {
        	std::cout << count.first << " : i: " << counterInsert.first << std::endl;
        	counterInsert.second.outPutInfo(std::cout, false);

        	for(const auto & a : counterInsert.second.alphabet_){
        		outFile << posing << "\t" << a << "\t" << counterInsert.second.chars_[a] << std::endl;
        	}
        	++posing;
        	//counterInsert.second.getBest(bestBase, bestQuality, std::round(seqBase_.cnt_));
          if (bestBase == ' ') {
            continue;
          } else {
            calculatedConsensus.push_back(bestBase);
            calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
          }
        }
      }
    	for(const auto & a : count.second.alphabet_){
    		outFile << posing << "\t" << a << "\t" << count.second.chars_[a] << std::endl;
    	}
      ++posing;
    	std::cout << count.first << std::endl;
    	count.second.outPutInfo(std::cout, false);
      //count.second.getBest(bestBase, bestQuality);
      if (bestBase == '-') {
        continue;
      }
      //calculatedConsensus.push_back(bestBase);
      //calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
    }
  }

  if (setToConsensus) {
    seqBase_.seq_ = calculatedConsensus;
    seqBase_.qual_ = calculatedConsensusQuality;
  }
  previousErrorChecks.clear();
  needToCalculateConsensus = false;
  // exit(1);
}
//////////// calculate the consensus
void baseCluster::calculateConsensus(aligner& alignerObj, bool setToConsensus) {
  // if the cluster is only one read, no need to create consensus
  if (seqBase_.cnt_ == 1 || reads_.size() <= 1) {
    return;
  }
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus) {
    return;
  }
  //std::cout << "creating consensus" << std::endl;
  // find the longest cluster
  // if the the current consensus is much shorter than the average amount of reads than choose the longest
  // read, this is to prevent consensus being built with only the a small portion of the sequences in the
  // the cluster
  auto averageSize = getAverageReadLength();
  readObject longestcluster;
  if (uAbsdiff(this->seqBase_.seq_.size(), averageSize) >
      0.1 * averageSize) {
  //if(true){
    uint64_t biggest = 0;
    for (const auto& read : reads_) {
      if (read.seqBase_.seq_.length() > biggest) {
        longestcluster = read;
        biggest = read.seqBase_.seq_.length();
      }
    }
  } else {
    longestcluster =
        readObject(seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_));
  }
  //this->seqBase_.outPutSeqAnsi(std::cout);
  //longestcluster.seqBase_.outPutSeqAnsi(std::cout);


  // clear the vectors in case of previous consensus calculations
  longestAlignments.clear();
  longestAlignmentsQualities.clear();
  longestAlingmentsRef.clear();
  calculatedConsensus.clear();
  calculatedConsensusQuality.clear();
  bool rStripReads = true;
  for (const auto & read : reads_) {
    alignerObj.alignVec(longestcluster, read, false);
    if(rStripReads){
    	seqUtil::rstripRead(alignerObj.alignObjectB_.seqBase_.seq_,
    			alignerObj.alignObjectB_.seqBase_.qual_, '-');
    }
    longestAlignments.push_back(alignerObj.alignObjectB_.seqBase_.seq_);
    longestAlignmentsQualities.push_back(
        alignerObj.alignObjectB_.seqBase_.qual_);
    longestAlingmentsRef.push_back(alignerObj.alignObjectA_.seqBase_.seq_);
    //alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
    //alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);

  }
  VecStr::iterator alignmentsItr;
  VecStr::iterator refAlignmentsItr;
  std::vector<std::vector<uint32_t>>::iterator alignmentsQaulitiesItr;
  std::vector<readObject>::iterator clusterIter;
  // create the map for letter counters for each position
  std::map<uint32_t, letterCounter> counter;
  // create a map in case of insertions
  std::map<uint32_t, std::map<uint32_t, letterCounter>> insertions;
  // counter of the reads
  //int numberOfReads = 0;
  // now iterator over the various clusters and alignments to count the bases
  // and their qualities
  // VecStr
  // stringsWithoutInsertions(longestAlingmentsRef.size(),"");
  int count = -1;
  for (alignmentsItr = longestAlignments.begin(),
      refAlignmentsItr = longestAlingmentsRef.begin(),
      alignmentsQaulitiesItr = longestAlignmentsQualities.begin(),
      clusterIter = reads_.begin();
       alignmentsItr != longestAlignments.end(); refAlignmentsItr++,
      alignmentsItr++, alignmentsQaulitiesItr++, clusterIter++) {
    ++count;
    // std::cout<<count+1<<"/"<<longestAlingmentsRef.size()<<std::endl;
    // the offset for the insertions
    uint32_t offSet = 0;
    // get the size of the cluster being processed
    //int sizeAdjustment = (int)clusterIter->seqBase_.cnt_;
    //numberOfReads += (int)clusterIter->seqBase_.cnt_;
    int currentOffset = 1;
    for (uint32_t i = 0; i < alignmentsItr->length(); i++) {
      // if the longest reference has an insertion in it put it in the
      // insertions letter counter map
      if ((*refAlignmentsItr)[i] == '-') {
        insertions[i - offSet][currentOffset]
            .letters[(*alignmentsItr)[i]] += clusterIter->seqBase_.cnt_;
        insertions[i - offSet][currentOffset]
            .qualities[(*alignmentsItr)[i]] +=
            (*alignmentsQaulitiesItr)[i] * clusterIter->seqBase_.cnt_;
        currentOffset++;
        offSet++;
        continue;
      }
      currentOffset = 1;
      counter[i - offSet].letters[(*alignmentsItr)[i]] +=
      		clusterIter->seqBase_.cnt_;
      counter[i - offSet].qualities[(*alignmentsItr)[i]] +=
          (*alignmentsQaulitiesItr)[i] * clusterIter->seqBase_.cnt_;
      addOtherVec(
          counter[i - offSet].allQualities[(*alignmentsItr)[i]],
          std::vector<uint32_t>(clusterIter->seqBase_.cnt_, (*alignmentsQaulitiesItr)[i]));
    }
  }
  // the iterators to over the letter counter maps
  //std::cout << "num of reads " << numberOfReads << std::endl;
  std::map<uint32_t, letterCounter>::iterator counterIter;
  std::map<uint32_t, letterCounter>::iterator insertionsIter;
  for (counterIter = counter.begin(); counterIter != counter.end();
       counterIter++) {
    uint32_t bestQuality = 0;
    char bestBase = ' ';
    // if there is an insertion look at those if there is a majority of reads
    // with that insertion
    if (insertions.find(counterIter->first) != insertions.end()) {
      for (insertionsIter = insertions[counterIter->first].begin();
           insertionsIter != insertions[counterIter->first].end();
           insertionsIter++) {
        insertionsIter->second.getBest(bestBase, bestQuality, std::round(seqBase_.cnt_));
        if (bestBase == ' ') {
          continue;
        } else {
          calculatedConsensus.push_back(bestBase);
          calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
        }
      }
    }
    counterIter->second.getBest(bestBase, bestQuality);
    //printOutMapContents(counterIter->second.letters, "\t", std::cout);
    //std::cout << calculatedConsensus.size() << ":" << bestBase << std::endl;
    double fortyPercent = 0.40 * this->seqBase_.cnt_;
    if (bestBase == '-' || counterIter->second.letters[bestBase] < fortyPercent) {
      continue;
    }
    calculatedConsensus.push_back(bestBase);
    calculatedConsensusQuality.push_back(bestQuality / counterIter->second.getTotalCount());
  }

  if (setToConsensus) {
    seqBase_.seq_ = calculatedConsensus;
    seqBase_.qual_ = calculatedConsensusQuality;
  }
  previousErrorChecks.clear();
  needToCalculateConsensus = false;
  //this->seqBase_.outPutSeqAnsi(std::cout);
  //std::cout << std::endl;
  /*if(this->seqBase_.cnt_ > 50){
  	std::ofstream tempfile("tempFilealns.fasta");
  	for(uint64_t pos : iter::range(longestAlignments.size())){
  		tempfile << ">" << pos << std::endl;
  		tempfile << longestAlignments[pos] << std::endl;
  		tempfile << ">" << pos << "_ref" << std::endl;
  		tempfile << longestAlingmentsRef[pos] << std::endl;
  	}
  }*/
  // exit(1);
}

void baseCluster::writeOutLongestAlignments(
    const std::string& workingDirectory) {

  std::stringstream longestAlignmentFilename;
  std::stringstream longestAlignmentQualFilename;
  longestAlignmentFilename << workingDirectory << seqBase_.name_
                           << "_long.fasta";
  longestAlignmentQualFilename << longestAlignmentFilename.str() << ".qual";
  std::ofstream longestAlignmentFile(longestAlignmentFilename.str().c_str());
  std::ofstream longestAlignmentQualFile(
      longestAlignmentQualFilename.str().c_str());
  VecStr::iterator aligns = longestAlignments.begin();
  VecStr::iterator refAligns = longestAlingmentsRef.begin();
  std::vector<readObject>::iterator clusIter = reads_.begin();
  for (; aligns != longestAlignments.end(); ++aligns, ++refAligns, ++clusIter) {
    longestAlignmentFile << ">" << clusIter->seqBase_.name_ << std::endl;
    longestAlignmentFile << *aligns << std::endl;
    longestAlignmentFile << ">" << clusIter->seqBase_.name_ << "_"
                         << seqBase_.name_ << std::endl;
    longestAlignmentFile << *refAligns << std::endl;
  }
  clusIter = reads_.begin();
  for (std::vector<std::vector<uint32_t>>::iterator
           qualIter = longestAlignmentsQualities.begin();
       qualIter != longestAlignmentsQualities.end(); ++qualIter, ++clusIter) {
    longestAlignmentQualFile << ">" << clusIter->seqBase_.name_ << std::endl;
    longestAlignmentQualFile << vectorToString(*qualIter) << std::endl;
    ;
  }
}


/*
 * old consensus building
 * void baseCluster::calculateConsensus(aligner& alignerObj, bool setToConsensus) {
  // if the cluster is only one read, no need to create consensus
  if (seqBase_.cnt_ == 1 || reads_.size() <= 1) {
    return;
  }
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus) {
    return;
  }
  //std::cout << "creating consensus" << std::endl;
  // find the longest cluster
  // if the the current consensus is much shorter than the average amount of reads than choose the longest
  // read, this is to prevent consensus being built with only the a small portion of the sequences in the
  // the cluster
  auto averageSize = getAverageReadLength();
  readObject longestcluster;
  if (uAbsdiff(this->seqBase_.seq_.size(), averageSize) >
      0.1 * averageSize) {
  //if(true){
    uint64_t biggest = 0;
    for (const auto& read : reads_) {
      if (read.seqBase_.seq_.length() > biggest) {
        longestcluster = read;
        biggest = read.seqBase_.seq_.length();
      }
    }
  } else {
    longestcluster =
        readObject(seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_));
  }
  this->seqBase_.outPutSeqAnsi(std::cout);
  longestcluster.seqBase_.outPutSeqAnsi(std::cout);


  // clear the vectors in case of previous consensus calculations
  longestAlignments.clear();
  longestAlignmentsQualities.clear();
  longestAlingmentsRef.clear();
  calculatedConsensus.clear();
  calculatedConsensusQuality.clear();
  for (const auto & read : reads_) {
    alignerObj.alignVec(longestcluster, read, false);
    longestAlignments.push_back(alignerObj.alignObjectB_.seqBase_.seq_);
    longestAlignmentsQualities.push_back(
        alignerObj.alignObjectB_.seqBase_.qual_);
    longestAlingmentsRef.push_back(alignerObj.alignObjectA_.seqBase_.seq_);
    //alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
    //alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);

  }
  VecStr::iterator alignmentsItr;
  VecStr::iterator refAlignmentsItr;
  std::vector<std::vector<uint32_t>>::iterator alignmentsQaulitiesItr;
  std::vector<readObject>::iterator clusterIter;
  // create the map for letter counters for each position
  std::map<uint32_t, letterCounter> counter;
  // create a map in case of insertions
  std::map<uint32_t, std::map<uint32_t, letterCounter>> insertions;
  // counter of the reads
  //int numberOfReads = 0;
  // now iterator over the various clusters and alignments to count the bases
  // and their qualities
  // VecStr
  // stringsWithoutInsertions(longestAlingmentsRef.size(),"");
  int count = -1;
  for (alignmentsItr = longestAlignments.begin(),
      refAlignmentsItr = longestAlingmentsRef.begin(),
      alignmentsQaulitiesItr = longestAlignmentsQualities.begin(),
      clusterIter = reads_.begin();
       alignmentsItr != longestAlignments.end(); refAlignmentsItr++,
      alignmentsItr++, alignmentsQaulitiesItr++, clusterIter++) {
    ++count;
    // std::cout<<count+1<<"/"<<longestAlingmentsRef.size()<<std::endl;
    // the offset for the insertions
    uint32_t offSet = 0;
    // get the size of the cluster being processed
    //int sizeAdjustment = (int)clusterIter->seqBase_.cnt_;
    //numberOfReads += (int)clusterIter->seqBase_.cnt_;
    int currentOffset = 1;
    for (uint32_t i = 0; i < alignmentsItr->length(); i++) {
      // if the longest reference has an insertion in it put it in the
      // insertions letter counter map
      if ((*refAlignmentsItr)[i] == '-') {
        insertions[i - offSet][currentOffset]
            .letters[(*alignmentsItr)[i]] += clusterIter->seqBase_.cnt_;
        insertions[i - offSet][currentOffset]
            .qualities[(*alignmentsItr)[i]] +=
            (*alignmentsQaulitiesItr)[i] * clusterIter->seqBase_.cnt_;
        currentOffset++;
        offSet++;
        continue;
      }
      currentOffset = 1;
      counter[i - offSet].letters[(*alignmentsItr)[i]] +=
      		clusterIter->seqBase_.cnt_;
      counter[i - offSet].qualities[(*alignmentsItr)[i]] +=
          (*alignmentsQaulitiesItr)[i] * clusterIter->seqBase_.cnt_;
      addOtherVec(
          counter[i - offSet].allQualities[(*alignmentsItr)[i]],
          std::vector<uint32_t>(clusterIter->seqBase_.cnt_, (*alignmentsQaulitiesItr)[i]));
    }
  }
  // the iterators to over the letter counter maps
  //std::cout << "num of reads " << numberOfReads << std::endl;
  std::map<uint32_t, letterCounter>::iterator counterIter;
  std::map<uint32_t, letterCounter>::iterator insertionsIter;
  for (counterIter = counter.begin(); counterIter != counter.end();
       counterIter++) {
    uint32_t bestQuality = 0;
    char bestBase = ' ';
    // if there is an insertion look at those if there is a majority of reads
    // with that insertion
    if (insertions.find(counterIter->first) != insertions.end()) {
      for (insertionsIter = insertions[counterIter->first].begin();
           insertionsIter != insertions[counterIter->first].end();
           insertionsIter++) {
        insertionsIter->second.getBest(bestBase, bestQuality, std::round(seqBase_.cnt_));
        if (bestBase == ' ') {
          continue;
        } else {
          calculatedConsensus.push_back(bestBase);
          calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
        }
      }
    }
    counterIter->second.getBest(bestBase, bestQuality);
    printOutMapContents(counterIter->second.letters, "\t", std::cout);
    std::cout << calculatedConsensus.size() << ":" << bestBase << std::endl;
    if (bestBase == '-') {
      continue;
    }
    calculatedConsensus.push_back(bestBase);
    calculatedConsensusQuality.push_back(bestQuality / seqBase_.cnt_);
  }

  if (setToConsensus) {
    seqBase_.seq_ = calculatedConsensus;
    seqBase_.qual_ = calculatedConsensusQuality;
  }
  previousErrorChecks.clear();
  needToCalculateConsensus = false;
  this->seqBase_.outPutSeqAnsi(std::cout);
  std::cout << std::endl;
  if(this->seqBase_.cnt_ > 50){
  	std::ofstream tempfile("tempFilealns.fasta");
  	for(uint64_t pos : iter::range(longestAlignments.size())){
  		tempfile << ">" << pos << std::endl;
  		tempfile << longestAlignments[pos] << std::endl;
  		tempfile << ">" << pos << "_ref" << std::endl;
  		tempfile << longestAlingmentsRef[pos] << std::endl;
  	}
  }
  // exit(1);
}
 */
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

  readObjectIO reader = readObjectIO();
  reader.writeSimpleVector(alignments.first, alignMentsFileName.str(),
                           "fastaQual", false, false);
  reader.writeSimpleVector(alignments.second, alignmentsOnlyFileName.str(),
                           "fastaQual", false, false);
}

/// output the clusters currently clustered to the this cluster
void baseCluster::writeOutClusters(const std::string& directoryName) const {
  if (reads_.size() == 0) {

  } else {
    readObjectIO reader = readObjectIO();
    ;
    reader.writeVector(reads_, directoryName + seqBase_.name_, "fastaQual",
                       false, false);
  }
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
  readObject outRead(seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_));
  outRead.seqBase_.cnt_ = seqBase_.cnt_;
  outRead.seqBase_.frac_ = seqBase_.frac_;
  outRead.normalizedFraction = normalizedFraction;
  return outRead;
}

void baseCluster::printDescription(std::ostream& out, bool deep) const {
  readObject::printDescription(out, deep);
  out << "baseCluster{" << std::endl << "firstReadName:" << firstReadName
      << std::endl << "firstReadCount:" << firstReadCount << std::endl;
  if (deep) {
    out << "reads:" << std::endl;
    out << "std::vector<bib::readObject>" << std::endl;
    for (const auto& read : reads_) {
      read.printDescription(out, deep);
    }
    out << "previousErrorChecks" << std::endl;
    out << "std::map<std::string, errorProfile>" << std::endl;
    for (const auto& profile : previousErrorChecks) {
      out << "\t" << profile.first << ":";
      profile.second.printDescription(out, deep);
    }
  }
  out << "longestAlignments:" << longestAlignments << std::endl
      << "longestAlingmentsRef:" << longestAlingmentsRef << std::endl;
  out << "longestAlignmentsQualities:" << std::endl;
  out << "std::vector<std::vector<uint32_t>>" << std::endl;
  for (const auto& subVec : longestAlignmentsQualities) {
    out << "\t" << subVec << std::endl;
  }
  out << "needToCalculateConsensus:" << needToCalculateConsensus << std::endl
      << "calculatedConsensus:" << calculatedConsensus << std::endl
      << "calculatedConsensusQuality:" << calculatedConsensusQuality
      << std::endl;
}

}  // namespace bib
