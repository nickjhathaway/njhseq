#pragma once
/*

 * nucCompCluster.hpp
 *
 *  Created on: Apr 9, 2014
 *      Author: nickhathaway
 */
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
#include "bibseq/utils.h"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp"
#include "bibseq/IO/SeqIO/SeqInput.hpp"

namespace bibseq {

class nucCompCluster{
public:
	nucCompCluster(){};

	template<typename T>
	nucCompCluster(const T& firstRead, uint64_t readPos):
	counter_(charCounter(firstRead.counter_.alphabet_)),
	readCnt_(firstRead.seqBase_.cnt_) {
		counter_.increaseCountByString(firstRead.seqBase_.seq_, firstRead.seqBase_.cnt_);
		readPositions_.emplace_back(readPos);
		counter_.setFractions();
	}

	template<typename T>
	nucCompCluster(const T& firstRead, uint64_t readPos, uint64_t minLen):
	counter_(charCounter(firstRead.counter_.alphabet_)),
	readCnt_(firstRead.seqBase_.cnt_) {
		counter_.increasePortion(firstRead.seqBase_.seq_,minLen, firstRead.seqBase_.cnt_);
		readPositions_.emplace_back(readPos);
		counter_.setFractions();
	}

	charCounter counter_;
	std::vector<uint64_t> readPositions_;
	double readCnt_ = 0;
	uint32_t readBuffer_ = 0;

	static uint32_t readBufferMax_;

	template<typename T>
	bool compareRead(const T& read, double diffCutOff){
		double sum = counter_.getFracDifference(read.counter_, counter_.alphabet_);
		if(sum > diffCutOff){
			return false;
		}else{
			return true;
		}
	}

	template<typename T>
	void addRead(const T& read, uint64_t readPos){
		readCnt_ += read.seqBase_.cnt_;
		counter_.increaseCountByString(read.seqBase_.seq_, read.seqBase_.cnt_);
		readPositions_.emplace_back(readPos);
		++readBuffer_;
		if(readBuffer_ > readBufferMax_){
			readBuffer_ = 0;
			counter_.setFractions();
		}
	}
	template<typename T>
	void addRead(const T& read, uint64_t readPos, uint64_t minLen){
		readCnt_ += read.seqBase_.cnt_;
		counter_.increasePortion(read.seqBase_.seq_,minLen, read.seqBase_.cnt_);
		readPositions_.emplace_back(readPos);
		++readBuffer_;
		if(readBuffer_ > readBufferMax_){
			readBuffer_ = 0;
			counter_.setFractions();
		}
	}
	/**@brief Add another cluster that is similar
	 *
	 * @param otherCluster The other nucCompCluster to add
	 * @param setFraction Whether to reset the fraction counts again after adding other cluster
	 */
	void addOtherCluster(const nucCompCluster & otherCluster, bool setFraction);
	template<typename T>
	std::vector<T> getReads(const std::vector<T> & reads)const{
		std::vector<T> outReads;
		outReads.reserve(readPositions_.size());
		for(const auto & pos : readPositions_){
			outReads.emplace_back(reads[pos]);
		}
		return outReads;
	}
	template<typename T>
	std::vector<T> getReads(const SeqIOOptions &
			ioOptions) const{
		std::vector<T> outReads;
		outReads.reserve(readPositions_.size());
		SeqInput reader(ioOptions);
		reader.openIn();
		readObject read;
		for(const auto & pos : readPositions_){
			reader.seekgPri(pos);
			reader.readNextRead(read);
			outReads.emplace_back(read);
		}
	  return outReads;
	}

};
template<typename T>
void findBestNucComp(const T & read, uint64_t pos, const std::vector<char> & alphabet,
		double diffCutOff, bool findBest, std::vector<nucCompCluster> & comps){
	bool found = false;
	double lowestDifference = 1;
	uint64_t bestPos = std::numeric_limits<uint64_t>::max();
  for(const auto & compPos : iter::range(comps.size())){
  	if(comps[compPos].compareRead(read, diffCutOff)){
  		found = true;
  		if(findBest){
  			double currentDiff = comps[compPos].counter_.getFracDifference(read.counter_, alphabet);
  			if(currentDiff < lowestDifference){
  				lowestDifference = currentDiff;
  				bestPos = compPos;
  			}
  		}else{
  			comps[compPos].addRead(read, pos);
    		break;
  		}
  	}
  }
  if(!found){
  	comps.emplace_back(nucCompCluster(read, pos));
  }else if(findBest){
  	comps[bestPos].addRead(read, pos);
  }
}
template<typename T>
void findBestNucComp(const T & read, uint64_t pos, uint64_t minLen,const std::vector<char> & alphabet,
		double diffCutOff, bool findBest, std::vector<nucCompCluster> & comps){
	bool found = false;
	double lowestDifference = 1;
	uint64_t bestPos = std::numeric_limits<uint64_t>::max();
  for(const auto & compPos : iter::range(comps.size())){
  	if(comps[compPos].compareRead(read, diffCutOff)){
  		found = true;
  		if(findBest){
  			double currentDiff = comps[compPos].counter_.getFracDifference(read.counter_, alphabet);
  			if(currentDiff < lowestDifference){
  				lowestDifference = currentDiff;
  				bestPos = compPos;
  			}
  		}else{
  			comps[compPos].addRead(read, pos, minLen);
    		break;
  		}
  	}
  }
  if(!found){
  	comps.emplace_back(nucCompCluster(read, pos, minLen));
  }else if(findBest){
  	comps[bestPos].addRead(read, pos, minLen);
  }
}

void collapseSimilarNucCompClusters(std::vector<nucCompCluster>& comps,
		double diffCutOff, bool verbose);
void sortNucCompVec(std::vector<nucCompCluster>& comps);
table getInfoNucCompVec(std::vector<nucCompCluster>& comps);

template<typename T>
std::vector<nucCompCluster>  clusterOnNucComp(std::vector<T> & reads,
		const std::vector<char> & alphabet,
		double diffCutOff, bool findBest, bool preSort,
		bool verbose){
	std::vector<nucCompCluster> comps;
	uint64_t pos = 0;
  readVecSorter::sortReadVector(reads, "seqCondensed");
	for(auto & read : reads ){
		read.setLetterCount(alphabet);
  	if(pos % 5000 == 0 && verbose){
  		std::cout << "On " << pos << std::endl;
  	}
  	findBestNucComp(read, pos, alphabet, diffCutOff, findBest, comps);
  	++pos;
	}
  sortNucCompVec(comps);
  collapseSimilarNucCompClusters(comps, diffCutOff, verbose);
  return comps;
}
template<typename T>
std::vector<nucCompCluster>  clusterOnNucComp(std::vector<T> & reads,
		uint64_t minLen,
		const std::vector<char> & alphabet,
		double diffCutOff, bool findBest, bool preSort,
		bool verbose){
	std::vector<nucCompCluster> comps;
	uint64_t pos = 0;
  readVecSorter::sortReadVector(reads, "seqCondensed");
	for(auto & read : reads ){
		read.setLetterCount(alphabet);
  	if(pos % 5000 == 0 && verbose){
  		std::cout << "On " << pos << std::endl;
  	}
  	findBestNucComp(read, pos,minLen, alphabet, diffCutOff, findBest, comps);
  	++pos;
	}
  sortNucCompVec(comps);
  collapseSimilarNucCompClusters(comps, diffCutOff, verbose);
  return comps;
}

std::vector<nucCompCluster>  clusterOnNucComp(const SeqIOOptions &
		ioOptions,
		const std::vector<char> & alphabet,
		double diffCutOff, bool findBest, bool preSort,
		bool verbose);



} /* namespace bibseq */



