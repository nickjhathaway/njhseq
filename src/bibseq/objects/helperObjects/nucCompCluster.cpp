/*
 * nucCompCluster.cpp
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
#include "nucCompCluster.hpp"

namespace bibseq {
void nucCompCluster::addOtherCluster(const nucCompCluster & otherCluster,
		bool setFraction){
	counter_.addOtherCounts(otherCluster.counter_, setFraction);
	addOtherVec(readPositions_, otherCluster.readPositions_);
	readCnt_+= otherCluster.readCnt_;
}
uint32_t nucCompCluster::readBufferMax_ = 10;


void collapseSimilarNucCompClusters(std::vector<nucCompCluster>& comps,
		double diffCutOff, bool verbose){
  std::vector<std::tuple<uint32_t, uint32_t, double>> similarComps;
  bool keepRemoving = true;
  while(keepRemoving){
    std::vector<uint32_t> clusterPositions(comps.size());
    iota<uint32_t>(clusterPositions, 0);
    reverse(clusterPositions);
    std::vector<uint32_t> removeThese;
    for(const auto & compPos : clusterPositions){
  		double smallestDiff = std::numeric_limits<double>::max();
  		uint32_t closestPos = std::numeric_limits<uint32_t>::max();
    	for(const auto & compPosSecond : iter::range<uint32_t>(0, compPos)){
    		if(compPos == compPosSecond){
    			continue;
    		}
    		double currentDiff = comps[compPos].counter_.getFracDifference(comps[compPosSecond].counter_, comps[compPosSecond].counter_.alphabet_);
  			if(currentDiff < diffCutOff){
  				if(currentDiff < smallestDiff){
  					smallestDiff = currentDiff;
  					closestPos = compPosSecond;
  				}
  			}
    	}
    	if(closestPos != std::numeric_limits<uint32_t>::max()){
    		comps[closestPos].addOtherCluster(comps[compPos], true);
    		removeThese.emplace_back(compPos);
    		if(verbose){
    			std::cout << compPos << ":" << closestPos << " by : "
    					<< smallestDiff
    					<< "\n";
    		}
    	}
    }
    if(removeThese.empty()){
    	keepRemoving = false;
    }else{
    	if(verbose){
    		printVector(removeThese);
    	}
      for(const auto & pos : removeThese){
      	if(verbose){
        	std::cout << pos << "\n";
      	}
      	comps.erase(comps.begin() + pos);
      }
    }
  }
  sortNucCompVec(comps);
}
void sortNucCompVec(std::vector<nucCompCluster>& comps){
  sort(comps, [&](const nucCompCluster & comp1,
  		const nucCompCluster & comp2){ return comp1.readCnt_ > comp2.readCnt_;});
}
table getInfoNucCompVec(std::vector<nucCompCluster>& comps){
  table outInfo(VecStr{"compPos","readNum","readFrac","Afrac","Cfrac","Gfrac","Tfrac"});
  double totalCount = 0;
  for(const auto & compPos : iter::range(comps.size())){
  	totalCount += comps[compPos].readCnt_;
  }
  for(const auto & compPos : iter::range(comps.size())){
  	comps[compPos].counter_.setFractions();
  	comps[compPos].readBuffer_ = 0;
  	std::stringstream tempStream;
  	tempStream << compPos << "\t" << comps[compPos].readCnt_ << "\t" << comps[compPos].readCnt_/totalCount <<
  		 "\t" << comps[compPos].counter_.fractions_['A'] <<
  	   "\t" << comps[compPos].counter_.fractions_['C'] <<
  	   "\t" << comps[compPos].counter_.fractions_['G'] <<
  	   "\t" << comps[compPos].counter_.fractions_['T'] << "\n";
  	outInfo.content_.emplace_back(stringToVector<std::string>(tempStream.str()));
  }
  return outInfo;
}

std::vector<nucCompCluster>  clusterOnNucComp(const SeqIOOptions & ioOptions,
		const std::vector<char> & alphabet,
		double diffCutOff, bool findBest, bool preSort,
		bool verbose){
	std::vector<nucCompCluster> comps;
	SeqInput reader(ioOptions);
	reader.openIn();
	readObject read;
  uint64_t pos = reader.tellgPri();
  uint32_t count = 0;

  if (SeqIOOptions::inFormats::FASTQ != ioOptions.inFormat_
  		&& SeqIOOptions::inFormats::FASTA != ioOptions.inFormat_){
  	std::stringstream ss;
  	ss << "clusterOnNucComp\n";
  	ss << "Currently only works on fastq files, improper format argument given\n";
  	throw std::runtime_error{ss.str()};
  }
	while (reader.readNextRead(read)) {
		read.setLetterCount(alphabet);
		if (count % 5000 == 0 && verbose) {
			std::cout << "On " << count << "\r";
			std::cout.flush();
		}
		findBestNucComp(read, pos, alphabet, diffCutOff, findBest, comps);
		pos = reader.tellgPri();
		++count;
	}
	std::cout << std::endl;


  sortNucCompVec(comps);
  collapseSimilarNucCompClusters(comps, diffCutOff, verbose);
  return comps;
}
} /* namespace bib */
