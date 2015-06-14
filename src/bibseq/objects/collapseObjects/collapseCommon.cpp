//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//  collapseCommon.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/7/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "collapseCommon.hpp"

namespace bibseq {
namespace collapse {
void clusterSet::setSubClustersPositions() {
  uint32_t pos = 0;
  for (const auto& clus : clusters_) {
    for (const auto& read : clus.reads_) {
      subClustersPositions_[read.getStubName(true)] = pos;
    }
    ++pos;
  }
}

void clusterSet::setSetInfo() {
	info_.clear();
  for (const auto& read : clusters_) {
  	info_.updateInfo(read.reads_);
  }
}

// writing
void clusterSet::writeClusters(std::string filename, std::string format,
                               bool overWrite, bool exitOnFailure) const {
  readObjectIO::write(clusters_, filename, format, overWrite, exitOnFailure);
}

void clusterSet::writeClusters(std::string filename, const readObjectIOOptions & ioOptions) const{
	writeClusters(filename, ioOptions.outFormat_, ioOptions.overWriteFile_, ioOptions.exitOnFailureToWrite_);
}

table clusterSet::getReplicateInfo()const {
	VecStr colNames{"clusterName", "repName", "fraction"};
	table ret(colNames);
	for(const auto & clus : clusters_){
		for(const auto & s : info_.infos_){
			auto search = clus.sampInfos().find(s.first);
			if(search == clus.sampInfos().end()){
				ret.content_.emplace_back(toVecStr(clus.seqBase_.name_, s.first, 0));
			}else{
				ret.content_.emplace_back(toVecStr(clus.seqBase_.name_, s.first, search->second.fraction_));
			}
		}
	}
	return ret;
}

double clusterSet::getRMSE()const{
	double sumSquares = 0;

	for(const auto & clus : clusters_){
		auto sampInfosNames = getVectorOfMapKeys(clus.sampInfos());
		std::vector<double> differences;
		for(const auto & pos : iter::range(sampInfosNames.size())){
			for(const auto & subPos : iter::range(pos)){
				if(pos == subPos){
					continue;
				}
				differences.emplace_back(clus.sampInfos().at(sampInfosNames[pos]).fraction_ - clus.sampInfos().at(sampInfosNames[subPos]).fraction_);
			}
		}
		sumSquares += std::pow(vectorMean(differences), 2);
	}
	return std::sqrt(sumSquares/clusters_.size());
}

}  // collapse
}  // bib
