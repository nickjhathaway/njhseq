/*
 * PairedCluster.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: nick
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
#include "PairedCluster.hpp"
#include "bibseq/helpers/consensusHelper.hpp"
#include "bibseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp"
#include "bibseq/IO/SeqIO/SeqOutput.hpp"


namespace bibseq {

void PairedCluster::addRead(const PairedCluster & otherRead){
	addOtherVec(reads_, otherRead.reads_);
  seqBase_.cnt_ += otherRead.seqBase_.cnt_;
  seqBase_.frac_ = seqBase_.frac_ * reads_.size();
  seqBase_.frac_ = (seqBase_.frac_ + otherRead.seqBase_.frac_) / reads_.size();
  mateSeqBase_.cnt_ = seqBase_.cnt_;
  mateSeqBase_.frac_ = seqBase_.frac_;
  // needToCalculateConsensus = true;
  if (seqBase_.cnt_ / 2 > firstReadCount) {
    needToCalculateConsensus = true;
  }
}

void PairedCluster::calculateConsensus(aligner & alignerObj, bool setToConsensus){
	//std::cout << "calculateConsensus start" << std::endl;
	// if the cluster is only one read, no need to create consensus
	if (reads_.size() <= 1) {
		//std::cout << "calculateConsensus stop" << std::endl;
		return;
	}
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus) {
  	//std::cout << "calculateConsensus stop" << std::endl;
    return;
  }

	/**@todo still needs to be validated*/
  std::function<const seqInfo& (const std::shared_ptr<PairedRead>& )> getSeqBase =
			[](const std::shared_ptr<PairedRead>& read) ->const seqInfo& {return read->seqBase_;};
	std::function<const seqInfo& (const std::shared_ptr<PairedRead>& )> getSeqBaseMate =
			[](const std::shared_ptr<PairedRead> & read) ->const seqInfo& {return read->mateSeqBase_;};

  /*
  std::function<const seqInfo& (const PairedCluster& )> getSeqBase =
			[](const PairedCluster& read) ->const seqInfo& {return std::cref(read.seqBase_);};
	std::function<const seqInfo& (const PairedCluster& )> getSeqBaseMate =
			[](const PairedCluster & read) ->const seqInfo& {return std::cref(read.mateSeqBase_);};*/

  /*
	std::function<seqInfo (const PairedCluster& )> getSeqBase =
			[](const PairedCluster& read) -> seqInfo {return read.seqBase_;};
	std::function<seqInfo (const PairedCluster& )> getSeqBaseMate =
			[](const PairedCluster & read) ->seqInfo {return read.mateSeqBase_;};
			*/
	/*
	std::cout << "calculateConsensus inbetween1" << std::endl;
	for(const auto & read : reads_){
		const auto & rInfo = getSeqBase(read);
		std::cout << rInfo.seq_.size() << std::endl;
		std::cout << rInfo.name_ << std::endl;
		std::cout << rInfo.seq_ << std::endl;
		std::cout << rInfo.qual_ << std::endl;
		rInfo.outPutFastq(std::cout);
		read.seqBase_.outPutFastq(std::cout);
		read.mateSeqBase_.outPutFastq(std::cout);
	}
	std::cout << std::endl;
	const auto & rInfo = getSeqBase(*this);
	std::cout << rInfo.seq_.size() << std::endl;
	rInfo.outPutFastq(std::cout);
	std::cout << "calculateConsensus inbetween1.5" << std::endl;
	std::cout << alignerObj.parts_.maxSize_ << std::endl;
*/
	auto seqBaseCon = consensusHelper::buildConsensus(reads_, getSeqBase, alignerObj, seqBase_.name_);
	//std::cout << "calculateConsensus inbetween2" << std::endl;
	auto mateSeqBaseCon =  consensusHelper::buildConsensus(reads_, getSeqBaseMate, alignerObj, mateSeqBase_.name_);
	//std::cout << "calculateConsensus inbetween3" << std::endl;
	if(setToConsensus){
		seqBase_ = seqBaseCon;
		mateSeqBase_ = mateSeqBaseCon;
	}
	//std::cout << "calculateConsensus inbetween4" << std::endl;
	setLetterCount();
	//std::cout << "calculateConsensus inbetween5" << std::endl;
	needToCalculateConsensus = false;
	//std::cout << "calculateConsensus stop" << std::endl;
}


bool PairedCluster::compare(PairedCluster & otherRead,
		aligner & alignerObj,
		const comparison & errorThreshold,
		const collapserOpts & collapserOptsObj){
	return compareRead(otherRead, alignerObj, errorThreshold, collapserOptsObj)
			&& compareMate(otherRead, alignerObj, errorThreshold, collapserOptsObj);
}


bool PairedCluster::compareRead(PairedCluster & otherRead,
		aligner & alignerObj,
		const comparison & errorThreshold,
		const collapserOpts & collapserOptsObj){
	//std::cout << "compare read start" << std::endl;
  if(collapserOptsObj.noAlign_){
  	alignerObj.noAlignSetAndScore(seqBase_, otherRead.seqBase_);
  }else{
  	alignerObj.alignCacheGlobal(seqBase_, otherRead.seqBase_);
  }
	//alignerObj.alignVec(seqBase_, otherRead.seqBase_, collapserOptsObj.local_);
	alignerObj.profileAlignment(seqBase_, otherRead.seqBase_,
			collapserOptsObj.kLength_,collapserOptsObj.kmersByPosition_,
			collapserOptsObj.checkKmers_, true, false, collapserOptsObj.weighHomopolyer_);
	//std::cout << "compare read stop" << std::endl;
	return errorThreshold.passErrorProfile(alignerObj.comp_);
}

bool PairedCluster::compareMate(PairedCluster & otherRead,
		aligner & alignerObj,
		const comparison & errorThreshold,
		const collapserOpts & collapserOptsObj){
	/*
	std::cout << "compareMate start" << std::endl;
	std::cout << alignerObj.parts_.maxSize_ << std::endl;
	std::cout << mateSeqBase_.seq_.size( ) << std::endl;
	std::cout << otherRead.seqBase_.seq_.size( ) << std::endl;*/
  if(collapserOptsObj.noAlign_){
  	alignerObj.noAlignSetAndScore(mateSeqBase_, otherRead.mateSeqBase_);
  }else{
  	alignerObj.alignCacheGlobal(mateSeqBase_, otherRead.mateSeqBase_);
  }
	//alignerObj.alignVec(mateSeqBase_, otherRead.mateSeqBase_, collapserOptsObj.local_);
	//std::cout << "compareMate inbetween" << std::endl;
	alignerObj.profileAlignment(mateSeqBase_, otherRead.mateSeqBase_,
			collapserOptsObj.kLength_,collapserOptsObj.kmersByPosition_,
			collapserOptsObj.checkKmers_, true, false, collapserOptsObj.weighHomopolyer_);
	//std::cout << "compareMate stop" << std::endl;
	return errorThreshold.passErrorProfile(alignerObj.comp_);
}

bool PairedCluster::compareId(PairedCluster & otherRead,
		aligner & alignerObj,
		const comparison & errorThreshold,
		const collapserOpts & collapserOptsObj){
	return compareReadId(otherRead, alignerObj, errorThreshold, collapserOptsObj)
			&& compareMateId(otherRead, alignerObj, errorThreshold, collapserOptsObj);
}


bool PairedCluster::compareReadId(PairedCluster & otherRead,
		aligner & alignerObj,
		const comparison & errorThreshold,
		const collapserOpts & collapserOptsObj){
	//std::cout << "compare read start" << std::endl;
  if(collapserOptsObj.noAlign_){
  	alignerObj.noAlignSetAndScore(seqBase_, otherRead.seqBase_);
  }else{
  	alignerObj.alignCacheGlobal(seqBase_, otherRead.seqBase_);
  }
	//alignerObj.alignVec(seqBase_, otherRead.seqBase_, collapserOptsObj.local_);
	alignerObj.profileAlignment(seqBase_, otherRead.seqBase_,
			collapserOptsObj.kLength_,collapserOptsObj.kmersByPosition_,
			collapserOptsObj.checkKmers_, true, false, collapserOptsObj.weighHomopolyer_);
	//std::cout << "compare read stop" << std::endl;
	return errorThreshold.passIdThreshold(alignerObj.comp_);

}

bool PairedCluster::compareMateId(PairedCluster & otherRead,
		aligner & alignerObj,
		const comparison & errorThreshold,
		const collapserOpts & collapserOptsObj){
	/*
	std::cout << "compareMate start" << std::endl;
	std::cout << alignerObj.parts_.maxSize_ << std::endl;
	std::cout << mateSeqBase_.seq_.size( ) << std::endl;
	std::cout << otherRead.seqBase_.seq_.size( ) << std::endl;*/
  if(collapserOptsObj.noAlign_){
  	alignerObj.noAlignSetAndScore(mateSeqBase_, otherRead.mateSeqBase_);
  }else{
  	alignerObj.alignCacheGlobal(mateSeqBase_, otherRead.mateSeqBase_);
  }
	//alignerObj.alignVec(mateSeqBase_, otherRead.mateSeqBase_, collapserOptsObj.local_);
	//std::cout << "compareMate inbetween" << std::endl;
	alignerObj.profileAlignment(mateSeqBase_, otherRead.mateSeqBase_,
			collapserOptsObj.kLength_,collapserOptsObj.kmersByPosition_,
			collapserOptsObj.checkKmers_, true, false, collapserOptsObj.weighHomopolyer_);
	//std::cout << "compareMate stop" << std::endl;
	return errorThreshold.passIdThreshold(alignerObj.comp_);
}

void collapseIdenticalPairedClusters(std::vector<PairedRead> & clusters,
		std::string qualRep){
	readVecSorter::sortReadVector(clusters, "seq", true);
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<uint64_t>>> identicalPos;
	for(const auto & clusPos : iter::range(clusters.size())){
		const auto & clus = clusters[clusPos];
		identicalPos[clus.seqBase_.seq_][clus.mateSeqBase_.seq_].emplace_back(clusPos);
	}
	std::vector<uint64_t> removeThese;
	for(const auto & iden : identicalPos){
		for(const auto & sub : iden.second){
			if(sub.second.size() > 1){
				for(const auto & posItem : iter::enumerate(sub.second)){
					if(posItem.index == 0){
						continue;
					}
					clusters[sub.second.front()].seqBase_.cnt_ += clusters[posItem.element].seqBase_.cnt_;
					clusters[sub.second.front()].mateSeqBase_.cnt_ += clusters[posItem.element].mateSeqBase_.cnt_;
					//clusters[sub.second.front()].reads_.front().seqBase_.cnt_ += clusters[posItem.element].seqBase_.cnt_;
					//clusters[sub.second.front()].reads_.front().mateSeqBase_.cnt_ += clusters[posItem.element].mateSeqBase_.cnt_;
					removeThese.emplace_back(posItem.element);
				}
				if(qualRep == "max"){
					for(const auto & qualPos : iter::range(clusters[sub.second.front()].seqBase_.qual_.size())){
						std::vector<uint32_t> quals;
						quals.reserve(sub.second.size());
						for(const auto & pos : sub.second){
							quals.emplace_back(clusters[pos].seqBase_.qual_[qualPos]);
						}
						clusters[sub.second.front()].seqBase_.qual_[qualPos] = vectorMaximum(quals);
					}
				}else if (qualRep == "mean"){
					for(const auto & qualPos : iter::range(clusters[sub.second.front()].seqBase_.qual_.size())){
						std::vector<uint32_t> quals;
						quals.reserve(sub.second.size());
						for(const auto & pos : sub.second){
							quals.emplace_back(clusters[pos].seqBase_.qual_[qualPos]);
						}
						clusters[sub.second.front()].seqBase_.qual_[qualPos] = std::round(vectorMean(quals));
					}
				}else if (qualRep == "median"){
					for(const auto & qualPos : iter::range(clusters[sub.second.front()].seqBase_.qual_.size())){
						std::vector<uint32_t> quals;
						quals.reserve(sub.second.size());
						for(const auto & pos : sub.second){
							quals.emplace_back(clusters[pos].seqBase_.qual_[qualPos]);
						}
						clusters[sub.second.front()].seqBase_.qual_[qualPos] = vectorMedianRef(quals);
					}
				}else {
					std::cerr << bib::bashCT::red << bib::bashCT::bold
							<< "Unrecognized qualRep for collapseIdenticalPairedClusters\n"
							<< "should be mean, max, or median\n"
							<< "not " << qualRep << "\n";
					exit(1);
				}
			}
		}
	}
	bib::sort(removeThese);
	for(const auto & pos : iter::reverse(removeThese)){
		clusters.erase(clusters.begin() + pos);
	}
}

void PairedCluster::writeOutClusters(const std::string& directoryName,
		const SeqIOOptions & ioOptions) const {
	SeqOutput writer(ioOptions);
	writer.openWrite(reads_);
}

void PairedCluster::writeOutClustersWithConsensus(std::ostream & firstOut,std::ostream & secondOut )const{
	throw std::runtime_error{"PairedCluster::writeOutClustersWithConsensus not yet implemented"};
}


} /* namespace bibseq */
