#include "cluster.hpp"
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "njhseq/IO/SeqIO/SeqOutput.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp"
namespace njhseq {

cluster::cluster() : baseCluster() { rejected_ = false; }

cluster::cluster(const seqInfo &firstRead) : baseCluster(firstRead) {
  rejected_ = false;
  setLetterCount();
}






void cluster::addRead(const cluster& cr) {
  reads_.insert(reads_.end(), cr.reads_.begin(), cr.reads_.end());
  seqBase_.cnt_ += cr.seqBase_.cnt_;
  seqBase_.frac_ = seqBase_.frac_ * reads_.size();
  seqBase_.frac_ = (seqBase_.frac_ + cr.seqBase_.frac_) / reads_.size();
  // needToCalculateConsensus = true;
  if (seqBase_.cnt_ / 2.0 >= firstReadCount_) {
    needToCalculateConsensus_ = true;
  }
}

int cluster::getBiggestReadSize() {
  int biggestSize = 0;
  for (const auto& read : reads_) {
    if (read->seqBase_.cnt_ > biggestSize) {
      biggestSize = read->seqBase_.cnt_;
    }
  }
  return biggestSize;
}
double cluster::getAverageSizeDifference() {
  double averageSizeDifference = 0;
  for (const auto& read : reads_) {
    averageSizeDifference +=
        uAbsdiff(read->seqBase_.seq_.length(), seqBase_.seq_.length()) *
        read->seqBase_.cnt_;
  }
  averageSizeDifference = averageSizeDifference / seqBase_.cnt_;
  return averageSizeDifference;
}

size_t cluster::getLargestSizeDifference() {
  size_t largestSizeDifference = 0;
  for (const auto& read : reads_) {
    if (uAbsdiff(read->seqBase_.seq_.length(),seqBase_.seq_.length()) >
        largestSizeDifference) {
      largestSizeDifference =
      		uAbsdiff(read->seqBase_.seq_.length(), seqBase_.seq_.length());
    }
  }
  return largestSizeDifference;
}


void cluster::outputInfoComp(const std::string& workingDir) const {
  std::ofstream info(njh::files::make_path(workingDir, seqBase_.name_).string());
  if (!info) {
    std::cout << "Error in opening" << workingDir << seqBase_.name_
              << std::endl;
  }
  std::map<int, int> clusterCounter;
  for (const auto& rIter : reads_) {
    ++clusterCounter[rIter->seqBase_.cnt_];
  }
  info << "clusterSize\tFrequency" << std::endl;
  for (const auto& kv : clusterCounter) {
    info << kv.first << "\t" << kv.second << std::endl;
  }
  return;
}

std::vector<cluster> cluster::breakoutClustersBasedOnSnps(aligner & alignerObj,
		const snpBreakoutPars& pars) {
	std::vector<cluster>  ret;
	//log snp information
	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> mismatches;
	for (const auto & subReadPos : iter::range(reads_.size())) {
		const auto & subRead = reads_[subReadPos];
		alignerObj.alignCache(*this, subRead, false);
		//count gaps and mismatches and get identity
		alignerObj.profilePrimerAlignment(*this, subRead);
		for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
			if (m.second.highQualityJustSeq(pars.qScorePars)) {
				mismatches[m.second.refBasePos][m.second.seqBase] += subRead->seqBase_.cnt_;
			}
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "mismatches.size(): " << mismatches.size() << std::endl;
	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> mismatchesAboveCutOff;
	for (const auto & position : mismatches) {
		for (const auto & base : position.second) {
//			std::cout << position.first << "\t" << base.first << '\t' << base.second << std::endl;
			if (base.second > pars.hardCutOff) {
				mismatchesAboveCutOff[position.first][base.first] = base.second;
			}
		}
	}
//	std::cout << "mismatchesAboveCutOff.size(): " << mismatchesAboveCutOff.size() << std::endl;
	if (!mismatchesAboveCutOff.empty()) {
		std::unordered_map<std::string, std::vector<uint32_t>> readsSnpUids;
		for (const auto & subReadPos : iter::range(
				reads_.size())) {
			const auto & subRead = reads_[subReadPos];
			alignerObj.alignCache(*this, subRead, false);
			//count gaps and mismatches and get identity
			alignerObj.profilePrimerAlignment(*this, subRead);
			std::stringstream ss;
			uint32_t snpCount = 0;
			double freqSum = 0;
			for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
				if (njh::in(m.second.seqBase,
						mismatchesAboveCutOff[m.second.refBasePos])) {
					++snpCount;
					freqSum +=  mismatchesAboveCutOff[m.second.refBasePos][m.second.seqBase]/seqBase_.cnt_;
					ss << m.second.refBasePos << ":" << m.second.seqBase << ";";
				}
			}
			std::string snpProfileUid = ss.str();
			if ("" != snpProfileUid &&
					snpCount >= pars.minSnps &&
					freqSum >= pars.snpFreqCutOff) {
				readsSnpUids[snpProfileUid].emplace_back(subReadPos);
			}
		}
		std::vector<uint32_t> readsToErase;
		for (const auto & readsWithSnpUid : readsSnpUids) {
//			std::cout << readsWithSnpUid.first << " " << readsWithSnpUid.second.size() << std::endl;
			if (readsWithSnpUid.second.size() > pars.hardCutOff) {
//				std::cout << "\t" << readsWithSnpUid.first << " " << readsWithSnpUid.second.size() << std::endl;
				std::vector<std::shared_ptr<readObject>> splitSeqs;
				for(const auto & pos : readsWithSnpUid.second){
					splitSeqs.push_back(reads_[pos]);
					readsToErase.emplace_back(pos);
				}
				readVecSorter::sort(splitSeqs);
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				OutOptions outOpts(bfs::path("first_name_of_first_splitSeqs.txt"));
//				outOpts.append_ = true;
//				OutputStream out(outOpts);
//				out << splitSeqs.front()->seqBase_.name_ << roundDecPlaces(splitSeqs.front()->averageErrorRate, 4)<< std::endl;
				ret.emplace_back(getConsensus(splitSeqs, alignerObj, splitSeqs.front()->seqBase_.name_));
				ret.back().needToCalculateConsensus_ = true;
				ret.back().calculateConsensus(alignerObj, true);;
			}
		}
		if(!readsToErase.empty()){
			removeReads(readsToErase);
		}
		readVecSorter::sort(reads_);
		needToCalculateConsensus_ = true;
		calculateConsensus(alignerObj, true);
		//recalc
		needToCalculateConsensus_ = true;
		calculateConsensus(alignerObj, true);
	}
	/*
	std::unordered_map<uint32_t, std::vector<char>> mismatchesAboveCutOff;
	for (const auto & position : mismatches) {
		for (const auto & base : position.second) {
			if (base.second / seqBase_.cnt_ > pars.snpFreqCutOff
					&& base.second > pars.hardCutOff) {
				mismatchesAboveCutOff[position.first].emplace_back(base.first);
			}
		}
	}
	if (!mismatchesAboveCutOff.empty()) {
		std::unordered_map<std::string, std::vector<uint32_t>> readsSnpUids;
		for (const auto & subReadPos : iter::range(
				reads_.size())) {
			const auto & subRead = reads_[subReadPos];
			alignerObj.alignCache(*this, subRead, false);
			//count gaps and mismatches and get identity
			alignerObj.profilePrimerAlignment(*this, subRead);
			std::stringstream ss;
			uint32_t snpCount = 0;
			double freqSum = 0;
			for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
				if (njh::in(m.second.seqBase,
						mismatchesAboveCutOff[m.second.refBasePos])) {
					++snpCount;
					ss << m.second.refBasePos << ":" << m.second.seqBase << ";";
				}
			}
			std::string snpProfileUid = ss.str();
			if ("" != snpProfileUid &&
					snpCount >= pars.minSnps) {
				readsSnpUids[snpProfileUid].emplace_back(subReadPos);
			}
		}
		std::vector<uint32_t> readsToErase;
		for (const auto & readsWithSnpUid : readsSnpUids) {
			if (readsWithSnpUid.second.size() > pars.hardCutOff) {
				//std::cout << readsWithSnpUid.first << " " << readsWithSnpUid.second.size() << std::endl;
				std::vector<std::shared_ptr<readObject>> splitSeqs;
				for(const auto & pos : readsWithSnpUid.second){
					splitSeqs.push_back(reads_[pos]);
					readsToErase.emplace_back(pos);
				}
				readVecSorter::sort(splitSeqs);
				ret.emplace_back(getConsensus(splitSeqs, alignerObj, splitSeqs.front()->seqBase_.name_));
			}
		}
		if(!readsToErase.empty()){
			removeReads(readsToErase);
		}
	}*/
	return ret;
}



}  // namespace njh
