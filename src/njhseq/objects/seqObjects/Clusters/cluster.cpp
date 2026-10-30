#include "cluster.hpp"
#include <boost/math/distributions/binomial.hpp>
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

bool cluster::snpPassesStrandComposition(double snpReadCount, double snpCompReadCount,
		double compReadsInCluster, double readsInCluster, double pValueCutOff) {
	//if the entire cluster is one direction or the other (which might mean the extraction was in one direction anyways) then don't filter based on comp
	if (compReadsInCluster <= 0 || compReadsInCluster >= readsInCluster) {
		return true;
	}
	const auto trials = static_cast<uint64_t>(std::round(snpReadCount));
	if (0 == trials) {
		return false;
	}
	const auto compSuccesses = std::min(trials, static_cast<uint64_t>(std::round(snpCompReadCount)));
	//two sided exact binomial test, the summed probability of every outcome no more likely than the observed one
	boost::math::binomial_distribution<double> dist(static_cast<double>(trials), compReadsInCluster / readsInCluster);
	const double observedProb = boost::math::pdf(dist, static_cast<double>(compSuccesses));
	const double relErr = 1 + 1e-7;
	double pValue = 0;
	for (uint64_t successes = 0; successes <= trials; ++successes) {
		const double prob = boost::math::pdf(dist, static_cast<double>(successes));
		if (prob <= observedProb * relErr) {
			pValue += prob;
		}
	}
	return std::min(1.0, pValue) >= pValueCutOff;
}

double cluster::snpBackgroundErrorRate(std::vector<double> mismatchRates, uint64_t possibleMismatches,
		double quantile, double minErrorRate) {
	//mismatchRates holds only the observed mismatches, the rest of the possible mismatches have a rate of 0
	double rate = 0;
	if (possibleMismatches > 0 && !mismatchRates.empty()) {
		std::sort(mismatchRates.begin(), mismatchRates.end());
		const uint64_t total = std::max<uint64_t>(possibleMismatches, mismatchRates.size());
		const uint64_t zeros = total - mismatchRates.size();
		const auto quantilePos = static_cast<uint64_t>(std::floor(std::clamp(quantile, 0.0, 1.0) * static_cast<double>(total - 1)));
		if (quantilePos >= zeros) {
			rate = mismatchRates[quantilePos - zeros];
		}
	}
	return std::max(rate, minErrorRate);
}

bool cluster::snpAboveBackgroundError(double snpReadCount, double readCount, double errorRate, double pValueCutOff) {
	const auto trials = static_cast<uint64_t>(std::round(readCount));
	const auto successes = std::min(trials, static_cast<uint64_t>(std::round(snpReadCount)));
	if (0 == trials || 0 == successes) {
		return false;
	}
	if (errorRate >= 1) {
		return false;
	}
	//upper tail, the chance of seeing at least this many reads with the snp from errors alone
	boost::math::binomial_distribution<double> dist(static_cast<double>(trials), errorRate);
	const double pValue = boost::math::cdf(boost::math::complement(dist, static_cast<double>(successes - 1)));
	return pValue < pValueCutOff;
}

bool cluster::snpPassesBreakoutCutOffs(double snpReadCount, double readCount, const snpBreakoutPars & pars,
		double backgroundErrorRate, double errorModelPValueCutOff) {
	const bool passCounts = snpReadCount > pars.hardCutOff &&
			snpReadCount / readCount > pars.hardSnpFreqCutOff;
	if (!pars.useErrorModel) {
		return passCounts;
	}
	const bool passErrorModel = snpReadCount >= pars.errorModelMinCount &&
			snpAboveBackgroundError(snpReadCount, readCount, backgroundErrorRate, errorModelPValueCutOff);
	return pars.errorModelAsFilter ? passCounts && passErrorModel : passErrorModel;
}

bool cluster::readIsRevComp(const std::string & name) {
	std::string_view nameView(name);
	//strip one trailing count (_t) or fraction (_f) that collapsing can add after the extractor's _Comp, numbers can be in scientific notation (e.g. _t1e+06)
	const auto lastUnderscore = nameView.rfind('_');
	if (std::string_view::npos != lastUnderscore && lastUnderscore + 2 < nameView.size()
			&& ('t' == nameView[lastUnderscore + 1] || 'f' == nameView[lastUnderscore + 1])) {
		const auto number = nameView.substr(lastUnderscore + 2);
		const bool allNumberChars = std::all_of(number.begin(), number.end(), [](char c) {
			return std::isdigit(static_cast<unsigned char>(c)) || '.' == c || 'e' == c || 'E' == c || '+' == c || '-' == c;
		});
		const bool hasDigit = std::any_of(number.begin(), number.end(), [](char c) {
			return std::isdigit(static_cast<unsigned char>(c));
		});
		if (allNumberChars && hasDigit) {
			nameView = nameView.substr(0, lastUnderscore);
		}
	}
	const std::string_view compMark("_Comp");
	return nameView.size() >= compMark.size() && nameView.substr(nameView.size() - compMark.size()) == compMark;
}

std::vector<cluster> cluster::breakoutClustersBasedOnSnps(aligner & alignerObj,
		const snpBreakoutPars& pars) {
	std::vector<cluster>  ret;
	//log snp information
	// std::cout << seqBase_.name_ << std::endl;
	std::unordered_map<uint32_t, std::unordered_map<char, double>> mismatches;
	//the strand test uses every read with the snp regardless of quality since base quality can differ by strand and would make a real snp look strand biased
	std::unordered_map<uint32_t, std::unordered_map<char, double>> allQualMismatches;
	std::unordered_map<uint32_t, std::unordered_map<char, double>> allQualMismatches_comp;
	for (const auto subReadPos : iter::range(reads_.size())) {
		const auto & subRead = reads_[subReadPos];
		alignerObj.alignCache(*this, subRead, false);
		//count gaps and mismatches and get identity
		alignerObj.profilePrimerAlignment(*this, subRead);
		const bool subReadIsRevComp = readIsRevComp(subRead->seqBase_.name_);
		for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
			allQualMismatches[m.second.refBasePos][m.second.seqBase] += subRead->seqBase_.cnt_;
			if (subReadIsRevComp) {
				allQualMismatches_comp[m.second.refBasePos][m.second.seqBase] += subRead->seqBase_.cnt_;
			}
			if (m.second.highQualityJustSeq(pars.qScorePars)) {
				mismatches[m.second.refBasePos][m.second.seqBase] += subRead->seqBase_.cnt_;
			}
		}
	}
	//read counts rather than sub read counts so they're on the same scale as the snp counts
	double compReadCount = 0;
	double readCount = 0;
	for (const auto & subRead : reads_) {
		readCount += subRead->seqBase_.cnt_;
		if (readIsRevComp(subRead->seqBase_.name_)) {
			compReadCount += subRead->seqBase_.cnt_;
		}
	}
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	// std::cout << "pars.hardCutOff: " << pars.hardCutOff << std::endl;
	// std::cout << "pars.hardSnpFreqCutOff: " << pars.hardSnpFreqCutOff << std::endl;
	// std::cout << "pars.snpFreqCutOff: " << pars.snpFreqCutOff << std::endl;
	// std::cout << "mismatches.size(): " << mismatches.size() << std::endl;

	//with the error model, the background error rate is a high quantile of every position's and alternative base's mismatch rate (most are 0)
	double backgroundErrorRate = 0;
	double errorModelPValueCutOff = 0;
	if (pars.useErrorModel) {
		std::vector<double> mismatchRates;
		for (const auto & position : mismatches) {
			for (const auto & base : position.second) {
				mismatchRates.emplace_back(base.second / seqBase_.cnt_);
			}
		}
		const uint64_t possibleMismatches = 3 * len(seqBase_);
		backgroundErrorRate = snpBackgroundErrorRate(mismatchRates, possibleMismatches,
				pars.errorModelQuantile, pars.errorModelMinErrorRate);
		//bonferroni correction for testing every position and alternative base
		errorModelPValueCutOff = pars.errorModelPValueCutOff / static_cast<double>(std::max<uint64_t>(1, possibleMismatches));
	}

	std::unordered_map<uint32_t, std::unordered_map<char, double>> mismatchesAboveCutOff;
	for (const auto & position : mismatches) {
		for (const auto & base : position.second) {

			const double snpAllQualCount = allQualMismatches[position.first][base.first];
			double snpCompCount = 0;
			auto compPosSearch = allQualMismatches_comp.find(position.first);
			if (allQualMismatches_comp.end() != compPosSearch) {
				auto compBaseSearch = compPosSearch->second.find(base.first);
				if (compPosSearch->second.end() != compBaseSearch) {
					snpCompCount = compBaseSearch->second;
				}
			}
			bool pass_comp = snpPassesStrandComposition(snpAllQualCount, snpCompCount,
					compReadCount, readCount, pars.strandBiasPValueCutOff);
			const bool passCountAndFreq = snpPassesBreakoutCutOffs(base.second, seqBase_.cnt_, pars,
					backgroundErrorRate, errorModelPValueCutOff);
			if (passCountAndFreq && pass_comp) {
				// std::cout << position.first << "\t" << base.first << '\t' << base.second << '\t' << base.second/seqBase_.cnt_<< std::endl;
				mismatchesAboveCutOff[position.first][base.first] = base.second;
			}
		}
	}
	// std::cout << "mismatchesAboveCutOff.size(): " << mismatchesAboveCutOff.size() << std::endl;

	if (!mismatchesAboveCutOff.empty()) {
		std::unordered_map<std::string, std::vector<uint32_t>> readsSnpUids;
		for (const auto subReadPos : iter::range(reads_.size())) {
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

			if (!snpProfileUid.empty() &&
					snpCount >= pars.minSnps
//          && freqSum >= pars.snpFreqCutOff
          ) {
				// std::cout << '\t' << "snpCount:      " << snpCount << std::endl;
				// std::cout << '\t' << "freqSum:       " << freqSum << std::endl;
				// std::cout << "\t" << "snpProfileUid: " << snpProfileUid << std::endl;
				readsSnpUids[snpProfileUid].emplace_back(subReadPos);
			}
		}
		std::vector<uint32_t> readsToErase;
		// std::cout << "Counts for reads with high freq snps: " << std::endl;
		for (const auto & readsWithSnpUid : readsSnpUids) {
      double countOfReadsWithSnpUID = 0;

      for(const auto & seqPos : readsWithSnpUid.second){
        countOfReadsWithSnpUID += reads_[seqPos]->seqBase_.cnt_;
      }
			double freqOfReadsWithSnpUID = countOfReadsWithSnpUID/seqBase_.cnt_;
			// std::cout << "\t" << readsWithSnpUid.first << " " << readsWithSnpUid.second.size() << " " << freqOfReadsWithSnpUID << std::endl;

			if (countOfReadsWithSnpUID > pars.hardCutOff && freqOfReadsWithSnpUID > pars.clusterFreqCutOff) {
//        if (readsWithSnpUid.second.size() > pars.hardCutOff) {
				// std::cout << njh::bashCT::red;
				// std::cout << "\t" << readsWithSnpUid.first << " " << readsWithSnpUid.second.size() << " " << freqOfReadsWithSnpUID << std::endl;
				// std::cout << njh::bashCT::reset;
				std::vector<std::shared_ptr<readObject>> splitSeqs;
				for(const auto & pos : readsWithSnpUid.second){
					splitSeqs.push_back(reads_[pos]);
					readsToErase.emplace_back(pos);
				}
				readVecSorter::sort(splitSeqs);
				// OutOptions outOpts(bfs::path("first_name_of_first_splitSeqs.txt"));
				// outOpts.append_ = true;
				// OutputStream out(outOpts);
				// out << splitSeqs.front()->seqBase_.name_ << roundDecPlaces(splitSeqs.front()->averageErrorRate, 4)<< std::endl;
				ret.emplace_back(getConsensus(splitSeqs, alignerObj, splitSeqs.front()->seqBase_.name_));
				ret.back().needToCalculateConsensus_ = true;
				ret.back().calculateConsensus(alignerObj, true);;
			}
		}
		// std::cout << "readsToErase: " << readsToErase.size()	<< std::endl;
		// std::cout << __FILE__ << " " << __LINE__ << std::endl << std::endl;

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

Json::Value cluster::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["super"] = baseCluster::toJson();
	ret["rejected_"] = njh::json::toJson(rejected_);
	ret["chimeras"] = njh::json::toJson(chimeras);
	ret["endChimeras"] = njh::json::toJson(endChimeras);
	ret["allInputClusters"] = njh::json::toJson(allInputClusters);
	return ret;
}

}  // namespace njh
