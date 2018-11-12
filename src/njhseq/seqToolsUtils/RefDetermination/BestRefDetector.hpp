#pragma once
/*
 * BestRefDetector.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: nick
 */
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
#include "njhseq/concurrency/pools/AlignerPool.hpp"
#include "njhseq/objects/kmer/kmerInfo.hpp"

namespace njhseq {


class BestRefDetector {
public:

	struct FindBestRefPars {
		double kmerCutOff_ = 0.80;
		bool forceMatch_ = false;
		uint32_t kmerLen_ = 5;
		uint32_t numThreads = 1;
		bool scoreBased_ = false;
		uint32_t batchAmount_ = 1;
	};

	template<typename SEQTYPE,typename REFTYPE>
	static std::vector<comparison> findBestRef(
			const SEQTYPE & inputSeq,
			const std::vector<REFTYPE> & refSeqs,
			const std::vector<kmerInfo> & refInfos,
			aligner & alignerObj,
			const FindBestRefPars & pars){
		kmerInfo inputKInfo(getSeqBase(inputSeq).seq_, pars.kmerLen_, false);
    double bestScore = std::numeric_limits<double>::lowest();
    std::vector<uint32_t> bestRefs;
    double currentKmerCutOff = pars.kmerCutOff_;
    bool run = true;
    while(run){
	    for (const auto& refPos : iter::range(refSeqs.size())) {
	      const auto & ref = refSeqs[refPos];
	      if (getSeqBase(ref).name_ == getSeqBase(inputSeq).name_) {
	        continue;
	      }
	      if(refInfos[refPos].compareKmers(inputKInfo).second < currentKmerCutOff){
	       	continue;
	      }
				alignerObj.alignCacheGlobal(ref, inputSeq);
				double currentScore = 0;
				if(!pars.scoreBased_) {
					alignerObj.profileAlignment(ref, inputSeq, false, true, false);
					currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
				} else {
					currentScore = alignerObj.parts_.score_;
				}
				if (currentScore == bestScore) {
					bestRefs.push_back(refPos);
				}
				if (currentScore > bestScore) {
					bestRefs.clear();
					bestRefs.push_back(refPos);
					bestScore = currentScore;
				}
			}
	    run = false;
	    if(bestRefs.empty() && pars.forceMatch_ && currentKmerCutOff > 0){
	    		run = true;
	    }
	    currentKmerCutOff -= 0.1;
    }
		std::vector<comparison> comps;
		for (const auto& bestPos : bestRefs) {
			const auto & best = refSeqs[bestPos];
			alignerObj.alignCacheGlobal(best, inputSeq);
			alignerObj.profileAlignment(best, inputSeq, false, true, false);
			comps.emplace_back(alignerObj.comp_);
		}
		return comps;
	}

	template<typename SEQTYPE,typename REFTYPE>
	static std::unordered_map<std::string, std::vector<comparison>> findBestRef(
			const std::vector<SEQTYPE> & inputSeqs,
			const std::vector<REFTYPE> & refSeqs,
			concurrent::AlignerPool & alnPool,
			const FindBestRefPars & pars){


		std::unordered_map<std::string, std::vector<comparison>> ret;
	  //set up queue
	  std::vector<uint32_t> positions(inputSeqs.size());
	  njh::iota<uint32_t>(positions, 0);
	  //njh::reverse(positions);
	  njh::concurrent::LockableQueue<uint32_t> posQueue(positions);
	  std::mutex mut;


	  auto compareInput = [&alnPool,&inputSeqs,&refSeqs,&posQueue,&pars,&mut,&ret](){
	  	std::vector<uint32_t> subPositions;
	  	auto curAligner = alnPool.popAligner();
	  	std::vector<kmerInfo> refInfos;
	  	for (const auto& refPos : iter::range(refSeqs.size())){
	  		refInfos.emplace_back(getSeqBase(refSeqs[refPos]).seq_, pars.kmerLen_, false);
	  	}
	  	while(posQueue.getVals(subPositions, pars.batchAmount_)){
	  		std::unordered_map<uint32_t, std::vector<uint32_t>> bestRefsForPos;
			for(const auto pos : iter::reversed(subPositions)){
					const auto & input = inputSeqs[pos];
					kmerInfo inputKInfo(getSeqBase(input).seq_, pars.kmerLen_, false);
			    double bestScore = std::numeric_limits<double>::lowest();
			    std::vector<uint32_t> bestRefs;
			    double currentKmerCutOff = pars.kmerCutOff_;
			    bool run = true;
			    while(run){
				    for (const auto& refPos : iter::range(refSeqs.size())) {
				      const auto & ref = refSeqs[refPos];
				      if (getSeqBase(ref).name_ == getSeqBase(input).name_) {
				        continue;
				      }
				      if(refInfos[refPos].compareKmers(inputKInfo).second < currentKmerCutOff){
				       	continue;
				      }
							curAligner->alignCacheGlobal(ref, input);
							double currentScore = 0;
							if(!pars.scoreBased_) {
								curAligner->profileAlignment(ref, input, false, true, false);
								currentScore = curAligner->comp_.distances_.eventBasedIdentity_;
							} else {
								currentScore = curAligner->parts_.score_;
							}
							if (currentScore == bestScore) {
								bestRefs.push_back(refPos);
							}
							if (currentScore > bestScore) {
								bestRefs.clear();
								bestRefs.push_back(refPos);
								bestScore = currentScore;
							}
						}
				    bestRefsForPos[pos] = bestRefs;

				    run = false;
				    if(bestRefs.empty() && pars.forceMatch_ && currentKmerCutOff > 0){
				    		run = true;
				    }
				    currentKmerCutOff -= 0.1;
			    }
				}
				{
					std::lock_guard<std::mutex> lock(mut);
					for(const auto & bestRefs : bestRefsForPos) {
						const auto & input = inputSeqs[bestRefs.first];
						std::vector<comparison> comps;
						for (const auto& bestPos : bestRefs.second) {
							const auto & best = refSeqs[bestPos];
							curAligner->alignCacheGlobal(best, input);
							curAligner->profileAlignment(best, input, false, true, false);
							comps.emplace_back(curAligner->comp_);
						}
						ret[getSeqBase(input).name_] = comps;
					}
				}
	  		}
		};
		std::vector<std::thread> threads;
		for (uint32_t t = 0; t < pars.numThreads; ++t) {
			threads.emplace_back(std::thread(compareInput));
		}
		njh::concurrent::joinAllJoinableThreads(threads);
		return ret;
	}

};

}  // namespace njhseq
