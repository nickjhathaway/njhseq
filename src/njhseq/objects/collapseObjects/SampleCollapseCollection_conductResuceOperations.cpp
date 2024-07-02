/*
 * SampleCollapseCollection_conductResuceOperations.cpp
 *
 *  Created on: Mar 15, 2020
 *      Author: nick
 */



#include "SampleCollapseCollection.hpp"
#include "njhseq/objects/seqObjects/Clusters/clusterUtils.hpp"


namespace njhseq {

namespace collapse {


bool SampleCollapseCollection::performLowLevelFilters(const performLowLevelFiltersPars & filtPars, aligner & alignerObj,
		const collapser & collapserObj, const CollapseIterations & popColIters) {
	bool anyFilterPerformed = false;
	if(filtPars.removeCommonlyLowFreqHaplotypes_){
		while(excludeCommonlyLowFreqHaps(filtPars.lowFreqHaplotypeFracCutOff_)){
			anyFilterPerformed = true;
			//if excluded run pop clustering again
			doPopulationClustering(createPopInput(), alignerObj, collapserObj, popColIters);
		}
	}

	if(filtPars.removeOneSampOnlyOneOffHaps_){
		if(excludeOneSampOnlyOneOffHaps(filtPars.oneSampOnlyOneOffHapsFrac_, alignerObj)){
			anyFilterPerformed = true;
			//if excluded run pop clustering again
			doPopulationClustering(createPopInput(),alignerObj, collapserObj, popColIters);
		}
	}

	if(filtPars.removeOneSampOnlyHaps_){
		if(excludeOneSampOnlyHaps(filtPars.oneSampOnlyHapsFrac_)){
			anyFilterPerformed = true;
			//if excluded run pop clustering again
			doPopulationClustering(createPopInput(), alignerObj, collapserObj, popColIters);
		}
	}
	return anyFilterPerformed;
}



bool SampleCollapseCollection::conductResuceOperations(const conductResuceOperationsPars & pars, aligner & alignerObj,
		const collapser & collapserObj, const CollapseIterations & popColIters){
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	//first gather major haplotypes
	std::set<std::string> majorHaps;
	std::set<std::string> majorHapsForChi;
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	for(const auto & sampleName : passingSamples_){
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		if(njh::in(sampleName, popNames_.controlSamples_)){
			continue;
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampleName);
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		auto sampPtr = sampleCollapses_.at(sampleName);
		for(uint32_t clusPos = 0;  clusPos < sampPtr->collapsed_.clusters_.size(); ++clusPos){
			if(sampPtr->collapsed_.clusters_[clusPos].seqBase_.frac_ >= pars.majorHaplotypeFracForRescue_){
				majorHaps.emplace(popCollapse_->collapsed_.clusters_[popCollapse_->collapsed_.subClustersPositions_.at(sampPtr->collapsed_.clusters_[clusPos].getStubName(true))].seqBase_.name_);
				if(clusPos < 2){
					majorHapsForChi.emplace(popCollapse_->collapsed_.clusters_[popCollapse_->collapsed_.subClustersPositions_.at(sampPtr->collapsed_.clusters_[clusPos].getStubName(true))].seqBase_.name_);
				}
				// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			}
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		if(!keepSampleInfoInMemory_){
			dumpSample(sampleName);
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	}
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if(pars.debug_){
		std::cout << "majorHaps: " << njh::conToStr(majorHaps, ",") << std::endl;
	}
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	bool rescuedHaplotypes = false;
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	for(const auto & sampleName : passingSamples_){
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampleName);
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		auto sampPtr = sampleCollapses_.at(sampleName);
		std::vector<uint32_t> toBeRescued;
		//iterator over haplotypes, determine if they should be considered for rescue, if they should be then check to see if they match a major haplotype
		for(const auto excludedPos : iter::range(sampPtr->excluded_.clusters_.size())){
			const auto & excluded = sampPtr->excluded_.clusters_[excludedPos];
			if(excluded.nameHasMetaData()){
				MetaDataInName excludedMeta(excluded.seqBase_.name_);
				std::set<std::string> otherExcludedCriteria;
				bool chimeriaExcludedRescue = false;
				bool oneOffExcludedRescue = false;
				bool lowFreqExcludedRescue = false;
				// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				for(const auto & excMeta : excludedMeta.meta_){
					if(njh::beginsWith(excMeta.first, "Exclude") ){
						bool other = true;
						if(pars.rescueExcludedChimericHaplotypes && "ExcludeIsChimeric" == excMeta.first){
							chimeriaExcludedRescue = true;
							other = false;
						}
						if(pars.rescueExcludedOneOffLowFreqHaplotypes && "ExcludeFailedLowFreqOneOff" == excMeta.first){
							oneOffExcludedRescue = true;
							other = false;
						}
						if(pars.rescueExcludedLowFreqHaplotypes && "ExcludeFailedFracCutOff" == excMeta.first){
							lowFreqExcludedRescue = true;
							other = false;
						}
						if(other){
							otherExcludedCriteria.emplace(excMeta.first);
						}
					}
				}
				// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				//check if it should be considered for rescue
				// std::cout << excluded.seqBase_.name_ << " consider for rescue: " << njh::colorBool((chimeriaExcludedRescue || oneOffExcludedRescue) && otherExcludedCriteria.empty()) << std::endl;
				// std::cout << excluded.seqBase_.name_ << " consider for rescue: " << njh::colorBool((chimeriaExcludedRescue || oneOffExcludedRescue || lowFreqExcludedRescue) && otherExcludedCriteria.empty()) << std::endl;
				// std::cout << excluded.seqBase_.name_ << " \tchimeriaExcludedRescue: " << njh::colorBool(chimeriaExcludedRescue) << std::endl;
				// std::cout << excluded.seqBase_.name_ << " \toneOffExcludedRescue: " << njh::colorBool(oneOffExcludedRescue) << std::endl;
				// std::cout << excluded.seqBase_.name_ << " \tlowFreqExcludedRescue: " << njh::colorBool(lowFreqExcludedRescue) << std::endl;
				// std::cout << excluded.seqBase_.name_ << " \totherExcludedCriteria.empty(): " << njh::colorBool(otherExcludedCriteria.empty()) << std::endl;

				if((chimeriaExcludedRescue || oneOffExcludedRescue || lowFreqExcludedRescue) && otherExcludedCriteria.empty()){
					//see if it matches a major haplotype
					bool rescue = false;
					if(chimeriaExcludedRescue){
						for(const auto & popHap : popCollapse_->collapsed_.clusters_){
							if(popHap.seqBase_.seq_ == excluded.seqBase_.seq_ &&
									popHap.seqBase_.cnt_ > excluded.seqBase_.cnt_ &&
									njh::in(popHap.seqBase_.name_, majorHapsForChi)){
								rescue = true;
								break;
							}
						}
					}else{
						for(const auto & popHap : popCollapse_->collapsed_.clusters_){
							if(popHap.seqBase_.seq_ == excluded.seqBase_.seq_ &&
									popHap.seqBase_.cnt_ > excluded.seqBase_.cnt_ &&
									njh::in(popHap.seqBase_.name_, majorHaps)){
								rescue = true;
								break;
							}
						}
					}
					if(rescue){
						toBeRescued.emplace_back(excludedPos);
					}
				}
				// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			}
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		if(!toBeRescued.empty()){
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			rescuedHaplotypes = true;
			std::sort(toBeRescued.rbegin(), toBeRescued.rend());
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			// std::cout << "toBeRescued.size() : " << toBeRescued.size() << std::endl;
			// for(const auto & toRescue : toBeRescued) {
			// 	std::cout << "\t" << "sampleName: " << sampleName << std::endl;
			// 	std::cout << "\t\t" << sampPtr->excluded_.clusters_[toRescue].seqBase_.name_ << std::endl;
			//
			// }
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			for(const auto toRescue : toBeRescued){
				MetaDataInName excludedMeta(sampPtr->excluded_.clusters_[toRescue].seqBase_.name_);
				excludedMeta.addMeta("rescue", "TRUE");
				excludedMeta.resetMetaInName(sampPtr->excluded_.clusters_[toRescue].seqBase_.name_);

				if(njh::in(std::string("ExcludeIsChimeric"), excludedMeta.meta_)) {
					//unmarking so as not to mess up chimera numbers
					sampPtr->excluded_.clusters_[toRescue].seqBase_.unmarkAsChimeric();
					for (auto & subRead : sampPtr->excluded_.clusters_[toRescue].reads_) {
						subRead->seqBase_.unmarkAsChimeric();
					}
				}
				sampPtr->collapsed_.clusters_.emplace_back(sampPtr->excluded_.clusters_[toRescue]);
				sampPtr->excluded_.clusters_.erase(sampPtr->excluded_.clusters_.begin() + toRescue);
			}
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			sampPtr->updateAfterExclustion();
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			sampPtr->renameClusters("fraction");
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		if(!keepSampleInfoInMemory_){
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			dumpSample(sampleName);
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		}
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	}
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if(rescuedHaplotypes){
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		//if excluded run pop clustering again
		doPopulationClustering(createPopInput(),
				alignerObj, collapserObj, popColIters);
		// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		return true;
	}
	// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	return false;
}


}  // namespace collapse

}  // namespace njhseq


