#include "cluster.hpp"

namespace bibseq {

std::unordered_map<std::string, uint32_t> cluster::mutateConensus(const std::vector<std::vector<uint32_t>>& currentAlignQuals,
		std::unordered_map<double, double> bestLikelihood,
		simulation::mismatchProfile & eProfile,
		bib::randomGenerator & gen){
	//mutate the number of reads are in the cluster
	std::unordered_map<std::string, uint32_t> mutated;
	std::vector<uint32_t> randomQualPos = gen.unifRandVector<uint32_t>(0,currentAlignQuals.size(), std::round(seqBase_.cnt_));
	for(const auto & randQual : randomQualPos ){
		std::vector<double> currentLikelihood = likelihoodForBaseQ(currentAlignQuals[randQual] , bestLikelihood);
		std::string currentMutant = eProfile.mutateSeq(seqBase_.seq_, gen, eProfile.alphabet_, currentLikelihood);
		++mutated[currentMutant];
	}
	return mutated;
}

std::unordered_map<std::string, uint32_t> cluster::mutateConensus(double errorRate,
		simulation::mismatchProfile & eProfile,
		bib::randomGenerator & gen){
	//mutate the number of reads are in the cluster
	std::unordered_map<std::string, uint32_t> mutated;
	for(uint32_t i = 0; i < seqBase_.cnt_; ++i){
		std::string currentMutant = eProfile.mutateSeqSameErrorRate(seqBase_.seq_, gen, eProfile.alphabet_, errorRate);
		++mutated[currentMutant];
	}
	return mutated;
}


std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<cluster>>> cluster::createSimReadClusters(aligner & alignerObj ){
  std::unordered_map<std::string, std::vector<readObject>> similarReads;
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<cluster>>> similarReadsClusters;
  //std::cout << "here0" << std::endl;
  //seqBase_.printDescription(std::cout, true);
  createCondensedSeq();
  //std::cout << "here1" << std::endl;
  for (auto& read : reads_) {
  	//std::cout << "here2" << std::endl;
  	read.createCondensedSeq();
  	//std::cout << "here3" << std::endl;
  	bool sameCondensed = false;
  	//std::cout << "here4" << std::endl;


  	//std::cout << "here5" << std::endl;
  	if(sameCondensed){
  		similarReads["none"].push_back(read);
  	}else{
    	//align all reads to consensus and get only the mismatch differences
    	//to count clusters without regards to indels

      alignerObj.alignVec(*this, read, false);
      alignerObj.profileAlignment(*this, read, 11, true, false,
                                  true, false, true);

      std::map<uint32_t, char> currentMismatches;
      for (const auto& mis : alignerObj.mismatches_) {
        currentMismatches[mis.second.refBasePos] = mis.second.seqBase;
      }
      std::stringstream currentStream;
      for (const auto& mis : currentMismatches) {
        currentStream << mis.first << "_" << mis.second;
      }
      if ("" == currentStream.str()) {
        similarReads["none"].push_back(read);
      } else {
        similarReads[currentStream.str()].emplace_back(read);
      }
  	}

  }

  //convert the similar read vectors into cluster objects
  for (const auto& simReads : similarReads) {
  	uint32_t mismatchCount =  static_cast<uint32_t>(
    std::count(simReads.first.begin(), simReads.first.end(), '_'));
  	uint32_t clusSize = readVec::getTotalReadCount(simReads.second);
    similarReadsClusters[mismatchCount][clusSize].emplace_back(createClusterFromSimiliarReads(simReads.second, alignerObj));
    if(similarReadsClusters[mismatchCount][clusSize].back().seqBase_.cnt_ < 1){
    	similarReadsClusters[mismatchCount][clusSize].back().printDescription(std::cout, true);
    	std::cout << std::endl;
    	bib::for_each(simReads.second, [&](const readObject & r){ r.printDescription(std::cout, true) ;});
    	exit(1);
    }
  }
  return similarReadsClusters;
}

void cluster::removeClustersOnFDR(aligner & alignerObj, double fdrCutOff, std::vector<cluster> & rejectedClusters){
	auto similarReadsClusters = createSimReadClusters(alignerObj);
	for(const auto & m : similarReadsClusters){
		for(const auto & c : m.second){
  		if(getFDRValue(c.first, m.first, c.second.size()) <= fdrCutOff){
  			for(const auto & rejected : c.second){
  				rejectedClusters.emplace_back(rejected);
  				rejectedClusters.back().rejected_ = true;
  				removeReads(rejected.reads_);
  			}
  		}
  	}
  }
	return;
}

void cluster::simOnQual(const std::vector<identicalCluster> &initialClusters,
		aligner & alignerObj, uint32_t runTimes,
		std::unordered_map<double, double> bestLikelihood,
		simulation::mismatchProfile & eProfile, bib::randomGenerator & gen) {

	std::vector<readObject> collection = getAllBeginingClusters(initialClusters);
	std::vector<std::vector<uint32_t>> currentAlignQuals = alignOrigQuals(
			collection, alignerObj, false);
	auto allRunsCounts = simulate(runTimes, currentAlignQuals, bestLikelihood,
			eProfile, gen);
	postProcessSimulation(runTimes, allRunsCounts);
	return;
}

void cluster::postProcessSimulation(uint32_t runTimes,
		const std::unordered_map<uint32_t,
				std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>>& allRunsCounts) {
	pValueMap_.clear();
  amountAverageMap_.clear();
	std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> presentMap;
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> amountMap;
  for(const auto & runCount : allRunsCounts){
  	std::map<uint32_t, uint32_t> maxes;
  	uint32_t lastMax = 0;
  	for(const auto & m : runCount.second){
  		maxes[m.first] = vectorMaximum(getVectorOfMapKeys(m.second));
  	}
  	for(const auto & mFirst : iter::reverse(maxes)){
  		for(const auto & mSecond : iter::range(lastMax + 1,mFirst.second + 1)){
  			for(const auto & mis : iter::range<uint32_t>(1, mFirst.first + 1)){
  				++presentMap[mis][mSecond];
  			}
  		}
  		if(mFirst.second > lastMax){
  			lastMax = mFirst.second;
  		}
  	}
  	for(const auto & mis : runCount.second){
  		if(mis.first == 0){
  			continue;
  		}
  		for(const auto & clus : mis.second){
  			for(auto cMis : iter::range<uint32_t>(1,mis.first +1)){
  				for(auto cClus : iter::range<uint32_t>(1,clus.first +1)){
  					amountMap[cMis][cClus]+=clus.second;
  				}
  			}
  		}
  	}
  }
  //convert the present map and amount maps into averaged over the run times
  for(const auto & mis : presentMap){
  	for(const auto & clus : mis.second){
  		pValueMap_[mis.first][clus.first] = clus.second/static_cast<double>(runTimes);
  	}
  }
  for(const auto & mis : amountMap){
  	for(const auto & clus : mis.second){
  		amountAverageMap_[mis.first][clus.first] = clus.second/static_cast<double>(runTimes);
  	}
  }
  return;
}
std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>> cluster::simulate(uint32_t runTimes,
		const std::vector<std::vector<uint32_t>>& currentAlignQuals,
		std::unordered_map<double, double> bestLikelihood,
		simulation::mismatchProfile & eProfile,
		bib::randomGenerator & gen){
	// k1 = run, k2 = number of mismatches, k3 = cluster number, v  = occurrences
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>> allRunsCounts;
  bib::scopedStopWatch mutation(seqBase_.name_ + "_mutagenesis", true);
  for(uint32_t run = 0; run < runTimes; ++run){
  	std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> currentRunCount;
  	std::unordered_map<std::string, uint32_t> mutated = mutateConensus(currentAlignQuals, bestLikelihood, eProfile, gen);
    for(const auto & mut : mutated){
    	++currentRunCount[numberOfMismatches(seqBase_.seq_, mut.first)][mut.second];
    }
    allRunsCounts[run] = currentRunCount;
  }
  return allRunsCounts;
}

double cluster::getAverageFrequencyOfClusterVector() const {
  double sum = 0.00;
  for (const auto& read : reads_) {
    sum += read.seqBase_.frac_;
  }
  return sum / reads_.size();
}

//fdr calc
double cluster::getFDRValue(uint32_t clusSize, uint32_t mismatches, uint32_t observedAmount){
	if(mismatches == 0){
		return 1;
	}
	return amountAverageMap_[mismatches][clusSize]/static_cast<double>(observedAmount);
}

//pvalue calc
double cluster::getPValue(const seqInfo & seqBase, uint32_t numOfMismatches){
	if(numOfMismatches == 0){
		return 1;
	}
	return pValueMap_[numOfMismatches][::round( seqBase.cnt_)];
}

double cluster::getPValue(const seqInfo & seqBase, aligner & alignerObj, bool local){
	alignerObj.alignVec(seqBase_, seqBase, local);
	alignerObj.profileAlignment(seqBase_, seqBase, 11, true, false, false, false, true);
	if(static_cast<int>(alignerObj.comp_.largeBaseIndel_) > 0){
		return 0;
	}
	return getPValue(seqBase, alignerObj.comp_.hqMismatches_);
}

void cluster::removeRead(const std::string & stubName){
	uint32_t readPos = std::numeric_limits<uint32_t>::max();
	for(const auto & pos : iter::range(reads_.size())){
		if(reads_[pos].getStubName(true) == stubName){
			readPos = pos;
			break;
		}
	}
	if(readPos!=std::numeric_limits<uint32_t>::max()){
		seqBase_.cnt_ -= reads_[readPos].seqBase_.cnt_;
		reads_.erase(reads_.begin() + readPos);
	  if (seqBase_.cnt_ / 2 > firstReadCount_) {
	    needToCalculateConsensus_ = true;
	  }
	}
}
void cluster::removeReads(const std::vector<readObject> & vec){
	for(const auto & read : vec){
		removeRead(read.getStubName(true));
	}
	updateName();
}
void cluster::addRead(const baseCluster& cr) {
  // clusterVector.push_back(cr);

  reads_.insert(reads_.end(), cr.reads_.begin(), cr.reads_.end());
  seqBase_.cnt_ += cr.seqBase_.cnt_;
  seqBase_.frac_ = seqBase_.frac_ * reads_.size();
  seqBase_.frac_ = (seqBase_.frac_ + cr.seqBase_.frac_) / reads_.size();
  // needToCalculateConsensus = true;
  if (seqBase_.cnt_ / 2 > firstReadCount_) {
    needToCalculateConsensus_ = true;
  }
}

int cluster::getBiggestReadSize() {
  int biggestSize = 0;
  for (const auto& read : reads_) {
    if (read.seqBase_.cnt_ > biggestSize) {
      biggestSize = read.seqBase_.cnt_;
    }
  }
  return biggestSize;
}
double cluster::getAverageSizeDifference() {
  double averageSizeDifference = 0;
  for (const auto& read : reads_) {
    averageSizeDifference +=
        uAbsdiff(read.seqBase_.seq_.length(), seqBase_.seq_.length()) *
        read.seqBase_.cnt_;
  }
  averageSizeDifference = averageSizeDifference / seqBase_.cnt_;
  return averageSizeDifference;
}

int cluster::getLargestSizeDifference() {
  int largestSizeDifference = 0;
  for (const auto& read : reads_) {
    if (uAbsdiff(read.seqBase_.seq_.length(),seqBase_.seq_.length()) >
        largestSizeDifference) {
      largestSizeDifference =
      		uAbsdiff(read.seqBase_.seq_.length(), seqBase_.seq_.length());
    }
  }
  return largestSizeDifference;
}


void cluster::outputInfoComp(const std::string& workingDir) const {
  std::ofstream info(combineStrings({workingDir, seqBase_.name_}).c_str());
  if (!info) {
    std::cout << "Error in opening" << workingDir << seqBase_.name_
              << std::endl;
  }
  std::map<int, int> clusterCounter;
  // gatherInfoAboutIdenticalReadCompostion(clusterCounter);
  for (const auto& rIter : reads_) {
    clusterCounter[rIter.seqBase_.cnt_]++;
  }
  info << "clusterSize\tFrequency" << std::endl;
  for (const auto& kv : clusterCounter) {
    info << kv.first << "\t" << kv.second << std::endl;
  }
  return;
}

std::vector<readObject> cluster::getAllBeginingClusters(
    const std::vector<identicalCluster>& initialClusters) const {
  std::vector<readObject> collection;
  collection.reserve(seqBase_.cnt_);
  for (const auto& read : reads_) {
    addOtherVec(collection, readVec::getReadByName(initialClusters,
                                                   read.seqBase_.name_).reads_);
  }
  return collection;
}

std::vector<std::vector<uint32_t>> cluster::alignOrigQuals(const std::vector<readObject> & reads,
		  aligner & alignerObj, bool local)const{
	std::vector<std::vector<uint32_t>> ans;
	std::map<std::string, std::vector<uint32_t>> positions;
	for(const auto & read : reads){
		ans.emplace_back(getQualToConsensus(read, alignerObj, local, positions));
	}
	return ans;
}
std::vector<uint32_t> cluster::getQualToConsensus(
    const readObject& read, aligner& alignerObj, bool local,
    std::map<std::string, std::vector<uint32_t>>& positions) const {
  if (read.seqBase_.seq_ == seqBase_.seq_) {
    return read.seqBase_.qual_;
  } else if (positions.find(read.seqBase_.seq_) == positions.end()) {
    alignerObj.alignVec(*this, read, local);
    auto pos =
        seqUtil::getQualPositions(alignerObj.alignObjectA_.seqBase_.seq_,
                                  alignerObj.alignObjectB_.seqBase_.seq_);
    positions[read.seqBase_.seq_] = pos;
    return seqUtil::rearrangeQuals(read.seqBase_.qual_, pos);
  } else {
    return seqUtil::rearrangeQuals(read.seqBase_.qual_,
                                   positions[read.seqBase_.seq_]);
  }
}

void cluster::writeOutInputClusters(const std::string& workingDir) const {
	readObjectIOOptions options;
	options.outFilename_ = workingDir + seqBase_.name_;
	options.outFormat_ = "fastaQual";
  readObjectIO::write(allInputClusters,options);
}
}  // namespace bib
