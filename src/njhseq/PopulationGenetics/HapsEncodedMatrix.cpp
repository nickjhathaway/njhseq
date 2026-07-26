/*
 * HapsEncodedMatrix.cpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nick
 */

#include "HapsEncodedMatrix.hpp"


namespace njhseq {




void HapsEncodedMatrix::SetWithExternalPars::setDefaults(seqSetUp & setUp){
  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);

  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);

  setUp.setOption(exclude_targets, "--exclude_targets", "Exclude these targets from analysis");
  setUp.setOption(exclude_samples, "--exclude_samples", "Exclude these samples from analysis");
  setUp.setOption(selectTargets, "--select_targets", "Analyze only these targets");
  setUp.setOption(selectSamples, "--select_samples", "Analyze only these samples");

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
	setUp.setOption(majorOnly, "--calcMajorHapOnly", "calculate differences by major haplotype only");
  setUp.setOption(minNumOfTargets, "--minNumOfTargets", "min number of targets per sample");

}


HapsEncodedMatrix::HapsEncodedMatrix(const SetWithExternalPars & pars): pars_(pars){
	TableReader hapTab(TableIOOpts(InOptions(pars.tableFnp), "\t", true));
	hapTab.header_.checkForColumnsThrow(VecStr{pars.sampleCol, pars.targetNameCol, pars.popIDCol, pars.relAbundCol}, __PRETTY_FUNCTION__);


	//set all info
	{
		//read in first to gather information on the table
		VecStr row;
		while(hapTab.getNextRow(row)){

			const auto &  samp = row[hapTab.header_.getColPos(pars.sampleCol)];
			const auto &  tar = row[hapTab.header_.getColPos(pars.targetNameCol)];
			if(!pars.selectSamples.empty() && !njh::in(samp, pars.selectSamples)){
				continue;
			}
		  if(!pars.exclude_samples.empty() && njh::in(samp, pars.exclude_samples)){
		    continue;
		  }
			if(!pars.selectTargets.empty() && !njh::in(tar, pars.selectTargets)){
				continue;
			}
		  if(!pars.exclude_targets.empty() && njh::in(tar, pars.exclude_targets)){
		    continue;
		  }
			const auto & hapName = row[hapTab.header_.getColPos(pars.popIDCol)];
			addSampTarHapForEncoding(samp, tar, hapName);
		}
	}
	//set encoding
	setEncodeKeys();
	//re-read and encode
	{
		TableReader reReadHapTab(TableIOOpts(InOptions(pars.tableFnp), "\t", true));
		VecStr row;
		while(reReadHapTab.getNextRow(row)){
			const auto &  samp = row[reReadHapTab.header_.getColPos(pars.sampleCol)];
			const auto &  tar = row[reReadHapTab.header_.getColPos(pars.targetNameCol)];
			if(!pars.selectSamples.empty() && !njh::in(samp, pars.selectSamples)){
				continue;
			}
			if(!pars.selectTargets.empty() && !njh::in(tar, pars.selectTargets)){
				continue;
			}
			const auto &  hapName = row[reReadHapTab.header_.getColPos(pars.popIDCol)];
			//auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(pars.relAbundCol)]); //doing nothing right now with this
			encodeSampTarHap(samp, tar, hapName);
		}
	}

  if(std::numeric_limits<uint32_t>::max() != pars.minNumOfTargets){
    auto numTargetsPerSample = getNumberTargetsPerSample();
    std::unordered_set<std::string> filteredSamples;
    for (const auto &count: numTargetsPerSample) {
      if (count.second >= pars_.minNumOfTargets) {
        filteredSamples.emplace(count.first);
      }
    }
    if(filteredSamples.empty()){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__  << " " << __LINE__ << "\n";
      ss << "Error, no filters above the sample min count: " << pars_.minNumOfTargets << "\n";
      throw std::runtime_error { ss.str() };
    }
    pars_.selectSamples = filteredSamples;
    resetEncoding();
    //set all info
    {
      //read in first to gather information on the table
      VecStr row;
      while(hapTab.getNextRow(row)){
        const auto &  samp = row[hapTab.header_.getColPos(pars_.sampleCol)];
        const auto &  tar = row[hapTab.header_.getColPos(pars_.targetNameCol)];
        if(!pars_.selectSamples.empty() && !njh::in(samp, pars_.selectSamples)){
          continue;
        }
        if(!pars_.exclude_samples.empty() && njh::in(samp, pars_.exclude_samples)){
          continue;
        }
        if(!pars_.selectTargets.empty() && !njh::in(tar, pars_.selectTargets)){
          continue;
        }
        if(!pars_.exclude_targets.empty() && njh::in(tar, pars_.exclude_targets)){
          continue;
        }
        const auto &  hapName = row[hapTab.header_.getColPos(pars_.popIDCol)];
        addSampTarHapForEncoding(samp, tar, hapName);
      }
    }
    //set encoding
    setEncodeKeys();
    //re-read and encode
    {
      TableReader reReadHapTab(TableIOOpts(InOptions(pars_.tableFnp), "\t", true));
      VecStr row;
      while(reReadHapTab.getNextRow(row)){
        const auto &  samp = row[reReadHapTab.header_.getColPos(pars_.sampleCol)];
        const auto &  tar = row[reReadHapTab.header_.getColPos(pars_.targetNameCol)];
        if(!pars_.selectSamples.empty() && !njh::in(samp, pars_.selectSamples)){
          continue;
        }
        if(!pars_.exclude_samples.empty() && njh::in(samp, pars_.exclude_samples)){
          continue;
        }
        if(!pars_.selectTargets.empty() && !njh::in(tar, pars_.selectTargets)){
          continue;
        }
        if(!pars_.exclude_targets.empty() && njh::in(tar, pars_.exclude_targets)){
          continue;
        }
        const auto &  hapName = row[reReadHapTab.header_.getColPos(pars_.popIDCol)];
        //auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(pars_.relAbundCol)]); //doing nothing right now with this
        encodeSampTarHap(samp, tar, hapName);
      }
    }

  }

	if(pars_.majorOnly){
	  add_relative_abundance();
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//now determine the major hap per sample per target
		for(const auto sampIndex : iter::range(sampNamesVec_.size())){
			for(const auto targetIndex : iter::range(tarNamesVec_.size())){
				//check if sample has target
				if(targetsEncodeBySamp_[sampIndex][targetIndex] > 0){
					auto maxAbundance = std::numeric_limits<double>::min();
					uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
					//determine the highest rel abundance
					for(const auto hapIndex : iter::range(tarStart_[targetIndex], tarStart_[targetIndex] + numberOfHapsPerTarget_[targetIndex]) ){
//							if(hapsEncodeBySampRelAbund_[sampIndex][hapIndex] > 0 && hapsEncodeBySampRelAbund_[sampIndex][hapIndex] < 1){
//								std::cout << hapsEncodeBySampRelAbund_[sampIndex][hapIndex] << std::endl;
//							}
						if(hapsEncodeBySampRelAbund_[sampIndex][hapIndex] > maxAbundance){
							maxAbundance = hapsEncodeBySampRelAbund_[sampIndex][hapIndex];
							bestIndex = hapIndex;
						}
					}
					//set all other indexes for the sample target hap presence to be 0
					for(const auto hapIndex : iter::range(tarStart_[targetIndex], tarStart_[targetIndex] + numberOfHapsPerTarget_[targetIndex]) ){
						if(hapIndex != bestIndex){
							hapsEncodeBySamp_[sampIndex][hapIndex] = 0;
						}
					}
				}
			}
		}
	}
}

HapsEncodedMatrix::SampTarHap::SampTarHap(const std::string &samp, const std::string &tar,
		const std::string &hap) :
		samp_(samp), tar_(tar), hap_(hap) {

}


void HapsEncodedMatrix::addSampTarHapForEncoding(const std::string &samp, const std::string &tar,
		const std::string &hap){
	sampNames_.emplace(samp);
	tarNames_.emplace(tar);
	++hapNamesForTars_[tar][hap];
}

void HapsEncodedMatrix::addSampTarHapForEncoding(const SampTarHap & adding){
	addSampTarHapForEncoding(adding.samp_, adding.tar_, adding.hap_);
}

void HapsEncodedMatrix::encodeSampTarHap(const std::string &samp, const std::string &tar,
		const std::string &hap){
	auto tKey = tarNameKey_[tar];
	auto hKey = hapNamesKey_[tar][hap];
//		std::cout << "sampNamesKey[samp]                        : " << sampNamesKey[samp] << std::endl;
//		std::cout << "tKey                                      : " << tKey << std::endl;
//		std::cout << "[tarStart[tKey]                           : " << tarStart[tKey] << std::endl;
//		std::cout << "hKey                                      : " << hKey << std::endl;
//		std::cout << "hKey                                      : " << hKey << std::endl;
//		std::cout << "hapsEncodeBySamp[sampNamesKey[samp]]size(): " << hapsEncodeBySamp[sampNamesKey[samp]].size() << std::endl;
	hapsEncodeBySamp_[sampNamesKey_[samp]][tarStart_[tKey] + hKey] = 1;
	targetsEncodeBySamp_[sampNamesKey_[samp]][tKey] = 1;
}

void HapsEncodedMatrix::encodeSampTarHap(const SampTarHap & adding){
	encodeSampTarHap(adding.samp_, adding.tar_, adding.hap_);
}


void HapsEncodedMatrix::calcHapProbs(){
	hapsProbs_ = std::vector<double>(totalHaps_, 0);
	for(uint32_t tarPos : iter::range(tarStart_.size())){
		std::vector<uint32_t> hapCounts(numberOfHapsPerTarget_[tarPos], 0);
		uint32_t totalHaps = 0;
		for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tarPos])){
			for(const auto sampPos : iter::range(hapsEncodeBySamp_.size())){
				hapCounts[hapPos] += hapsEncodeBySamp_[sampPos][tarStart_[tarPos] + hapPos];
				totalHaps += hapsEncodeBySamp_[sampPos][tarStart_[tarPos] + hapPos];
			}
		}
		for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tarPos])){
			hapsProbs_[tarStart_[tarPos] + hapPos] = hapCounts[hapPos]/static_cast<double>(totalHaps);
		}
	}
}

void HapsEncodedMatrix::add_relative_abundance() {
  //first fill the relative abundance vector with the input relative abundance
  hapsEncodeBySampRelAbund_ = std::vector<std::vector<double> >(sampNames_.size());
  for (const auto &samp: sampNames_) {
    hapsEncodeBySampRelAbund_[sampNamesKey_[samp]] = std::vector<double>(totalHaps_, 0);
  }
  TableReader reReadHapTab(TableIOOpts(InOptions(pars_.tableFnp), "\t", true));
  VecStr row;
  while (reReadHapTab.getNextRow(row)) {
    const auto &samp = row[reReadHapTab.header_.getColPos(pars_.sampleCol)];
    const auto &tar = row[reReadHapTab.header_.getColPos(pars_.targetNameCol)];
    if(!pars_.selectSamples.empty() && !njh::in(samp, pars_.selectSamples)){
      continue;
    }
    if(!pars_.exclude_samples.empty() && njh::in(samp, pars_.exclude_samples)){
      continue;
    }
    if(!pars_.selectTargets.empty() && !njh::in(tar, pars_.selectTargets)){
      continue;
    }
    if(!pars_.exclude_targets.empty() && njh::in(tar, pars_.exclude_targets)){
      continue;
    }
    const auto &hapName = row[reReadHapTab.header_.getColPos(pars_.popIDCol)];
    auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(pars_.relAbundCol)]);
    auto tKey = tarNameKey_[tar];
    auto hKey = hapNamesKey_[tar][hapName];
    hapsEncodeBySampRelAbund_[sampNamesKey_[samp]][tarStart_[tKey] + hKey] = rBund;
  }
  //now recalculate the relative abundance to be 0-1
  for (const auto pos: iter::range(sampNames_.size())) {
    for (const auto tpos: iter::range(tarNamesVec_.size())) {
      if (targetsEncodeBySamp_[pos][tpos] == 1) {
        double sum = 0;
        for (const auto hapPos: iter::range(numberOfHapsPerTarget_[tpos])) {
          sum += hapsEncodeBySampRelAbund_[pos][tarStart_[tpos] + hapPos];
        }
        for (const auto hapPos: iter::range(numberOfHapsPerTarget_[tpos])) {
          hapsEncodeBySampRelAbund_[pos][tarStart_[tpos] + hapPos] =
              hapsEncodeBySampRelAbund_[pos][tarStart_[tpos] + hapPos] / sum;
        }
      }
    }
  }
}


std::unordered_map<std::string, uint32_t> HapsEncodedMatrix::getNumberTargetsPerSample() const {
	std::unordered_map<std::string, uint32_t> ret;
	for (const auto row : iter::range(targetsEncodeBySamp_.size())) {
		ret.emplace(sampNamesVec_[row], vectorSum(targetsEncodeBySamp_[row]));
	}
	return ret;
}
std::unordered_map<std::string, double> HapsEncodedMatrix::getTargetCoveragePerSample() const {
	std::unordered_map<std::string, double> lociCoveragePerSample;
	for (const auto row : iter::range(targetsEncodeBySamp_.size())) {
		lociCoveragePerSample[sampNamesVec_[row]] = vectorSum(targetsEncodeBySamp_[row])/static_cast<double>(numberOfHapsPerTarget_.size());
	}
	return lociCoveragePerSample;
}


table HapsEncodedMatrix::exportEncodedTable() const{
	return exportEncodedTable(ExportEncodedTablePars{});
}

table HapsEncodedMatrix::exportEncodedTable(const ExportEncodedTablePars & exportPars) const{
	const bool hasRelAbund = !hapsEncodeBySampRelAbund_.empty();
	VecStr header{exportPars.sampleColName, exportPars.targetColName, exportPars.hapColName};
	if(hasRelAbund){
		header.emplace_back(exportPars.relAbundColName);
	}
	table ret(header);
	//invert the haplotype keys so the haplotype index within a target can be used to look up the haplotype name
	std::vector<std::vector<std::string>> hapNamesByTarKey(tarNamesVec_.size());
	for(const auto tarPos : iter::range(tarNamesVec_.size())){
		hapNamesByTarKey[tarPos] = std::vector<std::string>(numberOfHapsPerTarget_[tarPos]);
	}
	for(const auto & tarHaps : hapNamesKey_){
		const auto tarPos = tarNameKey_.at(tarHaps.first);
		for(const auto & hap : tarHaps.second){
			hapNamesByTarKey[tarPos][hap.second] = hap.first;
		}
	}
	for(const auto sampPos : iter::range(hapsEncodeBySamp_.size())){
		for(const auto tarPos : iter::range(tarNamesVec_.size())){
			if(0 == targetsEncodeBySamp_[sampPos][tarPos]){
				//sample has no data for this target
				continue;
			}
			for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tarPos])){
				if(0 == hapsEncodeBySamp_[sampPos][tarStart_[tarPos] + hapPos]){
					continue;
				}
				if(hasRelAbund){
					ret.addRow(sampNamesVec_[sampPos], tarNamesVec_[tarPos],
							hapNamesByTarKey[tarPos][hapPos],
							hapsEncodeBySampRelAbund_[sampPos][tarStart_[tarPos] + hapPos]);
				}else{
					ret.addRow(sampNamesVec_[sampPos], tarNamesVec_[tarPos],
							hapNamesByTarKey[tarPos][hapPos]);
				}
			}
		}
	}
	return ret;
}

table HapsEncodedMatrix::getTableNumberTargetsPerSample(double coverage_cut_off) const{
	VecStr header{"sample", "targetCount", "target_coverage"};
	if (coverage_cut_off != std::numeric_limits<double>::min()) {
		header.emplace_back("above_coverage_cut_off");
	}
	table ret(header);
	auto coverage = getTargetCoveragePerSample();
	auto count = getNumberTargetsPerSample();
	for (const auto & cov : coverage) {
		if (coverage_cut_off != std::numeric_limits<double>::min()) {
			ret.addRow(cov.first, count[cov.first], cov.second, cov.second >= coverage_cut_off);
		} else {
			ret.addRow(cov.first, count[cov.first], cov.second);
		}
	}
	ret.sortTable("sample", true);

	return ret;
}

void HapsEncodedMatrix::addMetaWithInputTab(const std::set<std::string> & metaFields){
	TableReader inputTab(TableIOOpts::genTabFileIn(pars_.tableFnp));
	std::vector<std::string> tableCheckerCols{pars_.sampleCol};
	njh::addConToVec(tableCheckerCols, metaFields);
	inputTab.header_.checkForColumnsThrow(tableCheckerCols, __PRETTY_FUNCTION__);
	table metaTab(tableCheckerCols);
	VecStr row;
	while(inputTab.getNextRow(row)){
		metaTab.addRow(inputTab.extractCols(row, tableCheckerCols));
	}
	//have to rename for the MultipleGroupMetaData class, it wants sample, Sample, samples, or Samples
	metaTab.columnNames_[metaTab.getColPos(pars_.sampleCol)] = "sample";
	metaTab.setColNamePositions();

	meta_ = std::make_shared<MultipleGroupMetaData>(metaTab.getUniqueRows(), njh::vecToSet(metaTab.getColumnLevels("sample")));
}


void HapsEncodedMatrix::addMeta(const bfs::path & metaFnp){
	if(encodeKeysSet_){
		meta_ = std::make_shared<MultipleGroupMetaData>(metaFnp, std::set<std::string>(sampNames_.begin(), sampNames_.end()));
	}else{
		meta_ = std::make_shared<MultipleGroupMetaData>(metaFnp);
	}
}

void HapsEncodedMatrix::resetEncoding(){
	sampNamesVec_.clear();
	sampNamesKey_.clear();
	tarNamesVec_.clear();
	tarNameKey_.clear();
	hapNamesKey_.clear();
	numberOfHapsPerTarget_.clear();
	tarStart_.clear();
	hapsEncodeBySamp_.clear();
	targetsEncodeBySamp_.clear();
	totalHaps_ = 0;

}
void HapsEncodedMatrix::setEncodeKeys(){
	if(encodeKeysSet_){
		//reset
		resetEncoding();
	}
	encodeKeysSet_ = true;
	sampNamesVec_ = std::vector<std::string> (sampNames_.begin(), sampNames_.end());
	tarNamesVec_ = std::vector<std::string> (tarNames_.begin(), tarNames_.end());
	njh::sort(sampNamesVec_);
	njh::sort(tarNamesVec_);

	for(const auto pos : iter::range(tarNamesVec_.size())){
		tarNameKey_[tarNamesVec_[pos]] = pos;
	}

	for(const auto pos : iter::range(sampNamesVec_.size())){
		sampNamesKey_[sampNamesVec_[pos]] = pos;
	}

	for(const auto & hapName : hapNamesForTars_){
		std::vector<std::string> hapNamesVec = getVectorOfMapKeys(hapName.second);
		//sort by hap sample abundance
		njh::sort(hapNamesVec, [&hapName](const std::string & hap1, const std::string & hap2){
			return hapName.second.at(hap1) > hapName.second.at(hap2);
		});
		for(const auto pos : iter::range(hapNamesVec.size())){
			hapNamesKey_[hapName.first][hapNamesVec[pos]] = pos;
			++totalHaps_;
		}
	}
	numberOfHapsPerTarget_ = std::vector<uint32_t> (tarNamesVec_.size());
	for(const auto & tarHaps : hapNamesForTars_){
		numberOfHapsPerTarget_[tarNameKey_[tarHaps.first]] = tarHaps.second.size();
	}
	tarStart_ = std::vector<uint32_t> (tarNamesVec_.size());

	{
		uint32_t pos = 0;
		for(const auto tarPos : iter::range(tarNamesVec_.size())){
			tarStart_[tarPos] = pos;
			pos += numberOfHapsPerTarget_[tarPos];
		}
	}
	hapsEncodeBySamp_ = std::vector<std::vector<uint8_t>> (sampNames_.size());
	for(const auto & samp : sampNames_){
		hapsEncodeBySamp_[sampNamesKey_[samp]] = std::vector<uint8_t>(totalHaps_, 0);
	}

	targetsEncodeBySamp_ = std::vector<std::vector<uint8_t>> (sampNames_.size());
	for(const auto & samp : sampNames_){
		targetsEncodeBySamp_[sampNamesKey_[samp]] = std::vector<uint8_t>(tarNamesVec_.size(), 0);
	}
}

HapsEncodedMatrix::IndexResults::IndexResults(const uint64_t numOfSamps){

	byAllHaps = std::vector<std::vector<double>>(numOfSamps, std::vector<double>(numOfSamps));
	byHapsTarShared = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));

	avgJacard = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));
	byTarget = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));

	byHapsTarSharedWeighted = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));
	avgJacardWeighted = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));
	targetsShared = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));



	for(uint32_t pos = 0; pos < numOfSamps; ++pos){
		//set diagonal
		byHapsTarShared[pos][pos] = 1;
		byAllHaps[pos][pos] = 1;

		avgJacard[pos][pos] = 1;
		byTarget[pos][pos] = 1;

		byHapsTarSharedWeighted[pos][pos] = 1;
		avgJacardWeighted[pos][pos] = 1;

		targetsShared[pos][pos] = 1;
	}
}
HapsEncodedMatrix::CCCRMSEResults::CCCRMSEResults(const uint64_t numOfSamps) {
  rmse = std::vector<std::vector<double>>(numOfSamps, std::vector<double>(numOfSamps,1.0));
  ccc =  std::vector<std::vector<double>>(numOfSamps, std::vector<double>(numOfSamps,0.0));
  targets_shared = std::vector<std::vector<uint32_t>>(numOfSamps, std::vector<uint32_t>(numOfSamps,0U));

  //set diagonal
  for(size_t pos = 0; pos < numOfSamps; ++pos){
    rmse[pos][pos] = 0;
    ccc[pos][pos]  = 1.0;
  }
}


void HapsEncodedMatrix::writeAbsoluteHapSharedPerSamplePerTar(const OutOptions & outOptions, bool verbose) const{
	OutputStream out(outOptions);
	out << "target\tsample";
	for(const auto & samp : sampNamesVec_){
		out << "\t" << samp;
	}
	out << std::endl;

	for(const auto tpos : iter::range(tarNamesVec_.size())){
		for(const auto samp1Pos : iter::range(sampNamesVec_.size())){
			out << tarNamesVec_[tpos] << "\t" << sampNamesVec_[samp1Pos];
			for(const auto samp2Pos : iter::range(sampNamesVec_.size())){
				uint8_t tarRes = targetsEncodeBySamp_[samp1Pos][tpos] + targetsEncodeBySamp_[samp2Pos][tpos];
				//should be 2 if both samples have this target
				uint32_t totalSharedForTar = 0;
				if(2 == tarRes){
					for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tpos])){
						//position in the encoded vector should be the target start ranged over the possible haplotypes for that target
						uint8_t res = hapsEncodeBySamp_[samp1Pos][tarStart_[tpos] + hapPos] + hapsEncodeBySamp_[samp2Pos][tarStart_[tpos] + hapPos];
						//if res is 2 then haps are shared
						if(2 == res) {
							++totalSharedForTar;
						}
					}
				}
				out << "\t" << totalSharedForTar;
			}
			out << std::endl;
		}
	}
}


HapsEncodedMatrix::CCCRMSEResults HapsEncodedMatrix::calc_ccc_rmse_measures(uint32_t bin_batch_size, bool verbose) const {
  CCCRMSEResults ret(sampNamesVec_.size());
  for(size_t pos = 0; pos < sampNamesVec_.size(); ++pos){
    ret.targets_shared[pos][pos] = vectorSum(targetsEncodeBySamp_[pos]);
  }
  PairwisePairFactory pFactor(sampNames_.size());
  njh::ProgressBar progpar(pFactor.totalCompares_);
  if(verbose){
    std::cout << "totalComps: " << pFactor.totalCompares_ << std::endl;
  }
  std::function<void()> compSamps = [&pFactor,
        &progpar,
        this,
        &ret,
        &verbose,
        bin_batch_size]() {
    PairwisePairFactory::PairwisePairVec pairVec;

    while(pFactor.setNextPairs(pairVec, bin_batch_size)) {
      if(verbose){
        progpar.outputProgAdd(std::cout, pairVec.pairs_.size(), true);
      }
      for(const auto pairPos : iter::range(pairVec.pairs_.size())) {
        const auto & pair = pairVec.pairs_[pairPos];
        std::vector<double> row_values;
        std::vector<double> col_values;
        // std::vector<double> rmses; could consider calculating a mean RMSE per target
        double sum = 0;
        uint32_t targets_shared = 0;
        uint32_t total_haps_shared = 0;
        for(const auto tpos : iter::range(tarNamesVec_.size())) {
          uint8_t tarRes = targetsEncodeBySamp_[pair.col_][tpos] + targetsEncodeBySamp_[pair.row_][tpos];
          //should be 2 if both samples have this target
          if(2 == tarRes) {
            ++targets_shared;
            // double current_sum = 0;
            double haps_shared_for_target = 0;
            for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tpos])) {
              uint8_t res = hapsEncodeBySamp_[pair.col_][tarStart_[tpos] + hapPos] + hapsEncodeBySamp_[pair.row_][tarStart_[tpos] + hapPos];
              if (res > 0) {
                //have to be > 0 for either sample to have this haplotype
                row_values.emplace_back(hapsEncodeBySampRelAbund_[pair.row_][tarStart_[tpos] + hapPos]);
                col_values.emplace_back(hapsEncodeBySampRelAbund_[pair.col_][tarStart_[tpos] + hapPos]);
                // current_sum += std::pow(hapsEncodeBySampRelAbund_[pair.col_][tarStart_[tpos] + hapPos] - hapsEncodeBySampRelAbund_[pair.row_][tarStart_[tpos] + hapPos],2);
                sum +=         std::pow(hapsEncodeBySampRelAbund_[pair.col_][tarStart_[tpos] + hapPos] - hapsEncodeBySampRelAbund_[pair.row_][tarStart_[tpos] + hapPos],2);
                ++haps_shared_for_target;
              }
            }
            total_haps_shared+= haps_shared_for_target;
            // rmses.emplace_back(std::sqrt(current_sum/haps_shared_for_target));
          }
        }

        //for when there are no targets shared between samples

        // auto ccc_calc = row_values.size() > 0 ? lins_concordance_correlation(row_values, col_values): 0.0;
        auto ccc_calc = row_values.size() > 0 ? ConcordanceCalculator::lins_ccc_point(row_values, col_values): 0.0;
        auto rmse_calc = total_haps_shared > 0 ? std::sqrt(sum/total_haps_shared): 1.0; //making this 1 for now but really should be NA as it's not calculable for non-sharing samples
        ret.rmse[pair.col_][pair.row_] = rmse_calc;
        ret.rmse[pair.row_][pair.col_] = rmse_calc;
        ret.ccc[pair.col_][pair.row_] = ccc_calc;
        ret.ccc[pair.row_][pair.col_] = ccc_calc;
        ret.targets_shared[pair.row_][pair.col_] = targets_shared;
        ret.targets_shared[pair.col_][pair.row_] = targets_shared;
      }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(compSamps, pars_.numThreads);
  return ret;
}

HapsEncodedMatrix::IndexResults HapsEncodedMatrix::genIndexMeasures(uint32_t bin_batch_size, bool verbose) const{
	IndexResults ret(sampNames_.size());

	PairwisePairFactory pFactor(sampNames_.size());


	njh::ProgressBar progpar(pFactor.totalCompares_);
	if(verbose){
		std::cout << "totalComps: " << pFactor.totalCompares_ << std::endl;
	}


	std::function<void()> compSamps = [&pFactor,&progpar,
																		 this,
																		 &ret, &bin_batch_size,
																		 &verbose](){
		PairwisePairFactory::PairwisePairVec pairVec;
		while(pFactor.setNextPairs(pairVec, bin_batch_size)){
		  if(verbose){
		    progpar.outputProgAdd(std::cout, pairVec.pairs_.size(), true);
		  }
			for(const auto pairPos : iter::range(pairVec.pairs_.size())){
				const auto & pair = pairVec.pairs_[pairPos];
				//for shared targets
				{
					uint32_t totalSet = 0;
					uint32_t totalShared = 0;
					double totalSetWeighted = 0;
					double totalSharedWeighted = 0;
					uint32_t totalTarsWithDataForBoth = 0;
					uint32_t totalTarsWithAtLeastOneHapShared = 0;
					std::vector<double> jacardsByTarget;
					std::vector<double> jacardsByTargetWeighted;

					for(const auto tpos : iter::range(tarNamesVec_.size())){
						uint8_t tarRes = targetsEncodeBySamp_[pair.col_][tpos] + targetsEncodeBySamp_[pair.row_][tpos];
						//should be 2 if both samples have this target
						if(2 == tarRes){
							++totalTarsWithDataForBoth;
							uint32_t totalSetForTar = 0;
							uint32_t totalSharedForTar = 0;
							double totalWeightedSetForTar = 0;
							double totalWeightedSharedForTar = 0;
							for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tpos])){
								//position in the encoded vector should be the target start ranged over the possible haplotypes for that target
								uint8_t res = hapsEncodeBySamp_[pair.col_][tarStart_[tpos] + hapPos] + hapsEncodeBySamp_[pair.row_][tarStart_[tpos] + hapPos];
								//if res is 2 then haps are shared
								if(2 == res) {
									++totalShared;
									++totalSharedForTar;
									totalWeightedSharedForTar += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
									totalSharedWeighted += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
								}
								//if results is either 1 or 2 then at least one of them has this hap
								if(res > 0) {
									++totalSet;
									++totalSetForTar;
									totalSetWeighted += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
									totalWeightedSetForTar += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
								}
							}
							if(totalSharedForTar > 0){
								//at least one haplotype is shared between samps
								++totalTarsWithAtLeastOneHapShared;
							}
							jacardsByTarget.emplace_back(totalSharedForTar/static_cast<double>(totalSetForTar));
							jacardsByTargetWeighted.emplace_back(totalWeightedSharedForTar/(totalWeightedSetForTar));
						}
					}
					ret.byHapsTarShared[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
					ret.byHapsTarShared[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
					ret.byHapsTarSharedWeighted[pair.row_][pair.col_] = totalSharedWeighted/(totalSetWeighted);
					ret.byHapsTarSharedWeighted[pair.col_][pair.row_] = totalSharedWeighted/(totalSetWeighted);
					auto meanJacard = vectorMean(jacardsByTarget);;
					auto meanJacardWeighted = vectorMean(jacardsByTargetWeighted);;
					ret.avgJacard[pair.row_][pair.col_] = meanJacard;
					ret.avgJacard[pair.col_][pair.row_] = meanJacard;
					ret.avgJacardWeighted[pair.row_][pair.col_] = meanJacardWeighted;
					ret.avgJacardWeighted[pair.col_][pair.row_] = meanJacardWeighted;
					ret.byTarget[pair.row_][pair.col_] = totalTarsWithAtLeastOneHapShared/static_cast<double>(totalTarsWithDataForBoth);
					ret.byTarget[pair.col_][pair.row_] = totalTarsWithAtLeastOneHapShared/static_cast<double>(totalTarsWithDataForBoth);

					ret.targetsShared[pair.row_][pair.col_] = totalTarsWithDataForBoth;
					ret.targetsShared[pair.col_][pair.row_] = totalTarsWithDataForBoth;

				}
				//for all
				{
					uint32_t totalSet = 0;
					uint32_t totalShared = 0;
					for(const auto pos : iter::range(totalHaps_)){
						uint8_t res = hapsEncodeBySamp_[pair.col_][pos] + hapsEncodeBySamp_[pair.row_][pos];
						if(2 == res){
							++totalShared;
						}
						if(res > 0){
							++totalSet;
						}
					}
//						if(("8034209115" == sampNamesVec_[pair.row_] && "8025874502" == sampNamesVec_[pair.col_]) ||
//							 ("8025874502" == sampNamesVec_[pair.row_] && "8034209115" == sampNamesVec_[pair.col_]) ){
//							std::cout << sampNamesVec_[pair.row_] << ":" << sampNamesVec_[pair.col_] << std::endl;
//							std::cout << "pair.row_: " << pair.row_ << " " << sampNamesVec_[pair.row_] << std::endl;
//							std::cout << "pair.col_: " << pair.col_ << " " << sampNamesVec_[pair.col_]<< std::endl;
//							std::cout << "totalShared: " << totalShared << std::endl;
//							std::cout << "totalSet: " << totalSet << std::endl;
//							std::cout << "totalShared/static_cast<double>(totalSet): " << totalShared/static_cast<double>(totalSet) << std::endl;
//						}
					ret.byAllHaps[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
					ret.byAllHaps[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
				}
			}
		}
	};


	for (const auto row : iter::range(targetsEncodeBySamp_.size())) {
		ret.targetsShared[row][row] = vectorSum(targetsEncodeBySamp_[row]);
	}

	njh::concurrent::runVoidFunctionThreaded(compSamps, pars_.numThreads);

	return ret;
}


}  // namespace njhseq
