/*
 * bamExtractUtils.cpp
 *
 *  Created on: Feb 2, 2017
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
#include "bamExtractUtils.hpp"
#include "njhseq/BamToolsUtils/BamAlnsCache.hpp"
#include "njhseq/BamToolsUtils/BamAlnsCacheWithRegion.hpp"
#include "njhseq/IO/SeqIO.h"
#include "njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp"

namespace njhseq {

uint32_t getAlnLen(const BamTools::BamAlignment bAln){
	if(bAln.IsMapped()){
		return bAln.GetEndPosition() - bAln.Position;
	}
	return 0;
}


uint32_t BamExtractor::ExtractCounts::getTotal(){

	return pairedReads_ + pairedReadsMateUnmapped_ + unpaiedReads_ + orphans_ + orphansUnmapped_ + mateFilteredOff_ + mateFilteredOffUnmapped_ + pairsUnMapped_ + unpairedUnMapped_ + discordant_ + inverse_ + pairedReadsBothFailedSoftClip_ + pairedReadsMateFailedSoftClip_ + unpairedFailedSoftClip_ + bothMatesFilteredOff_ + singlesFilteredOff_ + orphansFiltered_ + pairedReadsMateUnmappedFailedSoftClip_ + orphansFilteredSoftCip_;
}


void BamExtractor::ExtractCounts::log(std::ostream & out, const bfs::path & bamFnp) {
	out << "file\treadStatus\tcount\tpercent" << std::endl;
	uint32_t total = getTotal();
	auto logStats = [&total,&out,&bamFnp](const std::string & status, uint32_t count){
		out << bamFnp.string()
				<< "\t" << status
				<< "\t" << count
				<< "\t" << (0 == total ? 0 : count/static_cast<double>(total) )
				<< std::endl;
	};
	logStats("paired", pairedReads_);
	logStats("pairedMateUnmapped", pairedReadsMateUnmapped_);
	logStats("pairedReadsMateUnmappedFailedSoftClip", pairedReadsMateUnmappedFailedSoftClip_);
	logStats("pairedReadsMateFailedSoftClip", pairedReadsMateFailedSoftClip_);
	logStats("pairedReadsBothFailedSoftClip", pairedReadsBothFailedSoftClip_);
	logStats("singles", unpaiedReads_);
	logStats("singlesFailedSoftClip", unpairedFailedSoftClip_);
	logStats("singlesFilteredOff", singlesFilteredOff_);
	logStats("orphans", orphans_);
	logStats("orphansFiltered", orphansFiltered_);
	logStats("orphansFilteredSoftCip", orphansFilteredSoftCip_);
	logStats("orphansUnmapped", orphansUnmapped_);
	logStats("bothMatesFilteredOff", bothMatesFilteredOff_);
	logStats("mateFilteredOff", mateFilteredOff_);
	logStats("mateFilteredOffUnmapped", mateFilteredOffUnmapped_);
	logStats("unmappedPaired", pairsUnMapped_);
	logStats("unmappedSingles", unpairedUnMapped_);
	logStats("discordant", discordant_);
	logStats("inverse", inverse_);
	out << bamFnp.string()
			<< "\t" << "improperPairFiltered"
			<< "\t" << improperPairFiltered_
			<< "\t" << ""
			<< std::endl;
	out << bamFnp.string()
			<< "\t" << "markedDuplicateFiltered"
			<< "\t" << markedDuplicateFiltered_
			<< "\t" << ""
			<< std::endl;
}

void BamExtractor::ExtractedFilesOpts::removeAllInFiles(){
	auto removeIfExists = [](const bfs::path & fnp){
		if(bfs::exists(fnp)){
			bfs::remove(fnp);
		}
	};
	removeIfExists(inPairs_.firstName_);
	removeIfExists(inPairs_.secondName_);

	removeIfExists(inUnpaired_.firstName_);
	removeIfExists(inFilteredSingles_.firstName_);
	removeIfExists(inSoftClipFilteredSingles_.firstName_);

	removeIfExists(inPairsUnMapped_.firstName_);
	removeIfExists(inPairsUnMapped_.secondName_);

	removeIfExists(inUnpairedUnMapped_.firstName_);

	removeIfExists(inInverse_.firstName_);
	removeIfExists(inInverse_.secondName_);

	removeIfExists(inFilteredPairs_.firstName_);
	removeIfExists(inFilteredPairs_.secondName_);

	removeIfExists(inSoftClipFilteredPairs_.firstName_);
	removeIfExists(inSoftClipFilteredPairs_.secondName_);

	removeIfExists(inDiscordant_.firstName_);
	removeIfExists(inDiscordant_.secondName_);




}

BamExtractor::BamExtractor(bool verbose ):verbose_(verbose) {

}

BamExtractor::BamExtractSeqsResults BamExtractor::extractReadsFromBamRegion(
		const bfs::path & bamFnp, const GenomicRegion & region,
		double percInRegion) {
	BamExtractSeqsResults ret;
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	if (!bReader.SetRegion(refNameToId.at(region.chrom_), region.start_,
			refNameToId.at(region.chrom_), region.end_)) {
		std::stringstream ss;
		ss << bReader.GetErrorString() << std::endl;
		ss << "Failed to set region" << std::endl;
		ss << "Region: " << refNameToId.at(region.chrom_) << std::endl;
		ss << "Start: " << region.start_ << std::endl;
		ss << "Stop: " << region.end_ << std::endl;
		throw std::runtime_error { ss.str() };
	}

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//get only alignments that fall mostly in this region, setting percInRegion
		//would save any read that had any fall bases in this region
		if (region.getPercInRegion(bAln, refData) < percInRegion) {
			continue;
		}
		if (bAln.IsPaired()) {
			if (bAln.MateRefID != bAln.RefID) {
				// do non-concordant chromosome mapping operation
				continue;
			}
			//uncommenting this will make it so unmapping mate will now be recovered
//			if (!bAln.IsMapped() || !bAln.IsMateMapped()) {
//				// do non-mapping operation
//				continue;
//			}

			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this position
					//it should be in the cache therefore program was unable to find mate

					//do orphaned operation
					ret.singlets_.emplace_back(bamAlnToSeqInfo(bAln));
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					//need to check to see if there is an overlap
					if (search->GetEndPosition() > bAln.Position) {
						//overlap

						//do overlapped pairs operation

					} else {
						//no overlap, single read count both

						//do no overlap pairs operation
					}
					if (bAln.IsFirstMate()) {
						ret.pairs_.emplace_back(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						ret.pairs_.emplace_back(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			//unpaired read
			if (!bAln.IsMapped()) {
				// do non-mapping operation
			} else {
				//do unpaired read operation
				ret.singlets_.emplace_back(bamAlnToSeqInfo(bAln));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			ret.singlets_.emplace_back(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
	return ret;
}


void BamExtractor::writeExtractReadsFromBamRegion(
		const bfs::path & bamFnp,
		const GenomicRegion & region,
		double percInRegion,
		const OutOptions & outOpts) {

	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQPAIRED, outOpts);
	SeqOutput pairWriter(outOptsPaired);
	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQ, outOpts);
	SeqOutput writer(outOptsSingle);

	if (!bReader.SetRegion(refNameToId.at(region.chrom_), region.start_,
			refNameToId.at(region.chrom_), region.end_)) {
		std::stringstream ss;
		ss << bReader.GetErrorString() << std::endl;
		ss << "Failed to set region" << std::endl;
		ss << "Region: " << refNameToId.at(region.chrom_) << std::endl;
		ss << "Start: " << region.start_ << std::endl;
		ss << "Stop: " << region.end_ << std::endl;
		throw std::runtime_error { ss.str() };
	}

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//get only alignments that fall mostly in this region, setting percInRegion
		//would save any read that had any fall bases in this region
		if (region.getPercInRegion(bAln, refData) < percInRegion) {
			continue;
		}
		if (bAln.IsPaired()) {
			if (bAln.MateRefID != bAln.RefID) {
				// do non-concordant chromosome mapping operation
				continue;
			}
			//uncommenting this will make it so unmapping mate will now be recovered
//			if (!bAln.IsMapped() || !bAln.IsMateMapped()) {
//				// do non-mapping operation
//				continue;
//			}

			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this position
					//it should be in the cache therefore program was unable to find mate
					//do orphaned operation
					//this could be due mate not falling in this region;
					writer.openWrite(bamAlnToSeqInfo(bAln));
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					PairedRead outRead;
					if (bAln.IsFirstMate()) {
						outRead = PairedRead(bamAlnToSeqInfo(bAln),bamAlnToSeqInfo(*search), false);

					} else {
						outRead = PairedRead(bamAlnToSeqInfo(*search),bamAlnToSeqInfo(bAln), false);
					}
					outRead.seqBase_.name_.append("/1");
					outRead.mateSeqBase_.name_.append("/2");
					pairWriter.openWrite(outRead);
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			//unpaired read
			if (!bAln.IsMapped()) {
				// do non-mapping operation
			} else {
				//do unpaired read operation
				writer.openWrite(bamAlnToSeqInfo(bAln));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			writer.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
}

BamExtractor::ExtractedFilesWithStichingOpts BamExtractor::writeExtractReadsFromBamRegionStitch(
		const bfs::path & bamFnp, const GenomicRegion & region, double percInRegion,
		const OutOptions & outOpts, const std::string & extraFlashCms) {
	njh::sys::requireExternalProgramThrow("flash");

	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired = SeqIOOptions::genPairedOut(outOpts.outFilename_);
	outOptsPaired.out_.overWriteFile_ = outOpts.overWriteFile_;
	outOptsPaired.out_.append_ = outOpts.append_;
	SeqOutput pairWriter(outOptsPaired);
	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQ, outOpts);
	SeqOutput writer(outOptsSingle);


	ExtractedFilesWithStichingOpts ret;
	ret.inPairs_ = SeqIOOptions::genPairedIn(outOptsPaired.getPriamryOutName(),
			outOptsPaired.getSecondaryOutName());
	ret.inUnpaired_ = SeqIOOptions::genFastqIn(
			outOptsSingle.getPriamryOutName());
	if (!bReader.SetRegion(refNameToId.at(region.chrom_), region.start_,
			refNameToId.at(region.chrom_), region.end_)) {
		std::stringstream ss;
		ss << bReader.GetErrorString() << std::endl;
		ss << "Failed to set region" << std::endl;
		ss << "Region: " << refNameToId.at(region.chrom_) << std::endl;
		ss << "Start: " << region.start_ << std::endl;
		ss << "Stop: " << region.end_ << std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::vector<uint32_t> readLens;
	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//get only alignments that fall mostly in this region, setting percInRegion
		//would save any read that had any fall bases in this region
		if (region.getPercInRegion(bAln, refData) < percInRegion) {
			continue;
		}
		if (bAln.IsPaired()) {
			if (bAln.MateRefID != bAln.RefID) {
				// do non-concordant chromosome mapping operation
				continue;
			}
//			if (!bAln.IsMapped() || !bAln.IsMateMapped()) {
//				// do non-mapping operation
//				continue;
//			}

			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this position
					//it should be in the cache therefore program was unable to find mate
					//do orphaned operation
					//this could be due mate not falling in this region;
					writer.openWrite(bamAlnToSeqInfo(bAln));
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					//need to check to see if there is an overlap
					if (search->GetEndPosition() > bAln.Position) {
						//overlap

						//do overlapped pairs operation

					} else {
						//no overlap, single read count both

						//do no overlap pairs operation
					}
					readLens.emplace_back(bAln.QueryBases.length());
					readLens.emplace_back(search->QueryBases.length());
					if (bAln.IsFirstMate()) {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));

					} else {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			//unpaired read
			if (!bAln.IsMapped()) {
				// do non-mapping operation
			} else {
				//do unpaired read operation
				writer.openWrite(bamAlnToSeqInfo(bAln));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			writer.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
	writer.closeOut();
	pairWriter.closeOut();
	if (outOptsPaired.outExists()) {
		std::stringstream flashCmd;
		auto stitchFilesBase = njh::files::prependFileBasename(outOpts.outName(),
				"stiched_");
		flashCmd << "flash " << outOptsPaired.getPriamryOutName() << " "
				<< outOptsPaired.getSecondaryOutName() << " -o " << stitchFilesBase
				<< " " << extraFlashCms;
		if (std::string::npos != extraFlashCms.find("--max-overlap")) {
			flashCmd << " --max-overlap "
					<< estd::to_string(std::round(vectorMedianRef(readLens))) << " ";
		}
		if (std::string::npos != extraFlashCms.find("-x")) {
			flashCmd << " -x 0.02 ";
		}
		flashCmd << " > "
				<< njh::files::make_path(stitchFilesBase.parent_path(), "flashLog")
				<< " 2>&1\n";
		auto flashOut = njh::sys::run(VecStr { flashCmd.str() });
		if (!flashOut.success_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error in running flash to stitch reads"
					<< "\n";
			ss << "Mess: " << "\n";
			ss << flashOut.stdErr_ << "\n";
			throw std::runtime_error { ss.str() };
		}
		OutOptions flashRunLogOpts(
				njh::files::make_path(stitchFilesBase.parent_path(),
						"flashRunLog.json"));
		std::ofstream outFlashRunLog;
		flashRunLogOpts.openFile(outFlashRunLog);
		outFlashRunLog << flashOut.toJson() << std::endl;
		ret.stitched_ = SeqIOOptions::genFastqIn(
				njh::appendAsNeededRet(stitchFilesBase.string(),
						".extendedFrags.fastq"));

		ret.notCombinedPairs_ = SeqIOOptions::genPairedIn(
				njh::appendAsNeededRet(stitchFilesBase.string(),
						".notCombined_1.fastq"),
				njh::appendAsNeededRet(stitchFilesBase.string(),
						".notCombined_2.fastq"));
	}

	return ret;
}

//BamExtractor::BamExtractSeqsResults BamExtractor::extractReadsFromBam(
//		const bfs::path & bamFnp) {
//	BamExtractSeqsResults ret;
//	BamTools::BamReader bReader;
//	BamTools::BamAlignment bAln;
//	bReader.Open(bamFnp.string());
//	checkBamOpenThrow(bReader);
//	loadBamIndexThrow(bReader);
//	auto refs = bReader.GetReferenceData();
//	BamAlnsCache alnCache;
//	auto refData = bReader.GetReferenceData();
//	std::unordered_map<std::string, uint32_t> refNameToId;
//	for (auto pos : iter::range(refData.size())) {
//		refNameToId[refData[pos].RefName] = pos;
//	}
//
//	while (bReader.GetNextAlignment(bAln)) {
//		//skip secondary alignments
//		if (!bAln.IsPrimaryAlignment()) {
//			continue;
//		}
//
//		if (bAln.IsPaired()) {
//			if (bAln.MateRefID != bAln.RefID) {
//				// do non-concordant chromosome mapping operation
//				continue;
//			}
//			//uncommenting this will make it so unmapping mate will now be recovered
////			if (!bAln.IsMapped() || !bAln.IsMateMapped()) {
////				// do non-mapping operation
////				continue;
////			}
//
//			if (bAln.MatePosition == bAln.Position) {
//				if (!alnCache.has(bAln.Name)) {
//					//if mapped to the same place and the mate is yet to be encountered
//					//enter into cache for until mate is encountered
//					alnCache.add(bAln);
//					continue;
//				}
//			}
//			if (bAln.MatePosition <= bAln.Position) {
//				if (!alnCache.has(bAln.Name)) {
//					//since input should be sorted if matePosition is less than this position
//					//it should be in the cache therefore program was unable to find mate
//
//					//do orphaned operation
//					ret.singlets_.emplace_back(bamAlnToSeqInfo(bAln));
//					continue;
//				} else {
//					auto search = alnCache.get(bAln.Name);
//					//need to check to see if there is an overlap
//					if (search->GetEndPosition() > bAln.Position) {
//						//overlap
//
//						//do overlapped pairs operation
//
//					} else {
//						//no overlap, single read count both
//
//						//do no overlap pairs operation
//					}
//					if (bAln.IsFirstMate()) {
//						ret.pairs_.emplace_back(
//								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
//										false));
//					} else {
//						ret.pairs_.emplace_back(
//								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
//										false));
//					}
//					// now that operations have been computed, remove first mate found from cache
//					alnCache.remove(search->Name);
//				}
//			} else {
//				//enter into cache for until mate is encountered
//				alnCache.add(bAln);
//			}
//		} else {
//			//unpaired read
//			if (!bAln.IsMapped()) {
//				// do non-mapping operation
//			} else {
//				//do unpaired read operation
//				ret.singlets_.emplace_back(bamAlnToSeqInfo(bAln));
//			}
//		}
//	}
//	//save the orphans;
//	if (len(alnCache) > 0) {
//		auto names = alnCache.getNames();
//		for (const auto & name : names) {
//			auto search = alnCache.get(name);
//			ret.singlets_.emplace_back(bamAlnToSeqInfo(*search));
//			alnCache.remove(name);
//		}
//	}
//	return ret;
//}

void BamExtractor::writeExtractReadsFromBam(const bfs::path & bamFnp,
		const OutOptions & outOpts) {

	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired(outOpts.outFilename_, SeqIOOptions::outFormats::FASTQPAIREDGZ, outOpts);
	SeqOutput pairWriter(outOptsPaired);

	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_, SeqIOOptions::outFormats::FASTQGZ, outOpts);
	SeqOutput writer(outOptsSingle);

	//R1 orphans
	SeqIOOptions outOptsSingle_R1Orphans(outOpts.outFilename_.string() + "_R1Orphans", SeqIOOptions::outFormats::FASTQGZ, outOpts);
	SeqOutput writer_R1Orphans(outOptsSingle_R1Orphans);

	//R2 orphans
	SeqIOOptions outOptsSingle_R2Orphans(outOpts.outFilename_.string() + "_R2Orphans", SeqIOOptions::outFormats::FASTQGZ, outOpts);
	SeqOutput writer_R2Orphans(outOptsSingle_R2Orphans);

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		/**@todo consider skipping duplicates */
		if (bAln.IsPaired()) {
			if ((bAln.MateRefID != bAln.RefID)
					|| (!bAln.IsMapped() || !bAln.IsMateMapped())) {
				// do non-concordant chromosome mapping operation or non-mapping mates
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					if (bAln.IsFirstMate()) {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
					continue;
				}
			}
			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this position
					//it should be in the cache therefore program was unable to find mate

					//do orphaned operation
					writer.openWrite(bamAlnToSeqInfo(bAln));
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					//need to check to see if there is an overlap
					if (search->GetEndPosition() > bAln.Position) {
						//overlap

						//do overlapped pairs operation

					} else {
						//no overlap, single read count both

						//do no overlap pairs operation
					}
					if (bAln.IsFirstMate()) {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			//unpaired read
			if (!bAln.IsMapped()) {
				// do non-mapping operation
				writer.openWrite(bamAlnToSeqInfo(bAln));
			} else {
				//do unpaired read operation
				writer.openWrite(bamAlnToSeqInfo(bAln));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			if(search->IsFirstMate()){
				writer_R1Orphans.openWrite(bamAlnToSeqInfo(*search));
			}else{
				writer_R2Orphans.openWrite(bamAlnToSeqInfo(*search));
			}
			alnCache.remove(name);
		}
	}
}


BamExtractor::ExtractCounts BamExtractor::writeExtractReadsFromBamOnlyMapped(const bfs::path & bamFnp,
		const OutOptions & outOpts) {
	BamExtractor::ExtractCounts ret;
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQPAIRED, outOpts);
	SeqOutput pairWriter(outOptsPaired);

	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQ, outOpts);
	SeqOutput writer(outOptsSingle);

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//won't handled sequences with pairs that have a unmmaped mate, will make them orphan reads
		//should improve upon
		if(!bAln.IsMapped()){
			continue;
		}
		/**@todo consider skipping duplicates */
		if (bAln.IsPaired()) {
			if ((bAln.MateRefID != bAln.RefID)
					|| (!bAln.IsMapped() || !bAln.IsMateMapped())) {
				// do non-concordant chromosome mapping operation or non-mapping mates
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					if (bAln.IsFirstMate()) {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					++ret.pairedReads_;
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
					continue;
				}
			}
			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this position
					//it should be in the cache therefore program was unable to find mate
					//do orphaned operation
					writer.openWrite(bamAlnToSeqInfo(bAln));
					++ret.orphans_;
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					//need to check to see if there is an overlap
					if (search->GetEndPosition() > bAln.Position) {
						//overlap

						//do overlapped pairs operation

					} else {
						//no overlap, single read count both

						//do no overlap pairs operation
					}
					if (bAln.IsFirstMate()) {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					++ret.pairedReads_;
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			//unpaired read
			++ret.unpaiedReads_;
			if (!bAln.IsMapped()) {
				// do non-mapping operation
				writer.openWrite(bamAlnToSeqInfo(bAln));
			} else {
				//do unpaired read operation
				writer.openWrite(bamAlnToSeqInfo(bAln));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			++ret.orphans_;
			auto search = alnCache.get(name);
			writer.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
	return ret;
}


BamExtractor::BamExtractSeqsResultsAlns BamExtractor::extractReadsFromBamRegionAlns(
		const bfs::path & bamFnp, const GenomicRegion & region,
		double percInRegion) {
	BamExtractSeqsResultsAlns ret;
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	loadBamIndexThrow(bReader);

	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();

	setBamFileRegionThrow(bReader, region);

	while (bReader.GetNextAlignment(bAln)) {
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//get only alignments that fall mostly in this region, setting percInRegion
		//would save any read that had any fall bases in this region
//		if (bAln.IsMapped() && region.getPercInRegion(bAln, refData) < percInRegion) {
//			continue;
//		}
		if (!bAln.IsPaired()) {
			if (!bAln.IsMapped()) {
				++ret.unpairedUnMapped_;
				ret.singletsUnmapped_.emplace_back(bAln);
			} else {
				if(region.getPercInRegion(bAln, refData) >= percInRegion){
					++ret.unpaiedReads_;
					ret.singlets_.emplace_back(bAln);
				}
			}
		} else {
			if (!alnCache.has(bAln.Name)) {
				//pair hasn't been added to cache yet so add to cache
				//this only works if mate and first mate have the same name
				alnCache.add(bAln);
				continue;
			} else {
				auto search = alnCache.get(bAln.Name);
				if (!bAln.IsMapped() && !bAln.IsMateMapped()) {
					++ret.pairsUnMapped_;
					if (bAln.IsFirstMate()) {
						ret.pairsUnmapped_.emplace_back(bAln, *search);
					} else {
						ret.pairsUnmapped_.emplace_back(*search, bAln);
					}
				} else {
					bool bAlnIn = false;
					if(bAln.IsMapped()){
						bAlnIn = region.getPercInRegion(bAln, refData) >= percInRegion;
					}
					bool searchIn = false;
					if(search->IsMapped()){
						searchIn = region.getPercInRegion(*search, refData) >= percInRegion;
					}
					if(bAln.IsMapped() && bAln.IsMateMapped()){
						if(searchIn && bAlnIn){
							if (bAln.RefID == search->RefID
									&& bAln.IsReverseStrand() != search->IsReverseStrand()
									&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
								++ret.pairedReads_;
								if (bAln.IsFirstMate()) {
									ret.pairs_.emplace_back(bAln, *search);
								} else {
									ret.pairs_.emplace_back(*search, bAln);
								}
							} else {
								if (bAln.IsReverseStrand() != search->IsReverseStrand()
										&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
									++ret.inverse_;
									if (bAln.IsFirstMate()) {
										ret.inversePairs_.emplace_back(bAln, *search);
									} else {
										ret.inversePairs_.emplace_back(*search, bAln);
									}
								} else {
									++ret.discordant_;
									if (bAln.IsFirstMate()) {
										ret.discordantPairs_.emplace_back(bAln, *search);
									} else {
										ret.discordantPairs_.emplace_back(*search, bAln);
									}
								}
							}
						} else if (searchIn && !bAlnIn) {
							++ret.orphans_;
							ret.singlets_.emplace_back(*search);
							ret.thrownAwayMates_.emplace_back(bAln);
						} else if (!searchIn && bAlnIn) {
							++ret.orphans_;
							ret.singlets_.emplace_back(bAln);
							ret.thrownAwayMates_.emplace_back(*search);
						}
					}else{
						if(bAln.IsMapped() && bAlnIn){
							/**@todo should incorporate insert size if possible */
							//first check to see if un mapped mate is likely falling within the region of interest
							size_t posibleMatePosition = bAln.Position;
							size_t possibleMatePostionEnd = 0;
							bool reverseStrand = bAln.IsReverseStrand();
							if(bAln.IsReverseStrand()){
								possibleMatePostionEnd = bAln.Position > static_cast<int64_t>(bAln.QueryBases.size()) ? bAln.Position - bAln.QueryBases.size() : 0;
								posibleMatePosition = possibleMatePostionEnd > search->QueryBases.size() ? possibleMatePostionEnd - search->QueryBases.size() : 0;
							}else{
								posibleMatePosition += bAln.QueryBases.size();
								possibleMatePostionEnd = posibleMatePosition + search->QueryBases.size();
							}
							GenomicRegion possibleMateReg(search->Name, refData[search->RefID].RefName,posibleMatePosition, possibleMatePostionEnd, reverseStrand);
							double mateBases = search->QueryBases.size();

							//if(posibleMatePosition >= searchRegion->start_ && posibleMatePosition < searchRegion->end_){
							if(mateBases > 0 && possibleMateReg.getOverlapLen(region)/mateBases >= percInRegion){
								++ret.pairedReadsMateUnmapped_;
								if (bAln.IsFirstMate()) {
									ret.pairsMateUnmapped_.emplace_back(bAln, *search);
								} else {
									ret.pairsMateUnmapped_.emplace_back(*search, bAln);
								}
							}else{
								++ret.pairedReadsMateUnmapped_;
								ret.singlets_.emplace_back(bAln);
								ret.thrownAwayMates_.emplace_back(*search);
							}
						}else if(search->IsMapped() && searchIn){
							/**@todo should incorporate insert size if possible */
							//first check to see if un mapped mate is likely falling within the region of interest
							size_t posibleMatePosition = search->Position;
							size_t possibleMatePostionEnd = 0;
							bool reverseStrand = search->IsReverseStrand();
							if(search->IsReverseStrand()){
								possibleMatePostionEnd = search->Position > static_cast<int64_t>(search->QueryBases.size()) ? search->Position - search->QueryBases.size() : 0;
								posibleMatePosition = possibleMatePostionEnd > bAln.QueryBases.size() ? possibleMatePostionEnd - bAln.QueryBases.size() : 0;
							}else{
								posibleMatePosition += search->QueryBases.size();
								possibleMatePostionEnd = posibleMatePosition + bAln.QueryBases.size();
							}
							GenomicRegion possibleMateReg(bAln.Name, refData[search->RefID].RefName,posibleMatePosition, possibleMatePostionEnd, reverseStrand);
							double mateBases = bAln.QueryBases.size();

							//if(posibleMatePosition >= searchRegion->start_ && posibleMatePosition < searchRegion->end_){
							if(mateBases > 0 && possibleMateReg.getOverlapLen(region)/mateBases >= percInRegion){
								++ret.pairedReadsMateUnmapped_;
								if (bAln.IsFirstMate()) {
									ret.pairsMateUnmapped_.emplace_back(bAln, *search);
								} else {
									ret.pairsMateUnmapped_.emplace_back(*search, bAln);
								}
							}else{
								++ret.pairedReadsMateUnmapped_;
								ret.singlets_.emplace_back(*search);
								ret.thrownAwayMates_.emplace_back(bAln);
							}
						}
					}
				}
				// now that operations have been computed, remove ther other mate found from cache
				alnCache.remove(search->Name);
				continue;
			}
		}
	}

	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			if(search->IsPaired()){
				if (!search->IsMapped()) {
					++ret.orphansUnmapped_;
					ret.singletsUnmapped_.emplace_back(*search);
				} else {
					++ret.orphans_;
					ret.singlets_.emplace_back(*search);
				}
			}else{
				if (!search->IsMapped()) {
					++ret.unpairedUnMapped_;
					ret.singletsUnmapped_.emplace_back(*search);
				} else {
					++ret.unpaiedReads_;
					ret.singlets_.emplace_back(*search);
				}
			}

			alnCache.remove(name);
		}
	}
	if(verbose_){
		ret.log(std::cout, bamFnp);
	}
//	while (bReader.GetNextAlignment(bAln)) {
//		//skip secondary alignments
//		if (!bAln.IsPrimaryAlignment()) {
//			continue;
//		}
//		//get only alignments that fall mostly in this region, setting percInRegion
//		//would save any read that had any fall bases in this region
//		if (region.getPercInRegion(bAln, refData) < percInRegion) {
//			continue;
//		}
//		if (bAln.IsPaired()) {
//			if (bAln.MateRefID != bAln.RefID) {
//				// do non-concordant chromosome mapping operation
//				continue;
//			}
//			//uncommenting this will make it so unmapping mate will now be recovered
////			if (!bAln.IsMapped() || !bAln.IsMateMapped()) {
////				// do non-mapping operation
////				continue;
////			}
//
//			if (bAln.MatePosition == bAln.Position) {
//				if (!alnCache.has(bAln.Name)) {
//					//if mapped to the same place and the mate is yet to be encountered
//					//enter into cache for until mate is encountered
//					alnCache.add(bAln);
//					continue;
//				}
//			}
//			if (bAln.MatePosition <= bAln.Position) {
//				if (!alnCache.has(bAln.Name)) {
//					//since input should be sorted if matePosition is less than this position
//					//it should be in the cache therefore program was unable to find mate
//
//					//do orphaned operation
//					ret.singlets_.emplace_back(bAln);
//					continue;
//				} else {
//					auto search = alnCache.get(bAln.Name);
//					//need to check to see if there is an overlap
//					if (bAln.IsFirstMate()) {
//						ret.pairs_.emplace_back(bAln, *search);
//					} else {
//						ret.pairs_.emplace_back(*search, bAln);
//					}
//					// now that operations have been computed, remove first mate found from cache
//					alnCache.remove(search->Name);
//				}
//			} else {
//				//enter into cache for until mate is encountered
//				alnCache.add(bAln);
//			}
//		} else {
//			//unpaired read
//			if (!bAln.IsMapped()) {
//				// do non-mapping operation
//			} else {
//				//do unpaired read operation
//				ret.singlets_.emplace_back(bAln);
//			}
//		}
//	}
//	//save the orphans;
//	if (len(alnCache) > 0) {
//		auto names = alnCache.getNames();
//		for (const auto & name : names) {
//			auto search = alnCache.get(name);
//			ret.singlets_.emplace_back(*search);
//			alnCache.remove(name);
//		}
//	}
	return ret;
}

//BamExtractor::BamExtractSeqsResultsAlns BamExtractor::extractReadsFromBamAlns(
//		const bfs::path & bamFnp) {
//	BamExtractSeqsResultsAlns ret;
//	BamTools::BamReader bReader;
//	BamTools::BamAlignment bAln;
//	bReader.Open(bamFnp.string());
//	checkBamOpenThrow(bReader);
//	loadBamIndexThrow(bReader);
//	auto refs = bReader.GetReferenceData();
//	BamAlnsCache alnCache;
//	auto refData = bReader.GetReferenceData();
//	std::unordered_map<std::string, uint32_t> refNameToId;
//	for (auto pos : iter::range(refData.size())) {
//		refNameToId[refData[pos].RefName] = pos;
//	}
//	while (bReader.GetNextAlignment(bAln)) {
//		//skip secondary alignments
//		if (!bAln.IsPrimaryAlignment()) {
//			continue;
//		}
//		if (bAln.IsPaired()) {
//			if (bAln.MateRefID != bAln.RefID) {
//				// do non-concordant chromosome mapping operation
//				continue;
//			}
//			if (!bAln.IsMapped() || !bAln.IsMateMapped()) {
//				// do non-mapping operation
//				continue;
//			}
//
//			if (bAln.MatePosition == bAln.Position) {
//				if (!alnCache.has(bAln.Name)) {
//					//if mapped to the same place and the mate is yet to be encountered
//					//enter into cache for until mate is encountered
//					alnCache.add(bAln);
//					continue;
//				}
//			}
//			if (bAln.MatePosition <= bAln.Position) {
//				if (!alnCache.has(bAln.Name)) {
//					//since input should be sorted if matePosition is less than this position
//					//it should be in the cache therefore program was unable to find mate
//
//					//do orphaned operation
//					ret.singlets_.emplace_back(bAln);
//					continue;
//				} else {
//					auto search = alnCache.get(bAln.Name);
//					//need to check to see if there is an overlap
//					if (bAln.IsFirstMate()) {
//						ret.pairs_.emplace_back(bAln, *search);
//					} else {
//						ret.pairs_.emplace_back(*search, bAln);
//					}
//					// now that operations have been computed, remove first mate found from cache
//					alnCache.remove(search->Name);
//				}
//			} else {
//				//enter into cache for until mate is encountered
//				alnCache.add(bAln);
//			}
//		} else {
//			//unpaired read
//			if (!bAln.IsMapped()) {
//				// do non-mapping operation
//			} else {
//				//do unpaired read operation
//				ret.singlets_.emplace_back(bAln);
//			}
//		}
//	}
//	//save the orphans;
//	if (len(alnCache) > 0) {
//		auto names = alnCache.getNames();
//		for (const auto & name : names) {
//			auto search = alnCache.get(name);
//			ret.singlets_.emplace_back(*search);
//			alnCache.remove(name);
//		}
//	}
//	return ret;
//}


BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsFromBamWriteAsSingles(
		const SeqIOOptions & opts, bool referenceOrientation) {

	auto outUnpaired = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(opts.out_);
	auto outUnpairedUnMapped = SeqIOOptions::genFastqOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outUnpairedUnMapped.out_.transferOverwriteOpts(opts.out_);

	SeqOutput mappedSinglesWriter(outUnpaired);
	SeqOutput unmappedSinglesWriter(outUnpairedUnMapped);

	BamExtractor::ExtractedFilesOpts ret;

	ret.inUnpaired_ = SeqIOOptions::genFastqIn(outUnpaired.getPriamryOutName());
	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped.getPriamryOutName());

	BamTools::BamReader bReader;
	bReader.Open(opts.firstName_.string());
	checkBamOpenThrow(bReader, opts.firstName_.string());

	BamTools::BamAlignment bAln;


	while (bReader.GetNextAlignment(bAln)) {
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if(bAln.IsMapped()){
			++ret.unpaiedReads_;
			if(referenceOrientation){
				mappedSinglesWriter.openWrite(seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
						SangerQualOffset));
			}else{
				mappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
			}
		}else{
			++ret.unpairedUnMapped_;
			if(referenceOrientation){
				unmappedSinglesWriter.openWrite(seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
						SangerQualOffset));
			}else{
				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
			}
		}
	}
	if(verbose_){
		ret.log(std::cout, opts.firstName_);
	}
	return ret;
}

// old with center soft clipping filtering
//BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsFromBamToSameOrientationContigs(
//		const SeqIOOptions & opts,
//		bool throwAwayUnmmpaedMates) {
//
//	auto outPairsUnMapped = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
//	outPairsUnMapped.out_.transferOverwriteOpts(opts.out_);
//	auto outUnpairedUnMapped = SeqIOOptions::genFastqOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
//	outUnpairedUnMapped.out_.transferOverwriteOpts(opts.out_);
//
//	auto outPairs = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
//	outPairs.out_.transferOverwriteOpts(opts.out_);
//	auto outUnpaired = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
//	outUnpaired.out_.transferOverwriteOpts(opts.out_);
//
//	auto outPairsMateUnmapped = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "mateUnmapped_"));
//	outPairsMateUnmapped.out_.transferOverwriteOpts(opts.out_);
//
//	auto thrownAwayMate = SeqIOOptions::genFastqOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "thrownAwayMate_"));
//	thrownAwayMate.out_.transferOverwriteOpts(opts.out_);
//
//
//	auto outPairsInverse = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"));
//	outPairsInverse.out_.transferOverwriteOpts(opts.out_);
//
//	auto outPairsDiscordant = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"));
//	outPairsDiscordant.out_.transferOverwriteOpts(opts.out_);
//
//	SeqOutput unmappedPairWriter(outPairsUnMapped);
//	SeqOutput unmappedSinglesWriter(outUnpairedUnMapped);
//
//	SeqOutput inversePairWriter(outPairsInverse);
//
//	SeqOutput discordantPairWriter(outPairsDiscordant);
//
//	SeqOutput mappedPairWriter(outPairs);
//	SeqOutput mappedSinglesWriter(outUnpaired);
//
//	SeqOutput mateUnmappedPairWriter(outPairsMateUnmapped);
//	SeqOutput thrownAwayMateWriter(thrownAwayMate);
//
//
//
//	BamExtractor::ExtractedFilesOpts ret;
//	ret.inPairs_ = SeqIOOptions::genPairedIn(outPairs.getPriamryOutName(), outPairs.getSecondaryOutName());
//	ret.inPairsUnMapped_ = SeqIOOptions::genPairedIn(outPairsUnMapped.getPriamryOutName(), outPairsUnMapped.getSecondaryOutName());
//
//	ret.inUnpaired_ = SeqIOOptions::genFastqIn(outUnpaired.getPriamryOutName());
//	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped.getPriamryOutName());
//
//	ret.inDiscordant_ = SeqIOOptions::genPairedIn(outPairsDiscordant.getPriamryOutName(), outPairsDiscordant.getSecondaryOutName());
//	ret.inInverse_ = SeqIOOptions::genPairedIn(outPairsInverse.getPriamryOutName(), outPairsInverse.getSecondaryOutName());
//
//	ret.inPairsMateUnmapped_ = SeqIOOptions::genPairedIn(outPairsMateUnmapped.getPriamryOutName(), outPairsMateUnmapped.getSecondaryOutName());
//	ret.inThrownAwayUnmappedMate_ = SeqIOOptions::genFastqIn(thrownAwayMate.getPriamryOutName());
//
//
//
//	BamTools::BamReader bReader;
//	bReader.Open(opts.firstName_.string());
//	checkBamOpenThrow(bReader, opts.firstName_.string());
//	auto rData = bReader.GetReferenceData();
//
//	BamTools::BamWriter bWriter;
//	bWriter.Open(njh::files::prependFileBasename(bfs::path(opts.out_.outFilename_).replace_extension(".bam"), "mapped_").string(), bReader.GetHeader(), rData);
//
//
//	BamTools::BamAlignment bAln;
//	BamAlnsCache alnCache;
//
//	OutOptions inverseBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"), ".bed");
//	inverseBedOpts.transferOverwriteOpts(opts.out_);
//
//	OutOptions discordantBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"), ".bed");
//	discordantBedOpts.transferOverwriteOpts(opts.out_);
//
//	OutOptions mappedBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "mapped_"), ".bed");
//	mappedBedOpts.transferOverwriteOpts(opts.out_);
//	std::unique_ptr<OutputStream> inverseBedOut;
//	std::unique_ptr<OutputStream> discordantBedOut;
//	std::unique_ptr<OutputStream> mappedBedOut;
//	if (debug_) {
//		inverseBedOut = std::make_unique<OutputStream>(inverseBedOpts);
//		discordantBedOut = std::make_unique<OutputStream>(discordantBedOpts);
//		mappedBedOut = std::make_unique<OutputStream>(mappedBedOpts);
//	}
//
//	uint32_t centerClipCutOff = 10;
//	std::unordered_map<std::string, uint32_t> contigLengths;
//	for(const auto & r : rData){
//		contigLengths[r.RefName] = r.RefLength;
//	}
//	while (bReader.GetNextAlignment(bAln)) {
//		if (!bAln.IsPrimaryAlignment()) {
//			continue;
//		}
//		if (!bAln.IsPaired()) {
//			if (!bAln.IsMapped()) {
//				++ret.unpairedUnMapped_;
//				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
//			} else {
//				if((0 != bAln.Position
//						&& bAln.CigarData.front().Type == 'S'
//						&& bAln.CigarData.front().Length > centerClipCutOff) ||
//						(contigLengths[rData[bAln.RefID].RefName] != getEndPosition(bAln)
//						&& bAln.CigarData.back().Type == 'S'
//						&& bAln.CigarData.back().Length > centerClipCutOff)) {
//					++ret.unpairedUnMapped_;
//					unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
//				}else{
//					++ret.unpaiedReads_;
//					mappedSinglesWriter.openWrite(
//							seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
//									SangerQualOffset));
//					bWriter.SaveAlignment(bAln);
//				}
//			}
//		} else {
//			if (!alnCache.has(bAln.Name)) {
//				//pair hasn't been added to cache yet so add to cache
//				//this only works if mate and first mate have the same name
//				alnCache.add(bAln);
//				continue;
//			} else {
//				auto search = alnCache.get(bAln.Name);
//				bool balnMapped = bAln.IsMapped();
//				bool searchMapped = search->IsMapped();
//				if(balnMapped &&
//						((0 != bAln.Position
//						&& bAln.CigarData.front().Type == 'S'
//						&& bAln.CigarData.front().Length > centerClipCutOff) ||
//						(contigLengths[rData[bAln.RefID].RefName] != getEndPosition(bAln)
//						&& bAln.CigarData.back().Type == 'S'
//						&& bAln.CigarData.back().Length > centerClipCutOff))){
//					balnMapped = false;
//				}
//
//				if(searchMapped &&
//						((0 != search->Position
//						&& search->CigarData.front().Type == 'S'
//						&& search->CigarData.front().Length > centerClipCutOff) ||
//						(contigLengths[rData[search->RefID].RefName] != search->GetEndPosition()
//						&& search->CigarData.back().Type == 'S'
//						&& search->CigarData.back().Length > centerClipCutOff))){
//					searchMapped = false;
//				}
//
//				if (!balnMapped && !searchMapped) {
//					++ret.pairsUnMapped_;
//					if (bAln.IsFirstMate()) {
//						unmappedPairWriter.openWrite(
//								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
//										false));
//
//					} else {
//						unmappedPairWriter.openWrite(
//								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
//										false));
//					}
//				} else {
//					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
//							SangerQualOffset);
//					seqInfo searchSeq(search->Name, search->QueryBases,
//							search->Qualities, SangerQualOffset);
//
//					if (balnMapped && searchMapped) {
//						// test for inverse
//						if(bAln.IsReverseStrand() == search->IsReverseStrand()){
//							//inverse mates will there be written in technically the wrong orientation to each other but in the reference orientation
//							++ret.inverse_;
//							if (bAln.IsFirstMate()) {
//								inversePairWriter.openWrite(
//										PairedRead(bAlnSeq, searchSeq, false));
//							} else {
//								inversePairWriter.openWrite(
//										PairedRead(searchSeq, bAlnSeq, false));
//							}
//							if(debug_){
//								(*inverseBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
//								(*inverseBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
//							}
//
//						}else{
//							//test for concordant
//							if (bAln.RefID == search->RefID
//									&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
//								++ret.pairedReads_;
//								if (bAln.IsFirstMate()) {
//									mappedPairWriter.openWrite(
//											PairedRead(bAlnSeq, searchSeq, false));
//								} else {
//									mappedPairWriter.openWrite(
//											PairedRead(searchSeq, bAlnSeq, false));
//								}
//								bWriter.SaveAlignment(bAln);
//								bWriter.SaveAlignment(*search);
//								if(debug_){
//									(*mappedBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
//									(*mappedBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
//								}
//							} else {
//								if(debug_){
//									(*discordantBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
//									(*discordantBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
//								}
//								//discordant if mapping to different chromosome or very far away from each other
//								++ret.discordant_;
//								if (bAln.IsFirstMate()) {
//									discordantPairWriter.openWrite(
//											PairedRead(bAlnSeq, searchSeq, false));
//								} else {
//									discordantPairWriter.openWrite(
//											PairedRead(searchSeq, bAlnSeq, false));
//								}
//								bWriter.SaveAlignment(bAln);
//								bWriter.SaveAlignment(*search);
//							}
//						}
//					} else {
//						++ret.pairedReadsMateUnmapped_;
//						//when mate is unmapped the same operation is done to the mate as the other mate, not sure why, so to fix orientation (or at least keep the read in the orientation of it's mate)
//						if (balnMapped && !search->IsMapped()) {
//							searchSeq.reverseComplementRead(false, true);
//						} else if (!balnMapped && search->IsMapped()) {
//							bAlnSeq.reverseComplementRead(false, true);
//						} else {
//							//this shouldn't be happening....
//						}
//						if (throwAwayUnmmpaedMates) {
//							if (balnMapped && !searchMapped) {
//								mappedSinglesWriter.openWrite(bAlnSeq);
//								thrownAwayMateWriter.openWrite(searchSeq);
//								bWriter.SaveAlignment(bAln);
//							} else if (!balnMapped && searchMapped) {
//								mappedSinglesWriter.openWrite(searchSeq);
//								thrownAwayMateWriter.openWrite(bAlnSeq);
//								bWriter.SaveAlignment(*search);
//							}
//						}else{
//							if (bAln.IsFirstMate()) {
//								mateUnmappedPairWriter.openWrite(
//										PairedRead(bAlnSeq, searchSeq, false));
//							} else {
//								mateUnmappedPairWriter.openWrite(
//										PairedRead(searchSeq, bAlnSeq, false));
//							}
//							bWriter.SaveAlignment(bAln);
//							bWriter.SaveAlignment(*search);
//						}
//					}
//				}
//				// now that operations have been computed, remove ther other mate found from cache
//				alnCache.remove(search->Name);
//				continue;
//			}
//		}
//	}
//
//	//save the orphans;
//	if (len(alnCache) > 0) {
//		auto names = alnCache.getNames();
//		for (const auto & name : names) {
//			auto search = alnCache.get(name);
//			if (!search->IsMapped()) {
//				++ret.unpairedUnMapped_;
//				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
//			} else {
//				if((0 != search->Position
//						&& search->CigarData.front().Type == 'S'
//						&& search->CigarData.front().Length > centerClipCutOff) ||
//						(contigLengths[rData[search->RefID].RefName] != getEndPosition(*search)
//						&& search->CigarData.back().Type == 'S'
//						&& search->CigarData.back().Length > centerClipCutOff)){
//					++ret.unpairedUnMapped_;
//					unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
//				}else{
//					++ret.unpaiedReads_;
//					mappedSinglesWriter.openWrite(
//							seqInfo(search->Name, search->QueryBases, search->Qualities,
//									SangerQualOffset));
//					bWriter.SaveAlignment(*search);
//				}
//			}
//			alnCache.remove(name);
//		}
//	}
//	if(verbose_){
//		ret.log(std::cout, opts.firstName_);
//	}
//	return ret;
//}
//


//without center soft clipping filtering
//BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsFromBamToSameOrientationContigs(
//		const SeqIOOptions & opts,
//		bool throwAwayUnmmpaedMates) {
//
//	auto outPairsUnMapped = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
//	outPairsUnMapped.out_.transferOverwriteOpts(opts.out_);
//	auto outUnpairedUnMapped = SeqIOOptions::genFastqOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
//	outUnpairedUnMapped.out_.transferOverwriteOpts(opts.out_);
//
//	auto outPairs = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
//	outPairs.out_.transferOverwriteOpts(opts.out_);
//	auto outUnpaired = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
//	outUnpaired.out_.transferOverwriteOpts(opts.out_);
//
//	auto outPairsMateUnmapped = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "mateUnmapped_"));
//	outPairsMateUnmapped.out_.transferOverwriteOpts(opts.out_);
//
//	auto thrownAwayMate = SeqIOOptions::genFastqOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "thrownAwayMate_"));
//	thrownAwayMate.out_.transferOverwriteOpts(opts.out_);
//
//
//	auto outPairsInverse = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"));
//	outPairsInverse.out_.transferOverwriteOpts(opts.out_);
//
//	auto outPairsDiscordant = SeqIOOptions::genPairedOut(
//			njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"));
//	outPairsDiscordant.out_.transferOverwriteOpts(opts.out_);
//
//	SeqOutput unmappedPairWriter(outPairsUnMapped);
//	SeqOutput unmappedSinglesWriter(outUnpairedUnMapped);
//
//	SeqOutput inversePairWriter(outPairsInverse);
//
//	SeqOutput discordantPairWriter(outPairsDiscordant);
//
//	SeqOutput mappedPairWriter(outPairs);
//	SeqOutput mappedSinglesWriter(outUnpaired);
//
//	SeqOutput mateUnmappedPairWriter(outPairsMateUnmapped);
//	SeqOutput thrownAwayMateWriter(thrownAwayMate);
//
//
//
//	BamExtractor::ExtractedFilesOpts ret;
//	ret.inPairs_ = SeqIOOptions::genPairedIn(outPairs.getPriamryOutName(), outPairs.getSecondaryOutName());
//	ret.inPairsUnMapped_ = SeqIOOptions::genPairedIn(outPairsUnMapped.getPriamryOutName(), outPairsUnMapped.getSecondaryOutName());
//
//	ret.inUnpaired_ = SeqIOOptions::genFastqIn(outUnpaired.getPriamryOutName());
//	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped.getPriamryOutName());
//
//	ret.inDiscordant_ = SeqIOOptions::genPairedIn(outPairsDiscordant.getPriamryOutName(), outPairsDiscordant.getSecondaryOutName());
//	ret.inInverse_ = SeqIOOptions::genPairedIn(outPairsInverse.getPriamryOutName(), outPairsInverse.getSecondaryOutName());
//
//	ret.inPairsMateUnmapped_ = SeqIOOptions::genPairedIn(outPairsMateUnmapped.getPriamryOutName(), outPairsMateUnmapped.getSecondaryOutName());
//	ret.inThrownAwayUnmappedMate_ = SeqIOOptions::genFastqIn(thrownAwayMate.getPriamryOutName());
//
//
//
//	BamTools::BamReader bReader;
//	bReader.Open(opts.firstName_.string());
//	checkBamOpenThrow(bReader, opts.firstName_.string());
//	auto rData = bReader.GetReferenceData();
//
//	BamTools::BamWriter bWriter;
//	bWriter.Open(njh::files::prependFileBasename(bfs::path(opts.out_.outFilename_).replace_extension(".bam"), "mapped_").string(), bReader.GetHeader(), rData);
//
//
//	BamTools::BamAlignment bAln;
//	BamAlnsCache alnCache;
//
//	OutOptions inverseBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"), ".bed");
//	inverseBedOpts.transferOverwriteOpts(opts.out_);
//
//	OutOptions discordantBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"), ".bed");
//	discordantBedOpts.transferOverwriteOpts(opts.out_);
//
//	OutOptions mappedBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "mapped_"), ".bed");
//	mappedBedOpts.transferOverwriteOpts(opts.out_);
//	std::unique_ptr<OutputStream> inverseBedOut;
//	std::unique_ptr<OutputStream> discordantBedOut;
//	std::unique_ptr<OutputStream> mappedBedOut;
//	if (debug_) {
//		inverseBedOut = std::make_unique<OutputStream>(inverseBedOpts);
//		discordantBedOut = std::make_unique<OutputStream>(discordantBedOpts);
//		mappedBedOut = std::make_unique<OutputStream>(mappedBedOpts);
//	}
//
//	//uint32_t centerClipCutOff = 10;
//	std::unordered_map<std::string, uint32_t> contigLengths;
//	for(const auto & r : rData){
//		contigLengths[r.RefName] = r.RefLength;
//	}
//	while (bReader.GetNextAlignment(bAln)) {
//		if (!bAln.IsPrimaryAlignment()) {
//			continue;
//		}
//		if (!bAln.IsPaired()) {
//			if (!bAln.IsMapped()) {
//				++ret.unpairedUnMapped_;
//				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
//			} else {
//				++ret.unpaiedReads_;
//				mappedSinglesWriter.openWrite(
//						seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
//								SangerQualOffset));
//				bWriter.SaveAlignment(bAln);
//			}
//		} else {
//			if (!alnCache.has(bAln.Name)) {
//				//pair hasn't been added to cache yet so add to cache
//				//this only works if mate and first mate have the same name
//				alnCache.add(bAln);
//				continue;
//			} else {
//				auto search = alnCache.get(bAln.Name);
//				bool balnMapped = bAln.IsMapped();
//				bool searchMapped = search->IsMapped();
//
//				if (!balnMapped && !searchMapped) {
//					++ret.pairsUnMapped_;
//					if (bAln.IsFirstMate()) {
//						unmappedPairWriter.openWrite(
//								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
//										false));
//
//					} else {
//						unmappedPairWriter.openWrite(
//								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
//										false));
//					}
//				} else {
//					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
//							SangerQualOffset);
//					seqInfo searchSeq(search->Name, search->QueryBases,
//							search->Qualities, SangerQualOffset);
//
//					if (balnMapped && searchMapped) {
//						// test for inverse
//						if(bAln.IsReverseStrand() == search->IsReverseStrand()){
//							//inverse mates will there be written in technically the wrong orientation to each other but in the reference orientation
//							++ret.inverse_;
//							if (bAln.IsFirstMate()) {
//								inversePairWriter.openWrite(
//										PairedRead(bAlnSeq, searchSeq, false));
//							} else {
//								inversePairWriter.openWrite(
//										PairedRead(searchSeq, bAlnSeq, false));
//							}
//							if(debug_){
//								(*inverseBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
//								(*inverseBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
//							}
//
//						}else{
//							//test for concordant
//							if (bAln.RefID == search->RefID
//									&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
//								++ret.pairedReads_;
//								if (bAln.IsFirstMate()) {
//									mappedPairWriter.openWrite(
//											PairedRead(bAlnSeq, searchSeq, false));
//								} else {
//									mappedPairWriter.openWrite(
//											PairedRead(searchSeq, bAlnSeq, false));
//								}
//								bWriter.SaveAlignment(bAln);
//								bWriter.SaveAlignment(*search);
//								if(debug_){
//									(*mappedBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
//									(*mappedBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
//								}
//							} else {
//								if(debug_){
//									(*discordantBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
//									(*discordantBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
//								}
//								//discordant if mapping to different chromosome or very far away from each other
//								++ret.discordant_;
//								if (bAln.IsFirstMate()) {
//									discordantPairWriter.openWrite(
//											PairedRead(bAlnSeq, searchSeq, false));
//								} else {
//									discordantPairWriter.openWrite(
//											PairedRead(searchSeq, bAlnSeq, false));
//								}
//								bWriter.SaveAlignment(bAln);
//								bWriter.SaveAlignment(*search);
//							}
//						}
//					} else {
//						++ret.pairedReadsMateUnmapped_;
//						//when mate is unmapped the same operation is done to the mate as the other mate, not sure why, so to fix orientation (or at least keep the read in the orientation of it's mate)
//						if (balnMapped && !search->IsMapped()) {
//							searchSeq.reverseComplementRead(false, true);
//						} else if (!balnMapped && search->IsMapped()) {
//							bAlnSeq.reverseComplementRead(false, true);
//						} else {
//							//this shouldn't be happening....
//						}
//						if (throwAwayUnmmpaedMates) {
//							if (balnMapped && !searchMapped) {
//								mappedSinglesWriter.openWrite(bAlnSeq);
//								thrownAwayMateWriter.openWrite(searchSeq);
//								bWriter.SaveAlignment(bAln);
//							} else if (!balnMapped && searchMapped) {
//								mappedSinglesWriter.openWrite(searchSeq);
//								thrownAwayMateWriter.openWrite(bAlnSeq);
//								bWriter.SaveAlignment(*search);
//							}
//						}else{
//							if (bAln.IsFirstMate()) {
//								mateUnmappedPairWriter.openWrite(
//										PairedRead(bAlnSeq, searchSeq, false));
//							} else {
//								mateUnmappedPairWriter.openWrite(
//										PairedRead(searchSeq, bAlnSeq, false));
//							}
//							bWriter.SaveAlignment(bAln);
//							bWriter.SaveAlignment(*search);
//						}
//					}
//				}
//				// now that operations have been computed, remove ther other mate found from cache
//				alnCache.remove(search->Name);
//				continue;
//			}
//		}
//	}
//
//	//save the orphans;
//	if (len(alnCache) > 0) {
//		auto names = alnCache.getNames();
//		for (const auto & name : names) {
//			auto search = alnCache.get(name);
//			if (!search->IsMapped()) {
//				++ret.unpairedUnMapped_;
//				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
//			} else {
//				++ret.unpaiedReads_;
//				mappedSinglesWriter.openWrite(
//						seqInfo(search->Name, search->QueryBases, search->Qualities,
//								SangerQualOffset));
//				bWriter.SaveAlignment(*search);
//			}
//			alnCache.remove(name);
//		}
//	}
//	if(verbose_){
//		ret.log(std::cout, opts.firstName_);
//	}
//	return ret;
//}

BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsFromBamToSameOrientationContigs(
		const SeqIOOptions & opts,
		const extractReadsFromBamToSameOrientationContigsPars & pars) {

	auto outPairsUnMapped = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outPairsUnMapped.out_.transferOverwriteOpts(opts.out_);
	auto outUnpairedUnMapped = SeqIOOptions::genFastqOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outUnpairedUnMapped.out_.transferOverwriteOpts(opts.out_);

	auto outPairs = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
	outPairs.out_.transferOverwriteOpts(opts.out_);
	auto outUnpaired = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(opts.out_);

	auto outPairsMateUnmapped = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "mateUnmapped_"));
	outPairsMateUnmapped.out_.transferOverwriteOpts(opts.out_);

	auto thrownAwayMate = SeqIOOptions::genFastqOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "thrownAwayMate_"));
	thrownAwayMate.out_.transferOverwriteOpts(opts.out_);


	auto outPairsInverse = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"));
	outPairsInverse.out_.transferOverwriteOpts(opts.out_);

	auto outPairsDiscordant = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"));
	outPairsDiscordant.out_.transferOverwriteOpts(opts.out_);

	SeqOutput unmappedPairWriter(outPairsUnMapped);
	SeqOutput unmappedSinglesWriter(outUnpairedUnMapped);

	SeqOutput inversePairWriter(outPairsInverse);

	SeqOutput discordantPairWriter(outPairsDiscordant);

	SeqOutput mappedPairWriter(outPairs);
	SeqOutput mappedSinglesWriter(outUnpaired);

	SeqOutput mateUnmappedPairWriter(outPairsMateUnmapped);
	SeqOutput thrownAwayMateWriter(thrownAwayMate);

//	auto debugPairsOpts = SeqIOOptions::genPairedOut(njh::files::prependFileBasename(opts.out_.outFilename_, "debug_pairs_"));;
//	debugPairsOpts.out_.overWriteFile_ = true;
//	SeqOutput debugPairWriter(debugPairsOpts);
//
//	auto debugSinglsOpts = SeqIOOptions::genFastqOut(njh::files::prependFileBasename(opts.out_.outFilename_, "debug_singles_"));
//	debugSinglsOpts.out_.overWriteFile_ = true;
//	SeqOutput debugSinglesWriter(debugSinglsOpts);



	BamExtractor::ExtractedFilesOpts ret;
	ret.inPairs_ = SeqIOOptions::genPairedIn(outPairs.getPriamryOutName(), outPairs.getSecondaryOutName());
	ret.inPairsUnMapped_ = SeqIOOptions::genPairedIn(outPairsUnMapped.getPriamryOutName(), outPairsUnMapped.getSecondaryOutName());

	ret.inUnpaired_ = SeqIOOptions::genFastqIn(outUnpaired.getPriamryOutName());
	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped.getPriamryOutName());

	ret.inDiscordant_ = SeqIOOptions::genPairedIn(outPairsDiscordant.getPriamryOutName(), outPairsDiscordant.getSecondaryOutName());
	ret.inInverse_ = SeqIOOptions::genPairedIn(outPairsInverse.getPriamryOutName(), outPairsInverse.getSecondaryOutName());

	ret.inPairsMateUnmapped_ = SeqIOOptions::genPairedIn(outPairsMateUnmapped.getPriamryOutName(), outPairsMateUnmapped.getSecondaryOutName());
	ret.inThrownAwayUnmappedMate_ = SeqIOOptions::genFastqIn(thrownAwayMate.getPriamryOutName());



	BamTools::BamReader bReader;
	bReader.Open(opts.firstName_.string());
	checkBamOpenThrow(bReader, opts.firstName_.string());
	auto rData = bReader.GetReferenceData();

	BamTools::BamWriter bWriter;
	bWriter.Open(njh::files::prependFileBasename(bfs::path(opts.out_.outFilename_).replace_extension(".bam"), "mapped_").string(), bReader.GetHeader(), rData);


	BamTools::BamAlignment bAln;
	BamAlnsCache alnCache;

	OutOptions inverseBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"), ".bed");
	inverseBedOpts.transferOverwriteOpts(opts.out_);

	OutOptions discordantBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"), ".bed");
	discordantBedOpts.transferOverwriteOpts(opts.out_);

	OutOptions mappedBedOpts(njh::files::prependFileBasename(opts.out_.outFilename_, "mapped_"), ".bed");
	mappedBedOpts.transferOverwriteOpts(opts.out_);
	std::unique_ptr<OutputStream> inverseBedOut;
	std::unique_ptr<OutputStream> discordantBedOut;
	std::unique_ptr<OutputStream> mappedBedOut;
	if (debug_) {
		inverseBedOut = std::make_unique<OutputStream>(inverseBedOpts);
		discordantBedOut = std::make_unique<OutputStream>(discordantBedOpts);
		mappedBedOut = std::make_unique<OutputStream>(mappedBedOpts);
	}


	std::unordered_map<std::string, uint32_t> contigLengths;
	for(const auto & r : rData){
		contigLengths[r.RefName] = r.RefLength;
	}
	while (bReader.GetNextAlignment(bAln)) {
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (!bAln.IsPaired()) {
			if (!bAln.IsMapped()) {
				++ret.unpairedUnMapped_;
				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
			} else {
//				if((0 != bAln.Position
//						&& bAln.CigarData.front().Type == 'S'
//						&& bAln.CigarData.front().Length > centerClipCutOff) ||
//						(contigLengths[rData[bAln.RefID].RefName] != getEndPosition(bAln)
//						&& bAln.CigarData.back().Type == 'S'
//						&& bAln.CigarData.back().Length > centerClipCutOff)) {
				if(
						((bAln.Position > pars.forSoftClipFilterDistanceToEdges
						&& bAln.CigarData.front().Type == 'S'
						&& bAln.CigarData.front().Length > pars.centerClipCutOff) ||
						(contigLengths[rData[bAln.RefID].RefName] - getEndPosition(bAln) > pars.forSoftClipFilterDistanceToEdges
						&& bAln.CigarData.back().Type == 'S'
						&& bAln.CigarData.back().Length > pars.centerClipCutOff))
//						(bAln.Position > forSoftClipFilterDistanceToEdges &&
//												contigLengths[rData[bAln.RefID].RefName] - getEndPosition(bAln) > forSoftClipFilterDistanceToEdges)
//												&& ((bAln.Position > forSoftClipFilterDistanceToEdges
//												&& bAln.CigarData.front().Type == 'S'
//												&& bAln.CigarData.front().Length > centerClipCutOff) ||
//												(contigLengths[rData[bAln.RefID].RefName] - getEndPosition(bAln) > forSoftClipFilterDistanceToEdges
//												&& bAln.CigarData.back().Type == 'S'
//												&& bAln.CigarData.back().Length > centerClipCutOff))
						) {
					++ret.unpairedFailedSoftClip_;
					unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
//					debugSinglesWriter.openWrite(
//							seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
//									SangerQualOffset));
				}else{
					++ret.unpaiedReads_;
					mappedSinglesWriter.openWrite(
							seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
									SangerQualOffset));
					bWriter.SaveAlignment(bAln);
				}
			}
		} else {
			if (!alnCache.has(bAln.Name)) {
				//pair hasn't been added to cache yet so add to cache
				//this only works if mate and first mate have the same name
				alnCache.add(bAln);
				continue;
			} else {
				auto search = alnCache.get(bAln.Name);
				bool balnIsMapped = bAln.IsMapped();
				bool searchIsMapped = search->IsMapped();
				bool balnPassSoftClipTest = true;
				if(balnIsMapped &&
//						(bAln.Position > forSoftClipFilterDistanceToEdges &&
//												contigLengths[rData[bAln.RefID].RefName] - getEndPosition(bAln) > forSoftClipFilterDistanceToEdges)
//												&& ((bAln.Position > forSoftClipFilterDistanceToEdges
//												&& bAln.CigarData.front().Type == 'S'
//												&& bAln.CigarData.front().Length > centerClipCutOff) ||
//												(contigLengths[rData[bAln.RefID].RefName] - getEndPosition(bAln) > forSoftClipFilterDistanceToEdges
//												&& bAln.CigarData.back().Type == 'S'
//												&& bAln.CigarData.back().Length > centerClipCutOff))
							((bAln.Position > pars.forSoftClipFilterDistanceToEdges
													&& bAln.CigarData.front().Type == 'S'
													&& bAln.CigarData.front().Length > pars.centerClipCutOff) ||
													(contigLengths[rData[bAln.RefID].RefName] - getEndPosition(bAln) > pars.forSoftClipFilterDistanceToEdges
													&& bAln.CigarData.back().Type == 'S'
													&& bAln.CigarData.back().Length > pars.centerClipCutOff))
												){
					balnPassSoftClipTest = false;
				}
				bool searchPassSoftClipTest = true;
				if(searchIsMapped &&
//						(search->Position > forSoftClipFilterDistanceToEdges &&
//												contigLengths[rData[search->RefID].RefName] - getEndPosition(*search) > forSoftClipFilterDistanceToEdges)
//												&& ((search->Position > forSoftClipFilterDistanceToEdges
//												&& search->CigarData.front().Type == 'S'
//												&& search->CigarData.front().Length > centerClipCutOff) ||
//												(contigLengths[rData[search->RefID].RefName] - getEndPosition(*search) > forSoftClipFilterDistanceToEdges
//												&& search->CigarData.back().Type == 'S'
//												&& search->CigarData.back().Length > centerClipCutOff))
						((search->Position > pars.forSoftClipFilterDistanceToEdges
						&& search->CigarData.front().Type == 'S'
						&& search->CigarData.front().Length > pars.centerClipCutOff) ||
						(contigLengths[rData[search->RefID].RefName] - getEndPosition(*search) > pars.forSoftClipFilterDistanceToEdges
						&& search->CigarData.back().Type == 'S'
						&& search->CigarData.back().Length > pars.centerClipCutOff))

				) {
					searchPassSoftClipTest = false;
				}

				if (!balnIsMapped && !searchIsMapped) {
					++ret.pairsUnMapped_;
					if (bAln.IsFirstMate()) {
						unmappedPairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));

					} else {
						unmappedPairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
				} else {
					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities, SangerQualOffset);
					seqInfo searchSeq(search->Name, search->QueryBases, search->Qualities, SangerQualOffset);
					if (balnIsMapped && searchIsMapped) {
						//see if they passed
						if(!balnPassSoftClipTest || ! searchPassSoftClipTest){
							if (!balnPassSoftClipTest && !searchPassSoftClipTest) {
								++ret.pairedReadsBothFailedSoftClip_;
								//treat this like both are unammped
								if (bAln.IsFirstMate()) {
									unmappedPairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
													false));
								} else {
									unmappedPairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
													false));
								}
//								if (bAln.IsFirstMate()) {
//									debugPairWriter.openWrite(
//											PairedRead(bAlnSeq, searchSeq, false));
//								} else {
//									debugPairWriter.openWrite(
//											PairedRead(searchSeq, bAlnSeq, false));
//								}
							} else if (!balnPassSoftClipTest && searchPassSoftClipTest) {
								++ret.pairedReadsMateFailedSoftClip_;
								//treat this like baln is unmapped
								if (pars.throwAwayUnmmpaedMates) {
									mappedSinglesWriter.openWrite(searchSeq);
									thrownAwayMateWriter.openWrite(bAlnSeq);
									bWriter.SaveAlignment(*search);
								} else {
									//make sure pairs are oriented in the right direction
									//if they have both been reverse complemented,
									//then the mate that didn't pass needs to be reverse complemented again
									if(search->IsReverseStrand() == bAln.IsReverseStrand()){
										bAlnSeq.reverseComplementRead(false, true);
									}
									if (bAln.IsFirstMate()) {
										mateUnmappedPairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										mateUnmappedPairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
									bWriter.SaveAlignment(bAln);
									bWriter.SaveAlignment(*search);
								}
							} else if (balnPassSoftClipTest && !searchPassSoftClipTest) {
								++ret.pairedReadsMateFailedSoftClip_;
								//treat this like search is unmapped
								if (pars.throwAwayUnmmpaedMates) {
									mappedSinglesWriter.openWrite(bAlnSeq);
									thrownAwayMateWriter.openWrite(searchSeq);
									bWriter.SaveAlignment(bAln);
								}else{
									//make sure pairs are oriented in the right direction
									//if they have both been reverse complemented,
									//then the mate that didn't pass needs to be reverse complemented again
									if(search->IsReverseStrand() == bAln.IsReverseStrand()){
										searchSeq.reverseComplementRead(false, true);
									}
									if (bAln.IsFirstMate()) {
										mateUnmappedPairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										mateUnmappedPairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
									bWriter.SaveAlignment(bAln);
									bWriter.SaveAlignment(*search);
								}
							} else {
								//this shouldn't be happening...
							}
						} else {
							// test for inverse
							if(bAln.IsReverseStrand() == search->IsReverseStrand()){
								//inverse mates will there be written in technically the wrong orientation to each other but in the reference orientation
								/**@todo lol, why? check on this */
								++ret.inverse_;
//								if (bAln.IsFirstMate()) {
//									inversePairWriter.openWrite(
//											PairedRead(bAlnSeq, searchSeq, false));
//								} else {
//									inversePairWriter.openWrite(
//											PairedRead(searchSeq, bAlnSeq, false));
//								}

								//treat this like both are unammped, this will give them another chance to map correctly if things get extended
								if (bAln.IsFirstMate()) {
									unmappedPairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
													false));

								} else {
									unmappedPairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
													false));
								}
								if(debug_){
									(*inverseBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
									(*inverseBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
								}
							}else{
								//test for concordant
								if (bAln.RefID == search->RefID && std::abs(bAln.InsertSize) < insertLengthCutOff_) {
									++ret.pairedReads_;
									if (bAln.IsFirstMate()) {
										mappedPairWriter.openWrite(PairedRead(bAlnSeq, searchSeq, false));
									} else {
										mappedPairWriter.openWrite(PairedRead(searchSeq, bAlnSeq, false));
									}
									bWriter.SaveAlignment(bAln);
									bWriter.SaveAlignment(*search);
									if(debug_){
										(*mappedBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
										(*mappedBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
									}
								} else {
									if(debug_){
										(*discordantBedOut) << GenomicRegion(bAln,    rData).genBedRecordCore().toDelimStr() << std::endl;
										(*discordantBedOut) << GenomicRegion(*search, rData).genBedRecordCore().toDelimStr() << std::endl;
									}
//									std::cout << std::endl;
//									std::cout << __FILE__ << " " << __LINE__ << std::endl;

//									std::cout << "std::abs(bAln.InsertSize) < insertLengthCutOff_: " << njh::colorBool(std::abs(bAln.InsertSize) < insertLengthCutOff_) << std::endl;
//									std::cout << "std::abs(bAln.InsertSize) : " << std::abs(bAln.InsertSize)  << std::endl;
//									std::cout << "bAln.RefID == search->RefID: " << njh::colorBool(bAln.RefID == search->RefID) << std::endl;
//									std::cout << "bAln.RefID: " << bAln.RefID << std::endl;
//									std::cout << "search->RefID: " << search->RefID << std::endl;
//									std::cout << "bAln.IsReverseStrand(): "  << njh::colorBool(bAln.IsReverseStrand())  << std::endl;
//									std::cout << "search->IsReverseStrand(): "  << njh::colorBool(search->IsReverseStrand())  << std::endl;
//									std::cout << "search->IsReverseStrand() == bAln.IsReverseStrand(): " << njh::colorBool(search->IsReverseStrand() == bAln.IsReverseStrand()) << std::endl;
//									std::cout << "bAln.Position: " << bAln.Position << std::endl;
//									std::cout << "bAln.RefID length: " << rData[bAln.RefID].RefLength << std::endl;
//									std::cout << "search->Position: " << search->Position << std::endl;
//									std::cout << "search->RefID length: " << rData[search->RefID].RefLength << std::endl;
									//discordant if mapping to different chromosome or very far away from each other
									++ret.discordant_;
									if (bAln.IsFirstMate()) {
										discordantPairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										discordantPairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
									bWriter.SaveAlignment(bAln);
									bWriter.SaveAlignment(*search);
								}
							}
						}
					} else {
						bool failedSoftClip = false;
						if(search->IsMapped() && !searchPassSoftClipTest){
							failedSoftClip = true;
						}else if(bAln.IsMapped() && !balnPassSoftClipTest){
							failedSoftClip = true;
						}
						if(failedSoftClip){
							++ret.pairedReadsMateUnmappedFailedSoftClip_;
						} else {
							++ret.pairedReadsMateUnmapped_;
							//when mate is unmapped the same operation is done to the mate as the other mate, not sure why, so to fix orientation (or at least keep the read in the orientation of it's mate)
							//while the above is true, this can always just be handled by checking sure the appropriate actions have been taken
							//make sure pairs are oriented in the right direction
							//if they have both been reverse complemented,
							//then the mate that didn't pass needs to be reverse complemented again
							if (balnIsMapped && !searchIsMapped && bAln.IsReverseStrand() == search->IsReverseStrand()) {
								searchSeq.reverseComplementRead(false, true);
							} else if (!balnIsMapped && searchIsMapped && bAln.IsReverseStrand() == search->IsReverseStrand()) {
								bAlnSeq.reverseComplementRead(false, true);
							}
							if (pars.throwAwayUnmmpaedMates) {
								if (balnIsMapped && !searchIsMapped) {
									mappedSinglesWriter.openWrite(bAlnSeq);
									thrownAwayMateWriter.openWrite(searchSeq);
									bWriter.SaveAlignment(bAln);
								} else if (!balnIsMapped && searchIsMapped) {
									mappedSinglesWriter.openWrite(searchSeq);
									thrownAwayMateWriter.openWrite(bAlnSeq);
									bWriter.SaveAlignment(*search);
								}
							}else{
								if (bAln.IsFirstMate()) {
									mateUnmappedPairWriter.openWrite(PairedRead(bAlnSeq, searchSeq, false));
								} else {
									mateUnmappedPairWriter.openWrite(PairedRead(searchSeq, bAlnSeq, false));
								}
								bWriter.SaveAlignment(bAln);
								bWriter.SaveAlignment(*search);
							}
						}
					}
				}
				// now that operations have been computed, remove ther other mate found from cache
				alnCache.remove(search->Name);
				continue;
			}
		}
	}

	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			if (!search->IsMapped()) {
				++ret.unpairedUnMapped_;
				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
			} else {

				if(
						(search->Position > pars.forSoftClipFilterDistanceToEdges
						&& search->CigarData.front().Type == 'S'
						&& search->CigarData.front().Length > pars.centerClipCutOff) ||
						(contigLengths[rData[search->RefID].RefName] - getEndPosition(*search) > pars.forSoftClipFilterDistanceToEdges
						&& search->CigarData.back().Type == 'S'
						&& search->CigarData.back().Length > pars.centerClipCutOff)
//						(search->Position > forSoftClipFilterDistanceToEdges &&
//												contigLengths[rData[search->RefID].RefName] - getEndPosition(*search) > forSoftClipFilterDistanceToEdges)
//												&& ((search->Position > forSoftClipFilterDistanceToEdges
//												&& search->CigarData.front().Type == 'S'
//												&& search->CigarData.front().Length > centerClipCutOff) ||
//												(contigLengths[rData[search->RefID].RefName] - getEndPosition(*search) > forSoftClipFilterDistanceToEdges
//												&& search->CigarData.back().Type == 'S'
//												&& search->CigarData.back().Length > centerClipCutOff))
				){
					++ret.unpairedFailedSoftClip_;
					unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
//					debugSinglesWriter.openWrite(
//							seqInfo(search->Name, search->QueryBases, search->Qualities,
//									SangerQualOffset));
				}else{
					++ret.unpaiedReads_;
					mappedSinglesWriter.openWrite(
							seqInfo(search->Name, search->QueryBases, search->Qualities,
									SangerQualOffset));
					bWriter.SaveAlignment(*search);
				}
			}
			alnCache.remove(name);
		}
	}
	if(verbose_){
		ret.log(std::cout, opts.firstName_);
	}
	return ret;
}



BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsFromBamWrite(
		const SeqIOOptions & opts,
		bool referenceOrientation,
		bool throwAwayUnmmpaedMates) {

	auto outPairsUnMapped = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outPairsUnMapped.out_.transferOverwriteOpts(opts.out_);
	auto outUnpairedUnMapped = SeqIOOptions::genFastqOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outUnpairedUnMapped.out_.transferOverwriteOpts(opts.out_);

	auto outPairs = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
	outPairs.out_.transferOverwriteOpts(opts.out_);
	auto outUnpaired = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(opts.out_);

	auto outPairsMateUnmapped = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "mateUnmapped_"));
	outPairsMateUnmapped.out_.transferOverwriteOpts(opts.out_);

	auto thrownAwayMate = SeqIOOptions::genFastqOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "thrownAwayMate_"));
	thrownAwayMate.out_.transferOverwriteOpts(opts.out_);


	auto outPairsInverse = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "inverse_"));
	outPairsInverse.out_.transferOverwriteOpts(opts.out_);

	auto outPairsDiscordant = SeqIOOptions::genPairedOut(
			njh::files::prependFileBasename(opts.out_.outFilename_, "discordant_"));
	outPairsDiscordant.out_.transferOverwriteOpts(opts.out_);

	SeqOutput unmappedPairWriter(outPairsUnMapped);
	SeqOutput unmappedSinglesWriter(outUnpairedUnMapped);

	SeqOutput inversePairWriter(outPairsInverse);

	SeqOutput discordantPairWriter(outPairsDiscordant);

	SeqOutput mappedPairWriter(outPairs);
	SeqOutput mappedSinglesWriter(outUnpaired);

	SeqOutput mateUnmappedPairWriter(outPairsMateUnmapped);
	SeqOutput thrownAwayMateWriter(thrownAwayMate);





	BamExtractor::ExtractedFilesOpts ret;
	ret.inPairs_ = SeqIOOptions::genPairedIn(outPairs.getPriamryOutName(), outPairs.getSecondaryOutName());
	ret.inPairsUnMapped_ = SeqIOOptions::genPairedIn(outPairsUnMapped.getPriamryOutName(), outPairsUnMapped.getSecondaryOutName());

	ret.inUnpaired_ = SeqIOOptions::genFastqIn(outUnpaired.getPriamryOutName());
	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped.getPriamryOutName());

	ret.inDiscordant_ = SeqIOOptions::genPairedIn(outPairsDiscordant.getPriamryOutName(), outPairsDiscordant.getSecondaryOutName());
	ret.inInverse_ = SeqIOOptions::genPairedIn(outPairsInverse.getPriamryOutName(), outPairsInverse.getSecondaryOutName());

	ret.inPairsMateUnmapped_ = SeqIOOptions::genPairedIn(outPairsMateUnmapped.getPriamryOutName(), outPairsMateUnmapped.getSecondaryOutName());
	ret.inThrownAwayUnmappedMate_ = SeqIOOptions::genFastqIn(thrownAwayMate.getPriamryOutName());



	BamTools::BamReader bReader;
	bReader.Open(opts.firstName_.string());
	checkBamOpenThrow(bReader, opts.firstName_.string());

	BamTools::BamAlignment bAln;
	BamAlnsCache alnCache;

	while (bReader.GetNextAlignment(bAln)) {
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (!bAln.IsPaired()) {
			if (!bAln.IsMapped()) {
				++ret.unpairedUnMapped_;
				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
			} else {
				++ret.unpaiedReads_;
				if (referenceOrientation) {
					mappedSinglesWriter.openWrite(
							seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
									SangerQualOffset));
				} else {
					mappedSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
				}
			}
		} else {
			if (!alnCache.has(bAln.Name)) {
				//pair hasn't been added to cache yet so add to cache
				//this only works if mate and first mate have the same name
				alnCache.add(bAln);
				continue;
			} else {

				auto search = alnCache.get(bAln.Name);
				if (!bAln.IsMapped() && !bAln.IsMateMapped()) {
					++ret.pairsUnMapped_;
					if (bAln.IsFirstMate()) {
						unmappedPairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						unmappedPairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
				} else {

					if(referenceOrientation){
						seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
								SangerQualOffset);
						seqInfo searchSeq(search->Name, search->QueryBases,
								search->Qualities, SangerQualOffset);

						if (bAln.IsMapped() && search->IsMapped()) {
							//test for concordant
							if (bAln.RefID == search->RefID
									&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
								if(bAln.IsReverseStrand() != search->IsReverseStrand()){
									++ret.pairedReads_;
									if (bAln.IsFirstMate()) {
										mappedPairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										mappedPairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
								}else{
									//inverse mates will there be written in technically the wrong orientation to each other but in the reference orientation
									++ret.inverse_;
									if (bAln.IsFirstMate()) {
										inversePairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										inversePairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
								}
							} else {
								//discordant if mapping to different chromosome or very far away from each other
								++ret.discordant_;
								if (bAln.IsFirstMate()) {
									discordantPairWriter.openWrite(
											PairedRead(bAlnSeq, searchSeq, false));
								} else {
									discordantPairWriter.openWrite(
											PairedRead(searchSeq, bAlnSeq, false));
								}
							}
						} else {
							++ret.pairedReadsMateUnmapped_;
							//when mate is unmapped the same operation is done to the mate as the other mate, not sure why, so to fix orientation (or at least keep the read in the orientation of it's mate)
							if (bAln.IsMapped() && !search->IsMapped()) {
								searchSeq.reverseComplementRead(false, true);
							} else if (!bAln.IsMapped() && search->IsMapped()) {
								bAlnSeq.reverseComplementRead(false, true);
							} else {
								//this shouldn't be happening....
							}
							if (throwAwayUnmmpaedMates) {
								if (bAln.IsMapped() && !search->IsMapped()) {
									mappedSinglesWriter.openWrite(bAlnSeq);
									thrownAwayMateWriter.openWrite(searchSeq);
								} else if (!bAln.IsMapped() && search->IsMapped()) {
									mappedSinglesWriter.openWrite(searchSeq);
									thrownAwayMateWriter.openWrite(bAlnSeq);
								}
							}else{
								if (bAln.IsFirstMate()) {
									mateUnmappedPairWriter.openWrite(
											PairedRead(bAlnSeq, searchSeq, false));
								} else {
									mateUnmappedPairWriter.openWrite(
											PairedRead(searchSeq, bAlnSeq, false));
								}
							}

						}
						//to make it here at least one of alignments had to have mapped
					} else {
						if (bAln.RefID == search->RefID
								&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
							if(bAln.IsReverseStrand() != search->IsReverseStrand()){
								++ret.pairedReads_;
								if (bAln.IsFirstMate()) {
									mappedPairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
													false));
								} else {
									mappedPairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
													false));
								}
							}else{
								++ret.inverse_;
								if (bAln.IsFirstMate()) {
									inversePairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(bAln),
													bamAlnToSeqInfo(*search), false));
								} else {
									inversePairWriter.openWrite(
											PairedRead(bamAlnToSeqInfo(*search),
													bamAlnToSeqInfo(bAln), false));
								}
							}

						} else {
							++ret.discordant_;
							if (bAln.IsFirstMate()) {
								discordantPairWriter.openWrite(
										PairedRead(bamAlnToSeqInfo(bAln),
												bamAlnToSeqInfo(*search), false));
							} else {
								discordantPairWriter.openWrite(
										PairedRead(bamAlnToSeqInfo(*search),
												bamAlnToSeqInfo(bAln), false));
							}
						}
					}
				}
				// now that operations have been computed, remove ther other mate found from cache
				alnCache.remove(search->Name);
				continue;
			}
		}
	}

	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			if (!search->IsMapped()) {
				++ret.unpairedUnMapped_;
				unmappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
			} else {
				++ret.unpaiedReads_;
				if (!referenceOrientation) {
					mappedSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
				} else {
					mappedSinglesWriter.openWrite(
							seqInfo(search->Name, search->QueryBases, search->Qualities,
									SangerQualOffset));
				}
			}
			alnCache.remove(name);
		}
	}
	if(verbose_){
		ret.log(std::cout, opts.firstName_);
	}
	return ret;
}





BamExtractor::ExtractedFilesOpts BamExtractor::writeUnMappedSeqsAndSmallAlnsWithFilters(const SeqIOOptions & opts,
		const ReadCheckerQualCheck & qualChecker,
		const writeUnMappedSeqsAndSmallAlnsWithFiltersPars & alnSizeFilt,
		const std::vector<GenomicRegion> & skipRegions){
	auto outPairsUnMapped_ = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
	outPairsUnMapped_.out_.transferOverwriteOpts(opts.out_);
	auto outUnpairedUnMapped_ = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
	outUnpairedUnMapped_.out_.transferOverwriteOpts(opts.out_);
	SeqOutput pairWriter(outPairsUnMapped_);
	SeqOutput singlesWriter(outUnpairedUnMapped_);
	BamExtractor::ExtractedFilesOpts ret;
	ret.inPairsUnMapped_ = SeqIOOptions::genPairedIn(outPairsUnMapped_.getPriamryOutName(), outPairsUnMapped_.getSecondaryOutName());
	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped_.getPriamryOutName());


	BamTools::BamReader bReader;
	bReader.Open(opts.firstName_.string());
	checkBamOpenThrow(bReader, opts.firstName_.string());

//	if(!bReader.SetRegion(-1, 0, -1, 0)){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", error unable to set region to -1" << "\n";
//		ss << bReader.GetErrorString() << "\n";
//		throw std::runtime_error{ss.str()};
//	}

	BamTools::BamAlignment bAln;
	BamAlnsCache alnCache;

	bfs::path outBam(njh::appendAsNeededRet(opts.out_.outFilename_.string(), ".bam"));
	OutOptions outBamOpts(outBam);
	outBamOpts.transferOverwriteOpts(opts.out_);
	if(outBamOpts.outExists() && !outBamOpts.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outBamOpts.outName() << " already exits" << "\n";
		throw std::runtime_error{ss.str()};
	}

	BamTools::BamWriter bWriter;
	bWriter.Open(outBam.string(), bReader.GetHeader(), bReader.GetReferenceData());
	auto refData = bReader.GetReferenceData();

//	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	while(bReader.GetNextAlignmentCore(bAln)){
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}

		if(!bAln.IsPaired()){
			if(bAln.IsMapped()){
				if (getAlnLen(bAln) <= alnSizeFilt.maxAlnSize_) {
					bool foundInSkipRegion = false;
					for(const auto & region : skipRegions){
						if(region.fallsInThisRegion(bAln, refData)){
							foundInSkipRegion = true;
							break;
						}
					}
					if(!foundInSkipRegion){
						bAln.BuildCharData();
						auto balnSeq = bamAlnToSeqInfo(bAln);
						qualChecker.checkRead(balnSeq);
						bool pass = balnSeq.on_ && len(balnSeq) > alnSizeFilt.minQuerySize_;
						if(alnSizeFilt.filterOffLowEntropyShortAlnsRecruits_){
							kmerInfo kInfo(balnSeq.seq_, alnSizeFilt.entropyKlen_, false);
							if(kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_){
								pass=false;
							}
						}
						if (pass) {
							if(bAln.IsMapped()){
								++ret.unpaiedReads_;
							}else{
								++ret.unpairedUnMapped_;
							}
							singlesWriter.openWrite(balnSeq);
							bWriter.SaveAlignment(bAln);
						} else {
							++ret.singlesFilteredOff_;
						}
					}
				}
			} else {
				bAln.BuildCharData();
				auto balnSeq = bamAlnToSeqInfo(bAln);
				qualChecker.checkRead(balnSeq);
				if (balnSeq.on_) {
					++ret.unpairedUnMapped_;
					singlesWriter.openWrite(balnSeq);
					bWriter.SaveAlignment(bAln);
				} else {
					++ret.singlesFilteredOff_;
				}
			}
		}else{
			if((!bAln.IsMapped()  && !bAln.IsMateMapped()) ||
					(bAln.IsMapped()  && getAlnLen(bAln) <= alnSizeFilt.maxAlnSize_) ||
					(!bAln.IsMapped() && bAln.IsMateMapped() ) ){
				//if both mates are mapped and only one is less than the min aln size then you want find it here
				bAln.BuildCharData();

				if(bAln.IsMapped()){
					bool balnInSkipRegion = false;
					for(const auto & region : skipRegions){
						if(region.fallsInThisRegion(bAln, refData)){
							balnInSkipRegion = true;
							break;
						}
					}
					if(bAln.QueryBases.size() < alnSizeFilt.minQuerySize_){
						balnInSkipRegion = true;
					}
					if(!balnInSkipRegion){
						auto bAlnSeq = bamAlnToSeqInfo(bAln);
						qualChecker.checkRead(bAlnSeq);
						if(!bAlnSeq.on_){
							balnInSkipRegion = true;
						}
					}
					if(balnInSkipRegion){
						if(!bAln.IsMateMapped() && alnCache.has(bAln.Name)){
							alnCache.remove(bAln.Name);
						}
						continue;
					}
				}
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					if(!bAln.IsMapped() && bAln.IsMateMapped() && !bAln.IsFirstMate()){
						continue;
					}
					alnCache.add(bAln);
				} else {
					auto search = alnCache.get(bAln.Name);
					if(!search->IsMapped() && (bAln.IsMapped() && getAlnLen(bAln) > alnSizeFilt.maxAlnSize_)){
						alnCache.remove(search->Name);
						continue;
					}
					bool balnInSkipRegion = false;
					if(bAln.IsMapped()){
						for(const auto & region : skipRegions){
							if(region.fallsInThisRegion(bAln, refData)){
								balnInSkipRegion = true;
								break;
							}
						}
					}
					bool searchInSkipRegion = false;
					if(search->IsMapped()){
						for(const auto & region : skipRegions){
							if(region.fallsInThisRegion(*search, refData)){
								searchInSkipRegion = true;
								break;
							}
						}
					}

					if(!balnInSkipRegion && !searchInSkipRegion){
						auto bAlnSeq = bamAlnToSeqInfo(bAln);
						auto searchSeq = bamAlnToSeqInfo(*search);
						bool balnSizeCheck = bAln.QueryBases.size() > alnSizeFilt.minQuerySize_;
						bool searchSizeCheck = search->QueryBases.size() > alnSizeFilt.minQuerySize_;
						qualChecker.checkRead(bAlnSeq);
						qualChecker.checkRead(searchSeq);
						bool balnCheck = balnSizeCheck && bAlnSeq.on_;
						bool searchCheck = searchSizeCheck && searchSeq.on_;
						if (alnSizeFilt.filterOffLowEntropyShortAlnsRecruits_) {
							if(bAln.IsMapped()){
								kmerInfo kInfo(bAlnSeq.seq_, alnSizeFilt.entropyKlen_, false);
								if (kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_) {
									balnCheck = false;
								}
							}
							if(search->IsMapped()){
								kmerInfo kInfo(searchSeq.seq_, alnSizeFilt.entropyKlen_, false);
								if (kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_) {
									searchCheck = false;
								}
							}
						}
						if (balnCheck && searchCheck) {
							++ret.pairsUnMapped_;
							PairedRead pSeq;
							if (bAln.IsFirstMate()) {
								pSeq = PairedRead(bAlnSeq, searchSeq, false);
							} else {
								pSeq = PairedRead(searchSeq, bAlnSeq, false);
							}
							pairWriter.openWrite(pSeq);
							bWriter.SaveAlignment(bAln);
							bWriter.SaveAlignment(*search);
						} else if (balnCheck) {
							singlesWriter.openWrite(bAlnSeq);
							bWriter.SaveAlignment(bAln);
							++ret.mateFilteredOff_;
						} else if (searchCheck) {
							singlesWriter.openWrite(searchSeq);
							bWriter.SaveAlignment(*search);
							++ret.mateFilteredOff_;
						} else {
							++ret.bothMatesFilteredOff_;
						}
					}
					// now that operations have been computed, remove ther other mate found from cache
					alnCache.remove(search->Name);
				}
			}else{
				if(!bAln.IsMateMapped() && !bAln.IsFirstMate()){
					bAln.BuildCharData();
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << bAln.Name << std::endl;
//					std::cout << "alnCache.has(bAln.Name): " << njh::colorBool(alnCache.has(bAln.Name)) << std::endl;
//					exit(1);
					if(alnCache.has(bAln.Name)){
						alnCache.remove(bAln.Name);
					}
				}
			}
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	//save the orphans
	if (len(alnCache) > 0) {

		auto names = alnCache.getNames();
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(alnSizeFilt.tryToFindOrphanMate_){
			//find orphans' mates if possible;
			BamTools::BamReader bReaderMateFinder;
			bReaderMateFinder.Open(opts.firstName_.string());
			checkBamOpenThrow(bReaderMateFinder, opts.firstName_.string());
			loadBamIndexThrow(bReaderMateFinder);
			auto refData = bReaderMateFinder.GetReferenceData();
			//gather all the orphans regions
			std::unordered_map<std::string, std::set<uint32_t>> orphanPositions;
			for (const auto & name : names) {
				auto search = alnCache.get(name);
				if(search->IsPaired()){
					if(!search->IsMapped() && !search->IsMateMapped()){
//						if(4294967295 == search->MatePosition ){
//							std::cout << njh::bashCT::red;
//						}
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						std::cout << njh::json::toJson(*search) << std::endl;
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						if(4294967295 == search->MatePosition ){
//							std::cout << njh::bashCT::reset;
//
//						}
					} else {
						bool pass = true;
						if(alnSizeFilt.filterOffLowEntropyShortAlnsRecruits_){
							auto searchSeq = bamAlnToSeqInfo(*search);
							kmerInfo kInfo(searchSeq.seq_, alnSizeFilt.entropyKlen_, false);
							pass = kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_;
						}
						if(pass){
							orphanPositions[refData[search->MateRefID].RefName].emplace(search->MatePosition);
						}
					}
//					if (search->IsMateMapped()) {
//						orphanPositions[refData[search->MateRefID].RefName].emplace(search->MatePosition);
//					}
				}
			}
			std::vector<GenomicRegion> orphanMateRegions;
			for(const auto & orPos : orphanPositions){
				for(const auto & pos  : orPos.second){
					orphanMateRegions.emplace_back(GenomicRegion("", orPos.first, pos, pos + 1, false));
				}
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			sortGRegionsByStart(orphanMateRegions);
			for(const auto & reg : orphanMateRegions){
				setBamFileRegionThrow(bReaderMateFinder, reg);
				while (bReaderMateFinder.GetNextAlignmentCore(bAln)) {
					//skip secondary alignments
					if (!bAln.IsPrimaryAlignment()) {
						continue;
					}
					if(static_cast<uint32_t>(bAln.Position) < reg.start_){
						continue;
					}
					if (bAln.IsPaired()) {
						bAln.BuildCharData();
						if (alnCache.has(bAln.Name)) {
							bAln.BuildCharData();

							auto search = alnCache.get(bAln.Name);
							bool balnInSkipRegion = false;
							if(bAln.IsMapped()){
								for(const auto & region : skipRegions){
									if(region.fallsInThisRegion(bAln, refData)){
										balnInSkipRegion = true;
										break;
									}
								}
							}
							bool searchInSkipRegion = false;
							if(search->IsMapped()){
								for(const auto & region : skipRegions){
									if(region.fallsInThisRegion(*search, refData)){
										searchInSkipRegion = true;
										break;
									}
								}
							}

							if(!balnInSkipRegion && !searchInSkipRegion){
								auto bAlnSeq = bamAlnToSeqInfo(bAln);
								auto searchSeq = bamAlnToSeqInfo(*search);
								bool balnSizeCheck = bAln.QueryBases.size() > alnSizeFilt.minQuerySize_;
								bool searchSizeCheck = search->QueryBases.size() > alnSizeFilt.minQuerySize_;
								qualChecker.checkRead(bAlnSeq);
								qualChecker.checkRead(searchSeq);
								bool balnCheck = balnSizeCheck && bAlnSeq.on_;
								bool searchCheck = searchSizeCheck && searchSeq.on_;
								if (alnSizeFilt.filterOffLowEntropyShortAlnsRecruits_) {
									if(bAln.IsMapped()){
										kmerInfo kInfo(bAlnSeq.seq_, alnSizeFilt.entropyKlen_, false);
										if (kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_) {
											balnCheck = false;
										}
									}
									if(search->IsMapped()){
										kmerInfo kInfo(searchSeq.seq_, alnSizeFilt.entropyKlen_, false);
										if (kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_) {
											searchCheck = false;
										}
									}
								}
								if (balnCheck && searchCheck) {
									++ret.pairsUnMapped_;
									PairedRead pSeq;
									if (bAln.IsFirstMate()) {
										pSeq = PairedRead(bAlnSeq, searchSeq, false);
									} else {
										pSeq = PairedRead(searchSeq, bAlnSeq, false);
									}
									pairWriter.openWrite(pSeq);
									bWriter.SaveAlignment(bAln);
									bWriter.SaveAlignment(*search);
								} else if (balnCheck) {
									singlesWriter.openWrite(bAlnSeq);
									bWriter.SaveAlignment(bAln);
									++ret.mateFilteredOff_;
								} else if (searchCheck) {
									singlesWriter.openWrite(searchSeq);
									bWriter.SaveAlignment(*search);
									++ret.mateFilteredOff_;
								} else {
									++ret.bothMatesFilteredOff_;
								}
							}
							// now that operations have been computed, remove ther other mate found from cache
							alnCache.remove(search->Name);
						}
					}
				}
			}
		}
		names = alnCache.getNames();
//		OutOptions orphanNamesOpts(bfs::path("orphansNames.txt"));
//		orphanNamesOpts.overWriteFile_ = true;
//		OutputStream orphanNamesOut(orphanNamesOpts);
//		orphanNamesOut << njh::conToStr(names, "\n") << std::endl;
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			bool searchInSkipRegion = false;
			if(search->IsMapped()){
				for(const auto & region : skipRegions){
					if(region.fallsInThisRegion(*search, refData)){
						searchInSkipRegion = true;
						break;
					}
				}
			}

			if(!searchInSkipRegion){
				bool matePossiblyInSkipRegion = false;
				if(search->IsMateMapped()){
					GenomicRegion possibleMateRegion(search->Name, refData[search->RefID].RefName, search->MatePosition, search->MatePosition + 1, true);
					BedUtility::extendLeftRight(possibleMateRegion, search->QueryBases.size(), search->QueryBases.size() );
					for(const auto & region : skipRegions){
						if(region.fallsInThisRegion(possibleMateRegion)){
							matePossiblyInSkipRegion = true;
							break;
						}
					}
				}
				if(!matePossiblyInSkipRegion){
					auto searchSeq = bamAlnToSeqInfo(*search);

					qualChecker.checkRead(searchSeq);
					bool pass = searchSeq.on_ && len(searchSeq) > alnSizeFilt.minQuerySize_;

					if(alnSizeFilt.filterOffLowEntropyShortAlnsRecruits_){
						kmerInfo kInfo(searchSeq.seq_, alnSizeFilt.entropyKlen_, false);


						if(kInfo.computeKmerEntropy() < alnSizeFilt.filterOffLowEntropyShortAlnsRecruitsCutOff_){
							pass = false;
						}
					}
					if (pass) {
						++ret.orphans_;
						singlesWriter.openWrite(searchSeq);
						bWriter.SaveAlignment(*search);
					} else {
						++ret.orphansFiltered_;
					}
				}
			}
			alnCache.remove(name);
		}
	}
	if(verbose_){
		ret.log(std::cout, opts.firstName_);
	}
	return ret;
}




BamExtractor::ExtractedFilesOpts BamExtractor::writeUnMappedSeqs(const SeqIOOptions & opts){
	auto outPairsUnMapped_ = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
	outPairsUnMapped_.out_.transferOverwriteOpts(opts.out_);
	auto outUnpairedUnMapped_ = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
	outUnpairedUnMapped_.out_.transferOverwriteOpts(opts.out_);
	SeqOutput pairWriter(outPairsUnMapped_);
	SeqOutput singlesWriter(outUnpairedUnMapped_);
	BamExtractor::ExtractedFilesOpts ret;
	ret.inPairsUnMapped_ = SeqIOOptions::genPairedIn(outPairsUnMapped_.getPriamryOutName(), outPairsUnMapped_.getSecondaryOutName());
	ret.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(outUnpairedUnMapped_.getPriamryOutName());


	BamTools::BamReader bReader;
	bReader.Open(opts.firstName_.string());
	checkBamOpenThrow(bReader, opts.firstName_.string());

//	if(!bReader.SetRegion(-1, 0, -1, 0)){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", error unable to set region to -1" << "\n";
//		ss << bReader.GetErrorString() << "\n";
//		throw std::runtime_error{ss.str()};
//	}

	BamTools::BamAlignment bAln;
	BamAlnsCache alnCache;

	bfs::path outBam(njh::appendAsNeededRet(opts.out_.outFilename_.string(), ".bam"));
	OutOptions outBamOpts(outBam);
	outBamOpts.transferOverwriteOpts(opts.out_);
	if(outBamOpts.outExists() && !outBamOpts.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outBamOpts.outName() << " already exits" << "\n";
		throw std::runtime_error{ss.str()};
	}

	BamTools::BamWriter bWriter;
	bWriter.Open(outBam.string(), bReader.GetHeader(), bReader.GetReferenceData());


	while(bReader.GetNextAlignmentCore(bAln)){

		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if(!bAln.IsPaired() && !bAln.IsMapped()){
			bAln.BuildCharData();
			++ret.unpairedUnMapped_;
			singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
			bWriter.SaveAlignment(bAln);
		}else{
			if(!bAln.IsMapped() && !bAln.IsMateMapped()){
				bAln.BuildCharData();
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					++ret.pairsUnMapped_;
					if (bAln.IsFirstMate()) {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),
										false));
					} else {
						pairWriter.openWrite(
								PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),
										false));
					}
					bWriter.SaveAlignment(bAln);
					bWriter.SaveAlignment(*search);
					// now that operations have been computed, remove ther other mate found from cache
					alnCache.remove(search->Name);
					continue;
				}
			}
		}
	}

	if(verbose_){
		ret.log(std::cout, opts.firstName_);
	}
	return ret;
}

std::shared_ptr<GenomicRegion> findFirstRegionOfMate(
		const BamTools::BamAlignment & bAln,
		const std::vector<GenomicRegion> & regions,
		const BamTools::RefVector & refData) {
	//only a theoritical end of the mate position, can't be determine by it's mate
	GenomicRegion matePosition(bAln.Name, refData.at(bAln.MateRefID).RefName,
			bAln.MatePosition, bAln.MatePosition + bAln.QueryBases.size(), bAln.IsMateReverseStrand());
	for (const auto otherRegionPos : iter::range(regions.size())) {
		const auto & otherRegion = regions[otherRegionPos];
		if (otherRegion.startsInThisRegion(matePosition) || otherRegion.endsInThisRegion(matePosition)) {
			return std::make_shared<GenomicRegion>(otherRegion);
		}
	}
	return nullptr;
}

std::shared_ptr<GenomicRegion> findFirstRegion(
		const BamTools::BamAlignment & bAln, const BamTools::RefVector & refData,
		const std::vector<GenomicRegion> & regions, double percInRegion) {
	for (const auto otherRegionPos : iter::range(regions.size())) {
		const auto & otherRegion = regions[otherRegionPos];
		if (otherRegion.getPercInRegion(bAln, refData) >= percInRegion) {
			return std::make_shared<GenomicRegion>(otherRegion);
		}
	}
	return nullptr;
}




BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsWtihCrossRegionMappingAsSingles(
		const SeqIOOptions & inOutOpts,
		const std::vector<GenomicRegion> & regions,
		double percInRegion,
		bool originalOrientation) {
	//check to see if regions overlap
	bool overlapsFound = false;
	std::stringstream overalMessage;
	overalMessage << __PRETTY_FUNCTION__ << " error, found overlaps, overlap regions should be merged" << "\n";
	for(const auto regPos : iter::range(regions.size())){
		for(const auto secondPos : iter::range<size_t>(0, regPos)){
			if(regions[regPos].overlaps(regions[secondPos])){
				overlapsFound = true;
				overalMessage << regions[regPos].genBedRecordCore().toDelimStr() << " over laps " << regions[secondPos].genBedRecordCore().toDelimStr() << "\n";
			}
		}
	}
	if(overlapsFound){
		throw std::runtime_error{overalMessage.str()};
	}
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(inOutOpts.firstName_.string());
	checkBamOpenThrow(bReader, inOutOpts.firstName_.string());
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();

	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	auto outUnpaired = SeqIOOptions::genFastqOut(inOutOpts.out_.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(inOutOpts.out_);

	//non paired writer
	SeqOutput writer(outUnpaired);

	BamExtractor::ExtractedFilesOpts ret;

	ret.inUnpaired_ = SeqIOOptions::genFastqIn(outUnpaired.getPriamryOutName());

	for (const auto regionPos : iter::range(regions.size())) {
		const auto & region = regions[regionPos];
		if (verbose_) {
			std::cout << region.uid_ << std::endl;
		}
		setBamFileRegionThrow(bReader, region);

		while (bReader.GetNextAlignment(bAln)) {
			//skip secondary alignments
			if (!bAln.IsPrimaryAlignment()) {
				continue;
			}
			if (bAln.IsMapped()
					&& region.getPercInRegion(bAln, refData) >= percInRegion) {
				++ret.unpaiedReads_;
				if (originalOrientation) {
					writer.openWrite(bamAlnToSeqInfo(bAln));
				} else {
					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
							SangerQualOffset);
					if (region.reverseSrand_) {
						bAlnSeq.reverseComplementRead(false, true);
					}
					writer.openWrite(bAlnSeq);
				}
			}
		}
	}
	return ret;
}


Json::Value BamExtractor::extractReadsWtihCrossRegionMappingPars::toJson() const{
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["percInRegion_"] = njh::json::toJson(percInRegion_);
	ret["originalOrientation_"] = njh::json::toJson(originalOrientation_);
	ret["throwAwayUnmappedMate_"] = njh::json::toJson(throwAwayUnmappedMate_);
	ret["tryToFindOrphansMate_"] = njh::json::toJson(tryToFindOrphansMate_);
	ret["keepMarkedDuplicate_"] = njh::json::toJson(keepMarkedDuplicate_);
	ret["minAlnMapSize_"] = njh::json::toJson(minAlnMapSize_);
	ret["filterOffLowEntropyOrphansRecruitsCutOff_"] = njh::json::toJson(filterOffLowEntropyOrphansRecruitsCutOff_);
	ret["filterOffLowEntropyOrphansRecruits_"] = njh::json::toJson(filterOffLowEntropyOrphansRecruits_);
	ret["entropyKlen_"] = njh::json::toJson(entropyKlen_);
	ret["softClipPercentageCutOff_"] = njh::json::toJson(softClipPercentageCutOff_);
	ret["removeImproperPairs_"] = njh::json::toJson(removeImproperPairs_);
	ret["keepImproperMateUnmapped_"] = njh::json::toJson(keepImproperMateUnmapped_);
	ret["fivePrimeTrim_"] = njh::json::toJson(fivePrimeTrim_);
	ret["threePrimeTrim_"] = njh::json::toJson(threePrimeTrim_);
	ret["writeAll_"] = njh::json::toJson(writeAll_);
	ret["percentSubSample_"] = njh::json::toJson(percentSubSample_);
	return ret;
}




Json::Value BamExtractor::writeUnMappedSeqsAndSmallAlnsWithFiltersPars::toJson() const{
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["maxAlnSize_"] = njh::json::toJson(maxAlnSize_);
	ret["minSotClip_"] = njh::json::toJson(minSotClip_);
	ret["minQuerySize_"] = njh::json::toJson(minQuerySize_);
	ret["tryToFindOrphanMate_"] = njh::json::toJson(tryToFindOrphanMate_);

	return ret;
}


BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsWtihCrossRegionMapping(
		BamTools::BamReader & bReader,
		const OutOptions & outOpts,
		const std::vector<GenomicRegion> & regions,
		const extractReadsWtihCrossRegionMappingPars & extractPars){




	njh::randomGenerator rGen;
	std::function<bool()> subSamplingFunction;
	if(extractPars.percentSubSample_ == 1){
		subSamplingFunction = [](){
			return true;
		};
	}else{
		subSamplingFunction = [&rGen,&extractPars](){
			return rGen()  <= extractPars.percentSubSample_;
		};
	};
	BamTools::BamAlignment bAln;

	//std::cout << "percInRegion: " << percInRegion << std::endl;
	//check to see if regions overlap
	bool overlapsFound = false;
	std::stringstream overalMessage;
	overalMessage << __PRETTY_FUNCTION__ << " error, found overlaps, overlap regions should be merged" << "\n";
	for(const auto regPos : iter::range(regions.size())){
		for(const auto secondPos : iter::range<size_t>(0, regPos)){
			if(regions[regPos].overlaps(regions[secondPos])){
				overlapsFound = true;
				overalMessage << regions[regPos].genBedRecordCore().toDelimStr() << " over laps " << regions[secondPos].genBedRecordCore().toDelimStr() << "\n";
			}
		}
	}
	if(overlapsFound){
		throw std::runtime_error{overalMessage.str()};
	}
	/**@todo also need to do a check for if a sequence lands in more than 1 region */


	auto refs = bReader.GetReferenceData();
	BamAlnsCacheWithRegion alnCache;

	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	bfs::path outBam(njh::appendAsNeededRet(outOpts.outFilename_.string(), ".bam"));
	BamTools::BamWriter bWriter;
	if(extractPars.writeAll_){
		OutOptions outBamOpts(outBam);
		outBamOpts.transferOverwriteOpts(outOpts);
		if(outBamOpts.outExists() && !outBamOpts.overWriteFile_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << outBamOpts.outName() << " already exits" << "\n";
			throw std::runtime_error{ss.str()};
		}
		bWriter.Open(outBam.string(), bReader.GetHeader(), bReader.GetReferenceData());
	}

	auto outPairs = SeqIOOptions::genPairedOut(outOpts.outFilename_);
	outPairs.out_.transferOverwriteOpts(outOpts);

	auto outPairsUnmappedMate = SeqIOOptions::genPairedOut(njh::files::prependFileBasename(outOpts.outFilename_, "mateUnmapped_"));
	outPairsUnmappedMate.out_.transferOverwriteOpts(outOpts);

	auto thrownAwayUnammpedMateOpts = SeqIOOptions::genFastqOut(njh::files::prependFileBasename(outOpts.outFilename_, "thrownAwayMate_"));
	thrownAwayUnammpedMateOpts.out_.transferOverwriteOpts(outOpts);

	//inverse here is if when the pairs are put into the orientation of the target region of interest they end up being not being in the same orientation or regular inverse as well
	auto outPairsInverse = SeqIOOptions::genPairedOut(njh::files::prependFileBasename(outOpts.outFilename_, "inverse_"));
	outPairsInverse.out_.transferOverwriteOpts(outOpts);

	//filtered off seqs, the inverse and not completely in the region sequences
	auto outPairsFiltered = SeqIOOptions::genPairedOut(njh::files::prependFileBasename(outOpts.outFilename_, "filteredPairs_"));
	outPairsFiltered.out_.transferOverwriteOpts(outOpts);

	auto outPairsFilteredSoftClip = SeqIOOptions::genPairedOut(njh::files::prependFileBasename(outOpts.outFilename_, "filteredSoftClipPairs_"));
	outPairsFilteredSoftClip.out_.transferOverwriteOpts(outOpts);

	//filtered off single seqs
	auto filteredSinglesOpts = SeqIOOptions::genFastqOut(njh::files::prependFileBasename(outOpts.outFilename_, "filteredSingles_"));
	filteredSinglesOpts.out_.transferOverwriteOpts(outOpts);

	//filtered off single seqs due to soft clipping
	auto filteredSoftClipSinglesOpts = SeqIOOptions::genFastqOut(njh::files::prependFileBasename(outOpts.outFilename_, "filteredSoftClipSingles_"));
	filteredSoftClipSinglesOpts.out_.transferOverwriteOpts(outOpts);


	auto outUnpaired = SeqIOOptions::genFastqOut(outOpts.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(outOpts);
	//pair writers
	SeqOutput pairWriter(outPairs);
	SeqOutput filteredPairWriter(outPairsFiltered);
	SeqOutput filteredPairSoftClipWriter(outPairsFilteredSoftClip);

	SeqOutput mateUnmappedPairWriter(outPairsUnmappedMate);

	SeqOutput thrownAwayUnammpedMateWriter(thrownAwayUnammpedMateOpts);
	SeqOutput inversePairWriter(outPairsInverse);
	//non paired writer
	SeqOutput writer(outUnpaired);
	SeqOutput singleFilteredWriter(filteredSinglesOpts);
	SeqOutput singleFilteredSoftClipWriter(filteredSoftClipSinglesOpts);

	BamExtractor::ExtractedFilesOpts ret;

	ret.inPairs_ =                  SeqIOOptions::genPairedIn(outPairs.getPriamryOutName(), outPairs.getSecondaryOutName());
	ret.inPairsMateUnmapped_ =      SeqIOOptions::genPairedIn(outPairsUnmappedMate.getPriamryOutName(), outPairsUnmappedMate.getSecondaryOutName());
	ret.inThrownAwayUnmappedMate_ = SeqIOOptions::genFastqIn (thrownAwayUnammpedMateOpts.getPriamryOutName());
	ret.inInverse_ =                SeqIOOptions::genPairedIn(outPairsInverse.getPriamryOutName(), outPairsInverse.getSecondaryOutName());
	ret.inFilteredPairs_ =          SeqIOOptions::genPairedIn(outPairsFiltered.getPriamryOutName(), outPairsFiltered.getSecondaryOutName());
	ret.inSoftClipFilteredPairs_ =  SeqIOOptions::genPairedIn(outPairsFilteredSoftClip.getPriamryOutName(), outPairsFilteredSoftClip.getSecondaryOutName());
	//no disconcordant reads as this is aiming to grab those reads
	ret.inUnpaired_ =               SeqIOOptions::genFastqIn (outUnpaired.getPriamryOutName());
	ret.inFilteredSingles_ =        SeqIOOptions::genFastqIn (filteredSinglesOpts.getPriamryOutName());
	ret.inSoftClipFilteredSingles_ =SeqIOOptions::genFastqIn (filteredSoftClipSinglesOpts.getPriamryOutName());


	//debuging
	auto debugPairsOpts = SeqIOOptions::genPairedOut(njh::files::prependFileBasename(outOpts.outFilename_, "debug_pairs_"));;
	debugPairsOpts.out_.overWriteFile_ = true;
	SeqOutput debugPairWriter(debugPairsOpts);

	auto debugSinglsOpts = SeqIOOptions::genFastqOut(njh::files::prependFileBasename(outOpts.outFilename_, "debug_singles_"));
	debugSinglsOpts.out_.overWriteFile_ = true;
	SeqOutput debugSinglesWriter(debugSinglsOpts);


	auto writeMateFilteredOff = [&ret,&extractPars,&writer,&bWriter,&refData](const BamTools::BamAlignment & bAln, const GenomicRegion & region){
		//unpaired read
		++ret.mateFilteredOff_;
		if (extractPars.originalOrientation_) {
			auto outSeq = bamAlnToSeqInfo(bAln);
//			if(extractPars.fivePrimeTrim_ > 0 && len(outSeq) > extractPars.fivePrimeTrim_){
//				readVecTrimmer::trimOffForwardBases(outSeq, extractPars.fivePrimeTrim_);
//			}
//			if(extractPars.threePrimeTrim_ > 0 && len(outSeq) > extractPars.threePrimeTrim_){
//				readVecTrimmer::trimOffEndBases(outSeq, extractPars.threePrimeTrim_);
//			}
			writer.openWrite(outSeq);
		} else {
			seqInfo outSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
					SangerQualOffset);
			//put in the orientation of the output region
//			bool revComp = bAln.IsReverseStrand();
//			if(extractPars.fivePrimeTrim_ > 0 && len(outSeq) > extractPars.fivePrimeTrim_){
//				readVecTrimmer::trimOffForwardBases(outSeq, extractPars.fivePrimeTrim_);
//			}
//			if(extractPars.threePrimeTrim_ > 0 && len(outSeq) > extractPars.threePrimeTrim_){
//				readVecTrimmer::trimOffEndBases(outSeq, extractPars.threePrimeTrim_);
//			}
			if((len(outSeq) > region.getLen() || region.getLen() < 150)){
				seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
				GenomicRegion balnRegion(bAln, refData);
				uint32_t startRelative = region.start_ - balnRegion.start_;
				uint32_t endRelative = region.end_ - balnRegion.start_;

				seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
				auto alnInfo = bamAlnToAlnInfoLocal(bAln);
				alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
				alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);

				uint32_t startAln = 0;
				if(region.start_ > balnRegion.start_){
					startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
				}
				uint32_t endAln = len(holderSeq);
				if(region.end_ < balnRegion.end_){
					endAln =  getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
				}

				auto outSeqTrimmed = querySeq.getSubRead(startAln, endAln - startAln);
				outSeqTrimmed.removeGaps();
				outSeq = outSeqTrimmed;
			}
			if(region.reverseSrand_){
				outSeq.reverseComplementRead(false, true);
//				revComp = !revComp;
			}
//			MetaDataInName seqMeta;
//			if(MetaDataInName::nameHasMetaData(outSeq.name_)){
//				seqMeta = MetaDataInName(outSeq.name_);
//			}
//			seqMeta.addMeta("isFirstMate", bAln.IsFirstMate());
//			seqMeta.addMeta("RevComp", revComp);
//			seqMeta.resetMetaInName(outSeq.name_);
			writer.openWrite(outSeq);
		}
		if(extractPars.writeAll_){
			bWriter.SaveAlignment(bAln);
		}
	};
	auto writeMateFilteredOffSoftClip = [&ret,&extractPars,&writer,&bWriter,&refData](const BamTools::BamAlignment & bAln, const GenomicRegion & region){
		//unpaired read
		++ret.pairedReadsMateFailedSoftClip_;
		if (extractPars.originalOrientation_) {
			writer.openWrite(bamAlnToSeqInfo(bAln));
		} else {
			seqInfo outSeq(bAln.Name, bAln.QueryBases, bAln.Qualities, SangerQualOffset);

			if((len(outSeq) > region.getLen() || region.getLen() < 150)){
				seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
				GenomicRegion balnRegion(bAln, refData);
				uint32_t startRelative = region.start_ - balnRegion.start_;
				uint32_t endRelative = region.end_ - balnRegion.start_;

				seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
				auto alnInfo = bamAlnToAlnInfoLocal(bAln);
				alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
				alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);

				uint32_t startAln = 0;
				if(region.start_ > balnRegion.start_){
					startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
				}
				uint32_t endAln = len(holderSeq);
				if(region.end_ < balnRegion.end_){
					endAln =  getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
				}

				auto outSeqTrimmed = querySeq.getSubRead(startAln, endAln - startAln);
				outSeqTrimmed.removeGaps();
				outSeq = outSeqTrimmed;
			}

			//put in the orientation of the output region
			if(region.reverseSrand_){
				outSeq.reverseComplementRead(false, true);
			}
			writer.openWrite(outSeq);
		}
		if(extractPars.writeAll_){
			bWriter.SaveAlignment(bAln);
		}
	};



//	auto writeInversePair = [&ret,&extractPars,&inversePairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
//			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
//		++ret.inverse_;
//		if(extractPars.originalOrientation_){
//			if (bAln.IsFirstMate()) {
//				inversePairWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(searchAln)));
//			} else {
//				inversePairWriter.openWrite(PairedRead(bamAlnToSeqInfo(searchAln), bamAlnToSeqInfo(bAln)));
//			}
//		}else{
//			if (bAln.IsFirstMate()) {
//				inversePairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
//			} else {
//				inversePairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
//			}
//		}
//	if(extractPars.writeAll_){
		//		bWriter.SaveAlignment(bAln);
		//		bWriter.SaveAlignment(searchAln);
//	}
//		bWriter.SaveAlignment(bAln);
//		bWriter.SaveAlignment(searchAln);
//	};

	auto writeInverseFilteredPair = [&ret,&filteredPairWriter,&bWriter,&singleFilteredWriter,&extractPars](const BamTools::BamAlignment & bAln,
			const BamTools::BamAlignment & searchAln){
		++ret.inverse_;
		auto bAlnSeq = bamAlnToSeqInfo(bAln);
		auto searchSeq = bamAlnToSeqInfo(searchAln);

		bool bAlnSeqPass = true;
		bool searchSeqPass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			{
				kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
				if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
					bAlnSeqPass = false;
				}
			}
			{
				kmerInfo kInfo(searchSeq.seq_, extractPars.entropyKlen_, false);
				if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
					searchSeqPass= false;
				}
			}
		}
		if (extractPars.writeAll_) {
			if (bAlnSeqPass && searchSeqPass) {
				if (bAln.IsFirstMate()) {
					filteredPairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
				} else {
					filteredPairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
				}

			} else if (bAlnSeqPass) {
				singleFilteredWriter.openWrite(bAlnSeq);
			} else if (searchSeqPass) {
				singleFilteredWriter.openWrite(searchSeq);
			}
			bWriter.SaveAlignment(bAln);
			bWriter.SaveAlignment(searchAln);
		}
	};
//
	auto writeBothPairsFiltered = [&ret,&extractPars,&filteredPairWriter,&bWriter,&singleFilteredWriter](const BamTools::BamAlignment & bAln,
			const BamTools::BamAlignment & searchAln){
		++ret.bothMatesFilteredOff_;
		auto bAlnSeq = bamAlnToSeqInfo(bAln);
		auto searchSeq = bamAlnToSeqInfo(searchAln);

		bool bAlnSeqPass = true;
		bool searchSeqPass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			{
				kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
				if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
					bAlnSeqPass = false;
				}
			}
			{
				kmerInfo kInfo(searchSeq.seq_, extractPars.entropyKlen_, false);
				if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
					searchSeqPass= false;
				}
			}
		}
		if (extractPars.writeAll_) {
			if (bAlnSeqPass && searchSeqPass) {
				if (bAln.IsFirstMate()) {
					filteredPairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
				} else {
					filteredPairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
				}
			} else if (bAlnSeqPass) {
				singleFilteredWriter.openWrite(bAlnSeq);
			} else if (searchSeqPass) {
				singleFilteredWriter.openWrite(searchSeq);
			}
			bWriter.SaveAlignment(bAln);
			bWriter.SaveAlignment(searchAln);
		}
	};

	auto writeBothPairsFilteredSoftClip = [&ret,&extractPars,&filteredPairSoftClipWriter,&bWriter,&singleFilteredSoftClipWriter](const BamTools::BamAlignment & bAln,
			const BamTools::BamAlignment & searchAln){
		++ret.pairedReadsBothFailedSoftClip_;
		auto bAlnSeq = bamAlnToSeqInfo(bAln);
		auto searchSeq = bamAlnToSeqInfo(searchAln);

		bool bAlnSeqPass = true;
		bool searchSeqPass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			{
				kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
				if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
					bAlnSeqPass = false;
				}
			}
			{
				kmerInfo kInfo(searchSeq.seq_, extractPars.entropyKlen_, false);
				if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
					searchSeqPass= false;
				}
			}
		}
		if (extractPars.writeAll_) {
			if (bAlnSeqPass && searchSeqPass) {
				if (bAln.IsFirstMate()) {
					filteredPairSoftClipWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
				} else {
					filteredPairSoftClipWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
				}
			} else if (bAlnSeqPass) {
				singleFilteredSoftClipWriter.openWrite(bAlnSeq);
			} else if (searchSeqPass) {
				singleFilteredSoftClipWriter.openWrite(searchSeq);
			}
			bWriter.SaveAlignment(bAln);
			bWriter.SaveAlignment(searchAln);
		}
	};

	auto writeSingleFiltered = [&ret,&extractPars,&singleFilteredWriter,&bWriter](const BamTools::BamAlignment & bAln){
		++ret.singlesFilteredOff_;
		auto bAlnSeq = bamAlnToSeqInfo(bAln);

		bool pass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
			if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
				pass = false;
			}
		}
		if (extractPars.writeAll_) {
			if (pass) {
				singleFilteredWriter.openWrite(bAlnSeq);
			}
			bWriter.SaveAlignment(bAln);
		}
	};

	auto writeSingleFilteredSoftClip = [&ret,&extractPars,&singleFilteredSoftClipWriter,&bWriter](const BamTools::BamAlignment & bAln){
		++ret.unpairedFailedSoftClip_;
		auto bAlnSeq = bamAlnToSeqInfo(bAln);

		bool pass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
			if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
				pass = false;
			}
		}
		if (extractPars.writeAll_) {
			if (pass) {
				singleFilteredSoftClipWriter.openWrite(bAlnSeq);
			}
			bWriter.SaveAlignment(bAln);
		}
	};

	auto writeDiscordantPair = [&ret,&extractPars,&pairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.discordant_;
		if(extractPars.originalOrientation_){
			if (bAln.IsFirstMate()) {
				pairWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(searchAln)));
			} else {
				pairWriter.openWrite(PairedRead(bamAlnToSeqInfo(searchAln), bamAlnToSeqInfo(bAln)));
			}
		}else{
			if (bAln.IsFirstMate()) {
				pairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
			} else {
				pairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
			}
		}
		if (extractPars.writeAll_) {
			bWriter.SaveAlignment(bAln);
			bWriter.SaveAlignment(searchAln);
		}
	};

	auto writeRegPair = [&ret,&extractPars,&pairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.pairedReads_;
		if(extractPars.originalOrientation_){
			if (bAln.IsFirstMate()) {
				pairWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(searchAln)));
			} else {
				pairWriter.openWrite(PairedRead(bamAlnToSeqInfo(searchAln), bamAlnToSeqInfo(bAln)));
			}
		} else {
			if (bAln.IsFirstMate()) {
				pairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
			} else {
				pairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
			}
		}
		if (extractPars.writeAll_) {
			bWriter.SaveAlignment(bAln);
			bWriter.SaveAlignment(searchAln);
		}
	};
	auto writeMateUnmappedPair = [&ret,&extractPars,&mateUnmappedPairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.pairedReadsMateUnmapped_;
		if(extractPars.originalOrientation_){
			if (bAln.IsFirstMate()) {
				mateUnmappedPairWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(searchAln)));
			} else {
				mateUnmappedPairWriter.openWrite(PairedRead(bamAlnToSeqInfo(searchAln), bamAlnToSeqInfo(bAln)));
			}
		}else{
			if (bAln.IsFirstMate()) {
				mateUnmappedPairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
			} else {
				mateUnmappedPairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
			}
		}
		if (extractPars.writeAll_) {
			bWriter.SaveAlignment(bAln);
			bWriter.SaveAlignment(searchAln);
		}
	};

	auto writeThrowAwayUnmappedMate = [&ret,&extractPars,&writer,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		++ret.pairedReadsMateUnmapped_;
		if(extractPars.originalOrientation_){
			writer.openWrite(bamAlnToSeqInfo(bAln));
		}else{
			writer.openWrite(bAlnSeq);
		}
		if(extractPars.writeAll_){
			bWriter.SaveAlignment(bAln);
		}
	};

	auto writeTheThrownAwayUnmappedMate = [&extractPars,&thrownAwayUnammpedMateWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		bool pass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
			if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
				pass = false;
			}
		}
		if(pass){
			if(extractPars.originalOrientation_){
				thrownAwayUnammpedMateWriter.openWrite(bamAlnToSeqInfo(bAln));
			}else{
				thrownAwayUnammpedMateWriter.openWrite(bAlnSeq);
			}
		}
		if(extractPars.writeAll_){
			bWriter.SaveAlignment(bAln);
		}
	};

	auto writeTheThrownAwayMate = [&extractPars,&thrownAwayUnammpedMateWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		bool pass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
			if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
				pass = false;
			}
		}
		if(pass){
			if(extractPars.originalOrientation_){
				thrownAwayUnammpedMateWriter.openWrite(bamAlnToSeqInfo(bAln));
			}else{
				thrownAwayUnammpedMateWriter.openWrite(bAlnSeq);
			}
		}
		if(extractPars.writeAll_){
			bWriter.SaveAlignment(bAln);
		}
	};

	auto writeTheMateFailedSoftClip= [&extractPars,&singleFilteredSoftClipWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		bool pass = true;
		if(extractPars.filterOffLowEntropyOrphansRecruits_){
			kmerInfo kInfo(bAlnSeq.seq_, extractPars.entropyKlen_, false);
			if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
				pass = false;
			}
		}
		if (extractPars.writeAll_) {
			if (pass) {
				if (extractPars.originalOrientation_) {
					singleFilteredSoftClipWriter.openWrite(bamAlnToSeqInfo(bAln));
				} else {
					singleFilteredSoftClipWriter.openWrite(bAlnSeq);
				}
			}
			bWriter.SaveAlignment(bAln);
		}
	};

	auto writeUnmappedMateFilteredPair= [&ret,&extractPars,&writer,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		++ret.mateFilteredOffUnmapped_;
		if(extractPars.originalOrientation_){
			writer.openWrite(bamAlnToSeqInfo(bAln));
		}else{
			writer.openWrite(bAlnSeq);
		}
		if(extractPars.writeAll_){
			bWriter.SaveAlignment(bAln);
		}
	};


	auto writeUnmappedMateFilteredSoftClipPair= [&ret,&extractPars,&singleFilteredSoftClipWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		++ret.pairedReadsMateUnmappedFailedSoftClip_;
		if (extractPars.writeAll_) {
			if (extractPars.originalOrientation_) {
				singleFilteredSoftClipWriter.openWrite(bamAlnToSeqInfo(bAln));
			} else {
				singleFilteredSoftClipWriter.openWrite(bAlnSeq);
			}
			bWriter.SaveAlignment(bAln);
		}
	};
	try {
	for (const auto regionPos : iter::range(regions.size())) {
		const auto & region = regions[regionPos];
		if (verbose_) {
			std::cout << region.uid_ << std::endl;
		}
		setBamFileRegionThrow(bReader, region);

		while (bReader.GetNextAlignment(bAln)) {
			//skip secondary alignments
			if (!bAln.IsPrimaryAlignment()) {
				continue;
			}
			if(!extractPars.keepMarkedDuplicate_ && bAln.IsDuplicate()){
				++ret.markedDuplicateFiltered_;
				continue;
			}
			if(extractPars.removeImproperPairs_ && bAln.IsPaired() && !bAln.IsProperPair()){
				bool bothMapped = bAln.IsMapped() && bAln.IsMateMapped();
				//if keeping improper pairs due to one mate not being mapped
				if(!extractPars.keepImproperMateUnmapped_ || bothMapped) {
					++ret.improperPairFiltered_;
					continue;
				}
			}
			//handle non-mapping sequences
			if (bAln.IsPaired() && !bAln.IsMapped() && !bAln.IsMateMapped()) {
				++ret.pairsUnMapped_;
				//only interested in pairs with at least 1 pair mapping
				continue;
			} else if (!bAln.IsPaired() && !bAln.IsMapped()){
				++ret.unpairedUnMapped_;
				//only interested in seqs that are mapping
				continue;
			}
			//get only alignments that fall mostly in this region, setting percInRegion to 0
			//would save any read that had any fall bases in this region
//			if (bAln.IsMapped() && region.getPercInRegion(bAln, refData) < percInRegion) {
//				//if mate is unmapped and it came before this
//				if(!bAln.IsMateMapped() && alnCache.has(bAln.Name)){
//					alnCache.remove(bAln.Name);
//				}
//				continue;
//			}
//			bool print = false;
//			if("M_15515" == bAln.Name){
//				print = true;
//			}
//			if(print){
//				std::cout << njh::json::toJson(bAln) << std::endl;
//			}
			if (bAln.IsPaired()) {
				if (!alnCache.has(bAln.Name)) {
					//enter into cache for until mate is encountered
					alnCache.addWithRegion(bAln, region);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					auto searchRegion = alnCache.getRegion(bAln.Name);
					if (nullptr == search) {
						std::stringstream ss;
						ss << __FILE__ << "  " << __LINE__ << " "<< __PRETTY_FUNCTION__
								<< ", error search shouldn't be able to be nulltpr here"
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					if(nullptr == searchRegion){
						std::stringstream ss;
						ss << __FILE__ << "  " << __LINE__ << " "<< __PRETTY_FUNCTION__
								<< ", error region shouldn't be able to be nulltpr here, something has gone wrong"
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					bool bAlnIn = false;
					bool searchIn = false;

					bool bAlnPassAlnSize = getAlnLen(bAln) >= extractPars.minAlnMapSize_;
					bool searchPassAlnSize  = getAlnLen(*search) >= extractPars.minAlnMapSize_;

					bool bAlnPassSoftClipAmount = getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) < extractPars.softClipPercentageCutOff_;
					bool searchPassSoftClipAmount = getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_;

//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << "bAlnPassSoftClipAmount: " << njh::colorBool(bAlnPassSoftClipAmount) << std::endl;
//					std::cout << "searchPassSoftClipAmount: " << njh::colorBool(searchPassSoftClipAmount) << std::endl;
//					std::cout << region.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//					std::cout << "bAln.QueryBases.size()   : " << bAln.QueryBases.size() << std::endl;
//					std::cout << "search->QueryBases.size(): " << search->QueryBases.size() << std::endl;


					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities, SangerQualOffset);

					seqInfo searchSeq(search->Name, search->QueryBases,search->Qualities, SangerQualOffset);

					if (bAln.IsMapped()) {
						bAlnIn = region.getPercInRegion(bAln, refData) >= extractPars.percInRegion_;
					}
					if(search->IsMapped()){
						searchIn = searchRegion->getPercInRegion(*search, refData) >= extractPars.percInRegion_;
					}
					bool balnAllFilters = bAlnIn && bAlnPassAlnSize;
					bool searchAllFilters = searchIn && searchPassAlnSize;
					if((len(bAlnSeq) > region.getLen() || region.getLen() < 150) && bAlnIn){
						seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
						GenomicRegion balnRegion(bAln, refData);
						uint32_t startRelative = region.start_ - balnRegion.start_;
						uint32_t endRelative = region.end_ - balnRegion.start_;

						seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
						auto alnInfo = bamAlnToAlnInfoLocal(bAln);
						alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
						alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);

						uint32_t startAln = 0;
						if(region.start_ > balnRegion.start_){
							startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
						}
						uint32_t endAln = len(holderSeq);
						if(region.end_ < balnRegion.end_){
							endAln =  getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
						}

						auto outSeq = querySeq.getSubRead(startAln, endAln - startAln);
						outSeq.removeGaps();
						bAlnSeq = outSeq;
					}
					if((len(searchSeq) > searchRegion->getLen() || searchRegion->getLen() < 150) && searchIn){
						seqInfo querySeq = bamAlnToSeqInfo(*search, true);
						GenomicRegion balnRegion(*search, refData);
						uint32_t startRelative = searchRegion->start_ - balnRegion.start_;
						uint32_t endRelative = searchRegion->end_ - balnRegion.start_;
						seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
						auto alnInfo = bamAlnToAlnInfoLocal(*search);
						alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
						alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);
						uint32_t startAln = 0;
						if(region.start_ > balnRegion.start_){
							startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
						}
						uint32_t endAln = len(holderSeq);
						if(region.end_ < balnRegion.end_){
							endAln =  getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
						}
						auto outSeq = querySeq.getSubRead(startAln, endAln - startAln);
						outSeq.removeGaps();
						searchSeq = outSeq;
					}
					if (bAln.IsMapped()) {
						if (region.reverseSrand_) {
							bAlnSeq.reverseComplementRead(false, true);
						}
					}
					if (search->IsMapped()) {
						if (searchRegion->reverseSrand_) {
							searchSeq.reverseComplementRead(false, true);
						}
					}

//					bool balnAllFilters   = bAlnIn   && bAlnPassAlnSize && bAlnPassSoftClipAmount;
//					bool searchAllFilters = searchIn && searchPassAlnSize && searchPassSoftClipAmount;


					if(bAln.IsMapped() && search->IsMapped()){
						if(bAln.RefID == search->RefID && std::abs(bAln.InsertSize) < insertLengthCutOff_){
							//concordant mapping to the current region, re-orient
							if (searchAllFilters && balnAllFilters) {
								if(bAlnPassSoftClipAmount && searchPassSoftClipAmount){
									//check for inverse mapping
									if (bAln.IsReverseStrand() == search->IsReverseStrand()) {
										//writeInversePair(bAln, bAlnSeq, *search, searchSeq);
										if(subSamplingFunction()){
											writeInverseFilteredPair(bAln, *search);
										}
									} else {
										if(subSamplingFunction()){
											writeRegPair(bAln, bAlnSeq, *search, searchSeq);
										}
									}
								} else {
									if(!searchPassSoftClipAmount && !bAlnPassSoftClipAmount){
//										std::cout << __FILE__ << " " << __LINE__ << std::endl;
//										std::cout << "extractPars.softClipPercentageCutOff_: " << extractPars.softClipPercentageCutOff_ << std::endl;
//										std::cout << "getSoftClipAmount(bAln): " << getSoftClipAmount(bAln) << std::endl;
//										std::cout << "static_cast<double>(bAln.QueryBases.size(): " << static_cast<double>(bAln.QueryBases.size()) << std::endl;
//										std::cout << "getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()): " << getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) << std::endl;
//										std::cout << "getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) < extractPars.softClipPercentageCutOff_: " << njh::colorBool(getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) < extractPars.softClipPercentageCutOff_) << std::endl;
//										std::cout << "getSoftClipAmount(search): " << getSoftClipAmount(*search) << std::endl;
//										std::cout << "static_cast<double>(search.QueryBases.size(): " << static_cast<double>(search->QueryBases.size()) << std::endl;
//										std::cout << "getSoftClipAmount(search)/static_cast<double>(search.QueryBases.size()): " << getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) << std::endl;
//										std::cout << "getSoftClipAmount(search)/static_cast<double>(search.QueryBases.size()) < extractPars.softClipPercentageCutOff_: " << njh::colorBool(getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_) << std::endl;
										if(subSamplingFunction()){
											writeBothPairsFilteredSoftClip(bAln, *search);
										}
									}else if(searchPassSoftClipAmount){
										if(subSamplingFunction()){
											writeMateFilteredOffSoftClip(*search, *searchRegion);
											writeTheMateFailedSoftClip(bAln, bAlnSeq);
										}
									}else if(bAlnPassSoftClipAmount){
										if(subSamplingFunction()){
											writeMateFilteredOffSoftClip(bAln, region);
											writeTheMateFailedSoftClip(*search, searchSeq);
										}
									}
								}
							}else if(searchAllFilters){
								if(subSamplingFunction()){
									if(searchPassSoftClipAmount){
										writeMateFilteredOff(*search, *searchRegion);
									}else{
										writeUnmappedMateFilteredSoftClipPair(*search, searchSeq);
									}
									writeTheThrownAwayMate(bAln, bAlnSeq);
								}
							}else if(balnAllFilters){
								if(subSamplingFunction()){
									if(bAlnPassSoftClipAmount){
										writeMateFilteredOff(bAln, region);
									}else{
										writeUnmappedMateFilteredSoftClipPair(bAln, bAlnSeq);
									}
									writeTheThrownAwayMate(*search, searchSeq);
								}
							}else{
								//both searchAllFilters and balnAllFilters are false, which means they both partially map to this region but don't pass the criteria
								//for inclusion, would be good to count how often this happens
								//since we are keeping the thrown away mate why don't we keep both of these to try to map during recruitment
								/**@todo */
								if(subSamplingFunction()){
									writeBothPairsFiltered(bAln, *search);
								}
							}
						} else {
							if (searchAllFilters && balnAllFilters) {
								if(searchPassSoftClipAmount && bAlnPassSoftClipAmount){
									if(subSamplingFunction()){
										//if these checks end up being the same it means the seq is now in the original orientation
										//if false then they are now in the reverse complement of what they use to be
										bool bAlnCheck = region.reverseSrand_ == bAln.IsReverseStrand();
										bool searchCheck = searchRegion->reverseSrand_ == search->IsReverseStrand();
										//if the checks equal each other that means the mates are now in the opposite orientation from each other and therefore are inverse mapping
										if (bAlnCheck == searchCheck) {
											//writeInversePair(bAln, bAlnSeq, *search, searchSeq);
											writeInverseFilteredPair(bAln, *search);
										} else {
											writeDiscordantPair(bAln, bAlnSeq, *search, searchSeq);
										}
									}
								}else{
									if(subSamplingFunction()){
										if(!searchPassSoftClipAmount && !bAlnPassSoftClipAmount){
											writeBothPairsFilteredSoftClip(bAln, *search);
										}else if(searchPassSoftClipAmount){
											writeMateFilteredOffSoftClip(*search, *searchRegion);
											writeTheMateFailedSoftClip(bAln, bAlnSeq);
										}else if(bAlnPassSoftClipAmount){
											writeMateFilteredOffSoftClip(bAln, region);
											writeTheMateFailedSoftClip(*search, searchSeq);
										}
									}
								}
							}else if(searchAllFilters){
								if(subSamplingFunction()){
									if(searchPassSoftClipAmount){
										writeMateFilteredOff(*search, *searchRegion);
									}else{
										writeUnmappedMateFilteredSoftClipPair(*search, searchSeq);
									}
									writeTheThrownAwayMate(bAln, bAlnSeq);
								}
							}else if(balnAllFilters){
								if(subSamplingFunction()){
									if(bAlnPassSoftClipAmount){
										writeMateFilteredOff(bAln, region);
									}else{
										writeUnmappedMateFilteredSoftClipPair(bAln, bAlnSeq);
									}
									writeTheThrownAwayMate(*search, searchSeq);
								}
							}else{
								//both searchAllFilters and balnAllFilters are false, which means they both partially map to this region but don't pass the criteria
								//for inclusion, would be good to count how often this happens
								//since we are keeping the thrown away mate why don't we keep both of these to try to map during recruitment
								/**@todo */
								writeBothPairsFiltered(bAln, *search);
							}
						}
					}else if(bAln.IsMapped()){
						if (extractPars.throwAwayUnmappedMate_) {
							if(balnAllFilters){
								if(subSamplingFunction()){
									if(bAlnPassSoftClipAmount){
										writeThrowAwayUnmappedMate(bAln, bAlnSeq);
									}else{
										writeUnmappedMateFilteredSoftClipPair(bAln, bAlnSeq);
									}
									writeTheThrownAwayUnmappedMate(*search, searchSeq);
								}
							}
						} else {
							//first check to see if un mapped mate is likely falling within the region of interest
							/**@todo should incorporate insert size if possible */
							size_t posibleMatePosition = bAln.Position;
							size_t possibleMatePostionEnd = 0;
							bool reverseStrand = bAln.IsReverseStrand();
							if(bAln.IsReverseStrand()){
								//possibleMatePostionEnd = bAln.Position > static_cast<int64_t>(bAln.QueryBases.size()) ? bAln.Position - bAln.QueryBases.size() : 0;
								possibleMatePostionEnd = bAln.Position;
								posibleMatePosition = possibleMatePostionEnd > search->QueryBases.size() ? possibleMatePostionEnd - search->QueryBases.size() : 0;
							}else{
								posibleMatePosition += bAln.QueryBases.size();
								possibleMatePostionEnd = posibleMatePosition + search->QueryBases.size();
							}
							bool matePass = false;
							if(bAln.IsReverseStrand()){
								//include mate if the read points off the beginning of the region
								if(possibleMatePostionEnd < search->QueryBases.size()){
									matePass = true;
								}
							} else {
								GenomicRegion possibleMateReg(search->Name, refData[bAln.RefID].RefName, posibleMatePosition, possibleMatePostionEnd, reverseStrand);
								double mateBases = search->QueryBases.size();
								matePass = (mateBases > 0 && possibleMateReg.getOverlapLen(region)/mateBases >= extractPars.percInRegion_);
							}

							//if(mateBases > 0 && possibleMateReg.getOverlapLen(region)/mateBases >= percInRegion){
							//if(posibleMatePosition >= region.start_ && posibleMatePosition < region.end_){
							if (matePass) {
								if(!region.reverseSrand_){
									searchSeq.reverseComplementRead(false, true);
								}
								if(subSamplingFunction()){
									if(balnAllFilters){
										if(bAlnPassSoftClipAmount){
											writeMateUnmappedPair(bAln, bAlnSeq, *search, searchSeq);
										}else{
											writeUnmappedMateFilteredSoftClipPair(bAln, bAlnSeq);
											writeUnmappedMateFilteredPair(*search, searchSeq);
										}
									}else{
										writeUnmappedMateFilteredPair(*search, searchSeq);
									}
								}
							}else{
								if(subSamplingFunction()){
									if(balnAllFilters){
										if(bAlnPassSoftClipAmount){
											writeThrowAwayUnmappedMate(bAln, bAlnSeq);
										}else{
											writeUnmappedMateFilteredSoftClipPair(bAln, bAlnSeq);
										}
										writeTheThrownAwayUnmappedMate(*search, searchSeq);
									}
								}
							}
						}
					}else if(search->IsMapped()){
						if (extractPars.throwAwayUnmappedMate_) {
							if(searchAllFilters){
								if(subSamplingFunction()){
									if(searchPassSoftClipAmount){
										writeUnmappedMateFilteredSoftClipPair(*search, searchSeq);
									}else{
										writeThrowAwayUnmappedMate(*search, searchSeq);
									}
									writeTheThrownAwayUnmappedMate(bAln, bAlnSeq);
								}
							}
						} else {
							/**@todo should incorporate insert size if possible */
							//first check to see if un mapped mate is likely falling within the region of interest
							size_t posibleMatePosition = search->Position;
							size_t possibleMatePostionEnd = 0;
							bool reverseStrand = search->IsReverseStrand();
							if(search->IsReverseStrand()){
								//possibleMatePostionEnd = search->Position > static_cast<int64_t>(search->QueryBases.size()) ? search->Position - search->QueryBases.size() : 0;
								possibleMatePostionEnd = search->Position;
								posibleMatePosition = possibleMatePostionEnd > bAln.QueryBases.size() ? possibleMatePostionEnd - bAln.QueryBases.size() : 0;
							} else {
								posibleMatePosition += search->QueryBases.size();
								possibleMatePostionEnd = posibleMatePosition + bAln.QueryBases.size();
							}

							bool matePass = false;
							if(search->IsReverseStrand()){
								//include mate if the read points off the beginning of the region
								if(possibleMatePostionEnd < bAln.QueryBases.size()){
									matePass = true;
								}
							} else {
								GenomicRegion possibleMateReg(bAln.Name, refData[search->RefID].RefName, posibleMatePosition, possibleMatePostionEnd, reverseStrand);
								double mateBases = bAln.QueryBases.size();
								matePass = (mateBases > 0 && possibleMateReg.getOverlapLen(*searchRegion)/mateBases >= extractPars.percInRegion_);
							}
							//if(posibleMatePosition >= searchRegion->start_ && posibleMatePosition < searchRegion->end_){
							if(matePass){
								if(subSamplingFunction()){
									//since to make it here at least one of the mates had to align, just check which one did
									if(!searchRegion->reverseSrand_){
										bAlnSeq.reverseComplementRead(false, true);
									}
									if(searchAllFilters){
										if(searchPassSoftClipAmount){
											writeMateUnmappedPair(bAln, bAlnSeq, *search, searchSeq);
										}else{
											writeUnmappedMateFilteredSoftClipPair(*search, searchSeq);
											writeUnmappedMateFilteredPair(bAln, bAlnSeq);
										}
									}else{
										writeUnmappedMateFilteredPair(bAln, bAlnSeq);
									}
								}
							} else {
								if(subSamplingFunction()){
									if (searchAllFilters) {
										if(searchPassSoftClipAmount){
											writeThrowAwayUnmappedMate(*search, searchSeq);
										}else{
											writeUnmappedMateFilteredSoftClipPair(*search, searchSeq);
										}
										writeTheThrownAwayUnmappedMate(bAln, bAlnSeq);
									}
								}
							}
						}
					}else{
						std::stringstream ss;
						ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__
								<< ", error shouldn't be able to reach here"
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {

				if(region.getPercInRegion(bAln, refData) >= extractPars.percInRegion_ && getAlnLen(bAln) >= extractPars.minAlnMapSize_){
					if(getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) > extractPars.softClipPercentageCutOff_){
						if(subSamplingFunction()){
							writeSingleFilteredSoftClip(bAln);
						}
					}else{
						if(subSamplingFunction()){
							//unpaired read
							++ret.unpaiedReads_;
							if (extractPars.originalOrientation_) {
								writer.openWrite(bamAlnToSeqInfo(bAln));
							} else {
								seqInfo outSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
										SangerQualOffset);

								if((len(outSeq) > region.getLen() || region.getLen() < 150) ){
									seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
									GenomicRegion balnRegion(bAln, refData);
									uint32_t startRelative = region.start_ - balnRegion.start_;
									uint32_t endRelative = region.end_ - balnRegion.start_;
									seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
									auto alnInfo = bamAlnToAlnInfoLocal(bAln);
									alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
									alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);
									uint32_t startAln = 0;
									if(region.start_ > balnRegion.start_){
										startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
									}
									uint32_t endAln = len(holderSeq);
									if(region.end_ < balnRegion.end_){
										endAln = getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
									}
									auto outTrimmedSeq = querySeq.getSubRead(startAln, endAln - startAln);
									outTrimmedSeq.removeGaps();
									outSeq = outTrimmedSeq;
								}
								//put in the orientation of the output region
								if(region.reverseSrand_){
									outSeq.reverseComplementRead(false, true);
								}
								writer.openWrite(outSeq);
							}
						}
					}
				}else{
					//single filtered off
					if(subSamplingFunction()){
						writeSingleFiltered(bAln);
					}
				}
			}
		}
	}

	//save the orphans;
	/**@todo these will mostly be pairs that had mates that didn't fall in any regions,
	 *  should maybe investigate where they are falling or when remapping consider taking these now  */
	std::ofstream orphansSiblingsLocationFile;
	if(debug_){
		OutOptions orphansSiblingsLocationFileOpt(njh::files::make_path(outOpts.outFilename_.parent_path(), "orphansSiblingsStats.tab.txt"));
		orphansSiblingsLocationFileOpt.openFile(orphansSiblingsLocationFile);
		orphansSiblingsLocationFile << "OrphanName\tdistance\tclosetRegions\tRefId\tPosition\tIsMapped\tMateRefId\tMatePosition\tIsMateMapped" << "\n";
	}

	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		if (extractPars.tryToFindOrphansMate_) {
			//find orphans' mates if possible;
//			BamTools::BamReader bReaderMateFinder;
//			bReaderMateFinder.Open(inOutOpts.firstName_.string());
//			checkBamOpenThrow(bReaderMateFinder, inOutOpts.firstName_.string());
//			loadBamIndexThrow(bReaderMateFinder);
//			auto refData = bReaderMateFinder.GetReferenceData();
			//gather all the orphans regions
			std::unordered_map<std::string, std::set<uint32_t>> orphanPositions;
			for (const auto & name : names) {
				auto search = alnCache.get(name);
				if(search->IsPaired()){
					if (search->IsMateMapped()) {
						bool pass = true;
						if (extractPars.filterOffLowEntropyOrphansRecruits_) {
							seqInfo searchSeq(search->Name, search->QueryBases, search->Qualities, SangerQualOffset);
							kmerInfo kInfo(searchSeq.seq_, extractPars.entropyKlen_, false);
							if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
								pass = false;
							}
						}
						if(pass){
							orphanPositions[refData[search->MateRefID].RefName].emplace(search->MatePosition);
						}
					}
				}
			}
			std::vector<GenomicRegion> orphanMateRegions;
			for(const auto & orPos : orphanPositions){
				for(const auto & pos  : orPos.second){
					orphanMateRegions.emplace_back(GenomicRegion("", orPos.first, pos, pos + 1, false));
				}
			}
			for(const auto & reg : orphanMateRegions){
				setBamFileRegionThrow(bReader, reg);
				while (bReader.GetNextAlignment(bAln)) {
					//skip secondary alignments
					if (!bAln.IsPrimaryAlignment()) {
						continue;
					}
					if(static_cast<uint32_t>(bAln.Position) < reg.start_){
						continue;
					}
					if (bAln.IsPaired()) {
						if (alnCache.has(bAln.Name)) {
							auto search = alnCache.get(bAln.Name);
							if (bAln.IsFirstMate() != search->IsFirstMate()) {
								bool pass = true;
								auto balnSeq = seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities, SangerQualOffset);
								if (extractPars.filterOffLowEntropyOrphansRecruits_) {
									kmerInfo kInfo(balnSeq.seq_, extractPars.entropyKlen_, false);
									if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
										pass = false;
									}
								}
								if(pass){
									if(subSamplingFunction()){

										writeTheThrownAwayMate(bAln, balnSeq);
									}
								}
							}
						}
					}
				}
			}
		} else {
			for(const auto & region : regions){
//				BamTools::BamReader bReaderMateFinder;
//				bReaderMateFinder.Open(inOutOpts.firstName_.string());
//				checkBamOpenThrow(bReaderMateFinder, inOutOpts.firstName_.string());
//				loadBamIndexThrow(bReaderMateFinder);
				setBamFileRegionThrow(bReader, region);
				while (bReader.GetNextAlignment(bAln)) {
					//skip secondary alignments
					if (!bAln.IsPrimaryAlignment()) {
						continue;
					}
					if (bAln.IsPaired()) {
						if (alnCache.has(bAln.Name)) {
							auto search = alnCache.get(bAln.Name);
							if (bAln.IsFirstMate() != search->IsFirstMate()) {
								/**@todo look into this...*/
								if(subSamplingFunction()){
									writeTheThrownAwayMate(bAln,
											seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
													SangerQualOffset));
								}
							}
						}
					}
				}
			}
		}


		for (const auto & name : names) {
			auto search = alnCache.get(name);
			if(debug_){
				if(search->IsPaired()){
					if (search->IsMateMapped()) {
						int32_t smallestDistance = std::numeric_limits<int32_t>::max();
						std::vector<std::string> bestRegions;
						for (const auto & reg : regions) {
							if (reg.chrom_ != refData.at(search->MateRefID).RefName) {
								continue;
							}
							int32_t frontDiff = search->MatePosition
									- static_cast<int64_t>(reg.start_);
							if (std::abs(frontDiff) < std::abs(smallestDistance)) {
								smallestDistance = frontDiff;
								bestRegions.clear();
								bestRegions.emplace_back(reg.uid_);
							} else if (std::abs(frontDiff) == std::abs(smallestDistance)) {
								bestRegions.emplace_back(reg.uid_);
							}
						}
						orphansSiblingsLocationFile << search->Name
								<< "\t" << (bestRegions.empty() ? "none" : estd::to_string(smallestDistance))
								<< "\t" << (bestRegions.empty() ? "none" : njh::conToStr(bestRegions, ","))
								<< "\t" << (search->RefID < 0 ? std::string("*") : refData[search->RefID].RefName)
								<< "\t" << search->Position
								<< "\t" << njh::boolToStr(search->IsMapped())
								<< "\t" << (search->MateRefID < 0 ? std::string("*") : refData[search->MateRefID].RefName)
								<< "\t" << search->MatePosition
								<< "\t" << njh::boolToStr(search->IsMateMapped())
								<< "\n";
					} else {
						orphansSiblingsLocationFile << search->Name
								<< "\t" << "*"
								<< "\t" << "*"
								<< "\t" << (search->RefID < 0 ? std::string("*") : refData[search->RefID].RefName)
								<< "\t" << search->Position
								<< "\t" << njh::boolToStr(search->IsMapped())
								<< "\t" << (search->MateRefID < 0 ? std::string("*") : refData[search->MateRefID].RefName)
								<< "\t" << search->MatePosition
								<< "\t" << njh::boolToStr(search->IsMateMapped())
								<< "\n";
					}
				}
			}
			if(!search->IsMapped()){
				if(search->IsPaired()){
					++ret.orphansUnmapped_;
				}else{
					++ret.unpairedUnMapped_;
				}
				alnCache.remove(name);
				continue;
			}

			auto searchRegion = alnCache.getRegion(name);
			bool searchIn = searchRegion->getPercInRegion(*search, refData)>= extractPars.percInRegion_;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << search->Name << std::endl;
//			std::cout << "\tgetAlnLen(*search): " << getAlnLen(*search) << std::endl;
//			std::cout << "\tnjh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_): " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//			std::cout << "\tsearchIn: " << njh::colorBool(searchIn) << std::endl;
//			std::cout << "\tsearchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
			if (searchIn &&
					getAlnLen(*search) >= extractPars.minAlnMapSize_) {
//				std::cout << "\tgetSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_: " << njh::colorBool(getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_) << std::endl;
				if(getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_){
					if (search->IsPaired()) {
						++ret.orphans_;
					} else {
						++ret.unpaiedReads_;
					}
					bool pass = true;
					auto searchSeq = seqInfo(search->Name, search->QueryBases, search->Qualities, SangerQualOffset);
//					std::cout << searchSeq.name_ << std::endl;

					if (extractPars.filterOffLowEntropyOrphansRecruits_) {
						kmerInfo kInfo(searchSeq.seq_, extractPars.entropyKlen_, false);
//						std::cout << "\t" << kInfo.computeKmerEntropy() << std::endl;
//						std::cout << "\t" << extractPars.filterOffLowEntropyOrphansRecruitsCutOff_ << std::endl;
//						std::cout << "\t" << njh::colorBool(kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) << std::endl;
						if (kInfo.computeKmerEntropy() < extractPars.filterOffLowEntropyOrphansRecruitsCutOff_) {
							pass = false;
						}
					}

					if(pass){
						if (extractPars.originalOrientation_) {
							if(subSamplingFunction()){
								if((len(searchSeq) > searchRegion->getLen() || searchRegion->getLen() < 150) && searchIn){
									seqInfo querySeq = bamAlnToSeqInfo(*search, true);
									GenomicRegion balnRegion(*search, refData);
									uint32_t startRelative = searchRegion->start_ - balnRegion.start_;
									uint32_t endRelative = searchRegion->end_ - balnRegion.start_;
									seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
									auto alnInfo = bamAlnToAlnInfoLocal(*search);
									alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
									alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);
									uint32_t startAln = 0;
									if(searchRegion->start_ > balnRegion.start_){
										startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
									}
									uint32_t endAln = len(holderSeq);
									if(searchRegion->end_ < balnRegion.end_){
										endAln =  getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
									}
									auto outSeq = querySeq.getSubRead(startAln, endAln - startAln);
									outSeq.removeGaps();
									searchSeq = outSeq;
								}
								if(search->IsReverseStrand()){
									searchSeq.reverseComplementRead(false, true);
								}
								writer.openWrite(searchSeq);
								//writer.openWrite(bamAlnToSeqInfo(*search));
							}
						} else {
							seqInfo searchSeq(search->Name, search->QueryBases, search->Qualities, SangerQualOffset);
							if((len(searchSeq) > searchRegion->getLen() || searchRegion->getLen() < 150) && searchIn){
								seqInfo querySeq = bamAlnToSeqInfo(*search, true);
								GenomicRegion balnRegion(*search, refData);
								uint32_t startRelative = searchRegion->start_ - balnRegion.start_;
								uint32_t endRelative = searchRegion->end_ - balnRegion.start_;
								seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
								auto alnInfo = bamAlnToAlnInfoLocal(*search);
								alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
								alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);
								uint32_t startAln = 0;
								if(searchRegion->start_ > balnRegion.start_){
									startAln = getAlnPosForRealPos(holderSeq.seq_, startRelative);
								}
								uint32_t endAln = len(holderSeq);
								if(searchRegion->end_ < balnRegion.end_){
									endAln = getAlnPosForRealPos(holderSeq.seq_, endRelative - 1) + 1;
								}
								auto outSeq = querySeq.getSubRead(startAln, endAln - startAln);
								outSeq.removeGaps();
								searchSeq = outSeq;
							}
							if (searchRegion->reverseSrand_) {
								searchSeq.reverseComplementRead(false, true);
							}
							if(subSamplingFunction()){
								writer.openWrite(searchSeq);
							}
						}
						if(extractPars.writeAll_){
							bWriter.SaveAlignment(*search);
						}
					}else{
						++ret.orphansFiltered_;
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						std::cout << search->Name << std::endl;
//						std::cout << "\tgetAlnLen(*search): " << getAlnLen(*search) << std::endl;
//						std::cout << "\tnjh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_): " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//						std::cout << "\tsearchIn: " << njh::colorBool(searchIn) << std::endl;
//						std::cout << "\tsearchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						std::cout << "getAlnLen(*search) >= extractPars.minAlnMapSize_: " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//						std::cout << "getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_: " << njh::colorBool(getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_) << std::endl;
//						std::cout << "pass: " << njh::colorBool(pass) << std::endl;
//						std::cout <<  "searchIn: " <<njh::colorBool(searchIn) << std::endl;
						if(extractPars.writeAll_){
							if(subSamplingFunction()){
								singleFilteredWriter.openWrite(bamAlnToSeqInfo(*search));
							}
							bWriter.SaveAlignment(*search);
						}
					}
				}else{
					//write filterd orphan
					++ret.orphansFilteredSoftCip_;
					if(extractPars.writeAll_){
						if(subSamplingFunction()){
							singleFilteredSoftClipWriter.openWrite(bamAlnToSeqInfo(*search));
						}
						bWriter.SaveAlignment(*search);
					}
				}
			}else{
				//write filterd orphan
				++ret.orphansFiltered_;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << search->Name << std::endl;
//				std::cout << "\tgetAlnLen(*search): " << getAlnLen(*search) << std::endl;
//				std::cout << "\tnjh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_): " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//				std::cout << "\tsearchIn: " << njh::colorBool(searchIn) << std::endl;
//				std::cout << "\tsearchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << "getAlnLen(*search) >= extractPars.minAlnMapSize_: " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//				std::cout << "getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_: " << njh::colorBool(getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_) << std::endl;
//				std::cout <<  "searchIn: " <<njh::colorBool(searchIn)  << std::endl;
				if(extractPars.writeAll_){
					if(subSamplingFunction()){
						singleFilteredWriter.openWrite(bamAlnToSeqInfo(*search));
					}
					bWriter.SaveAlignment(*search);
				}
			}
//			std::cout << std::endl << std::endl;
			alnCache.remove(name);
		}
	}
	} catch (std::exception & e) {
		std::cout << e.what()	<< std::endl;
		exit(1);
	}
	return ret;
}




BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsWtihCrossRegionMapping(
		const SeqIOOptions & inOutOpts,
		const std::vector<GenomicRegion> & regions,
		const extractReadsWtihCrossRegionMappingPars & extractPars) {
	BamTools::BamReader bReader;
	bReader.Open(inOutOpts.firstName_.string());
	checkBamOpenThrow(bReader, inOutOpts.firstName_.string());
	loadBamIndexThrow(bReader);

	return extractReadsWtihCrossRegionMapping(bReader, inOutOpts.out_, regions, extractPars);
}

}  // namespace njhseq
