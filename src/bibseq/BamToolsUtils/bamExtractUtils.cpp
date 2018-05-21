/*
 * bamExtractUtils.cpp
 *
 *  Created on: Feb 2, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bamExtractUtils.hpp"
#include "bibseq/BamToolsUtils/BamAlnsCache.hpp"
#include "bibseq/BamToolsUtils/BamAlnsCacheWithRegion.hpp"
#include "bibseq/IO/SeqIO.h"

namespace bibseq {


uint32_t BamExtractor::ExtractCounts::getTotal(){
	return pairedReads_ + pairedReadsMateUnmapped_ + unpaiedReads_ + orphans_ + orphansUnmapped_ + pairFilteredOff_ + pairFilteredOffUnmapped_ + pairsUnMapped_ + unpairedUnMapped_ + discordant_ + inverse_;;
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
	logStats("singles", unpaiedReads_);
	logStats("orphans", orphans_);
	logStats("orphansUnmapped", orphansUnmapped_);
	logStats("pairFilteredOff", pairFilteredOff_);
	logStats("pairFilteredOffUnmapped", pairFilteredOffUnmapped_);
	logStats("unmappedPaired", pairsUnMapped_);
	logStats("unmappedSingles", unpairedUnMapped_);
	logStats("discordant", discordant_);
	logStats("inverse", inverse_);
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

	removeIfExists(inPairsUnMapped_.firstName_);
	removeIfExists(inPairsUnMapped_.secondName_);

	removeIfExists(inUnpairedUnMapped_.firstName_);

	removeIfExists(inInverse_.firstName_);
	removeIfExists(inInverse_.secondName_);

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
	bib::sys::requireExternalProgramThrow("flash");

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
		auto stitchFilesBase = bib::files::prependFileBasename(outOpts.outName(),
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
				<< bib::files::make_path(stitchFilesBase.parent_path(), "flashLog")
				<< " 2>&1\n";
		auto flashOut = bib::sys::run(VecStr { flashCmd.str() });
		if (!flashOut.success_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error in running flash to stitch reads"
					<< "\n";
			ss << "Mess: " << "\n";
			ss << flashOut.stdErr_ << "\n";
			throw std::runtime_error { ss.str() };
		}
		OutOptions flashRunLogOpts(
				bib::files::make_path(stitchFilesBase.parent_path(),
						"flashRunLog.json"));
		std::ofstream outFlashRunLog;
		flashRunLogOpts.openFile(outFlashRunLog);
		outFlashRunLog << flashOut.toJson() << std::endl;
		ret.stitched_ = SeqIOOptions::genFastqIn(
				bib::appendAsNeededRet(stitchFilesBase.string(),
						".extendedFrags.fastq"));

		ret.notCombinedPairs_ = SeqIOOptions::genPairedIn(
				bib::appendAsNeededRet(stitchFilesBase.string(),
						".notCombined_1.fastq"),
				bib::appendAsNeededRet(stitchFilesBase.string(),
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
			writer.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
}


void BamExtractor::writeExtractReadsFromBamOnlyMapped(const bfs::path & bamFnp,
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
			writer.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
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
								possibleMatePostionEnd = bAln.Position > bAln.QueryBases.size() ? bAln.Position - bAln.QueryBases.size() : 0;
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
								possibleMatePostionEnd = search->Position > search->QueryBases.size() ? search->Position - search->QueryBases.size() : 0;
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
			bib::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
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
BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsFromBamWrite(
		const SeqIOOptions & opts, bool referenceOrientation, bool throwAwayUnmmpaedMates) {

	auto outPairsUnMapped = SeqIOOptions::genPairedOut(
			bib::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outPairsUnMapped.out_.transferOverwriteOpts(opts.out_);
	auto outUnpairedUnMapped = SeqIOOptions::genFastqOut(
			bib::files::prependFileBasename(opts.out_.outFilename_, "unmapped_"));
	outUnpairedUnMapped.out_.transferOverwriteOpts(opts.out_);

	auto outPairs = SeqIOOptions::genPairedOut(opts.out_.outFilename_);
	outPairs.out_.transferOverwriteOpts(opts.out_);
	auto outUnpaired = SeqIOOptions::genFastqOut(opts.out_.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(opts.out_);

	auto outPairsMateUnmapped = SeqIOOptions::genPairedOut(
			bib::files::prependFileBasename(opts.out_.outFilename_, "mateUnmapped_"));
	outPairsMateUnmapped.out_.transferOverwriteOpts(opts.out_);

	auto thrownAwayMate = SeqIOOptions::genFastqOut(
			bib::files::prependFileBasename(opts.out_.outFilename_, "thrownAwayMate_"));
	thrownAwayMate.out_.transferOverwriteOpts(opts.out_);


	auto outPairsInverse = SeqIOOptions::genPairedOut(
			bib::files::prependFileBasename(opts.out_.outFilename_, "inverse_"));
	outPairsInverse.out_.transferOverwriteOpts(opts.out_);

	auto outPairsDiscordant = SeqIOOptions::genPairedOut(
			bib::files::prependFileBasename(opts.out_.outFilename_, "discordant_"));
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
									&& bAln.IsReverseStrand() != search->IsReverseStrand()
									&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
								++ret.pairedReads_;
								if (bAln.IsFirstMate()) {
									mappedPairWriter.openWrite(
											PairedRead(bAlnSeq, searchSeq, false));
								} else {
									mappedPairWriter.openWrite(
											PairedRead(searchSeq, bAlnSeq, false));
								}
							} else {
								//inverse if closely mapping but to same strand
								if (bAln.IsReverseStrand() != search->IsReverseStrand()
										&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
									//inverse mates will there be written in technically the wrong orientation to each other but in the reference orientation
									if (bAln.IsFirstMate()) {
										inversePairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										inversePairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
									++ret.inverse_;
								}else{
									//discordant if mapping to different chromosome or very far away from each other
									if (bAln.IsFirstMate()) {
										discordantPairWriter.openWrite(
												PairedRead(bAlnSeq, searchSeq, false));
									} else {
										discordantPairWriter.openWrite(
												PairedRead(searchSeq, bAlnSeq, false));
									}
									++ret.discordant_;
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
								&& bAln.IsReverseStrand() != search->IsReverseStrand()
								&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
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
						} else {
							if (bAln.IsReverseStrand() != search->IsReverseStrand()
									&& std::abs(bAln.InsertSize) < insertLengthCutOff_) {
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

	BamTools::BamAlignment bAln;
	BamAlnsCache alnCache;

	bfs::path outBam(bib::appendAsNeededRet(opts.out_.outFilename_.string(), ".bam"));
	OutOptions outBamOpts(outBam);
	outBamOpts.transferOverwriteOpts(opts.out_);
	if(outBamOpts.outExists() && !outBamOpts.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outBamOpts.outName() << " already exits" << "\n";
		throw std::runtime_error{ss.str()};
	}

	BamTools::BamWriter bWriter;
	bWriter.Open(outBam.string(), bReader.GetHeader(), bReader.GetReferenceData());


	while(bReader.GetNextAlignment(bAln)){
		if(!bAln.IsPaired() && !bAln.IsMapped()){
			++ret.unpairedUnMapped_;
			singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
			bWriter.SaveAlignment(bAln);
		}else{
			if(!bAln.IsMapped() && !bAln.IsMateMapped()){
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
	for(const auto & regPos : iter::range(regions.size())){
		for(const auto & secondPos : iter::range<size_t>(0, regPos)){
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

BamExtractor::ExtractedFilesOpts BamExtractor::extractReadsWtihCrossRegionMapping(
		const SeqIOOptions & inOutOpts,
		const std::vector<GenomicRegion> & regions,
		double percInRegion,
		bool originalOrientation,
		bool throwAwayUnmappedMate,
		bool tryToFindOrphansMate) {
	//std::cout << "percInRegion: " << percInRegion << std::endl;
	//check to see if regions overlap
	bool overlapsFound = false;
	std::stringstream overalMessage;
	overalMessage << __PRETTY_FUNCTION__ << " error, found overlaps, overlap regions should be merged" << "\n";
	for(const auto & regPos : iter::range(regions.size())){
		for(const auto & secondPos : iter::range<size_t>(0, regPos)){
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
	BamTools::BamAlignment bAln;
	BamTools::BamReader bReader;
	bReader.Open(inOutOpts.firstName_.string());
	checkBamOpenThrow(bReader, inOutOpts.firstName_.string());
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCacheWithRegion alnCache;

	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	bfs::path outBam(bib::appendAsNeededRet(inOutOpts.out_.outFilename_.string(), ".bam"));
	OutOptions outBamOpts(outBam);
	outBamOpts.transferOverwriteOpts(inOutOpts.out_);
	if(outBamOpts.outExists() && !outBamOpts.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outBamOpts.outName() << " already exits" << "\n";
		throw std::runtime_error{ss.str()};
	}

	BamTools::BamWriter bWriter;
	bWriter.Open(outBam.string(), bReader.GetHeader(), bReader.GetReferenceData());


	auto outPairs = SeqIOOptions::genPairedOut(inOutOpts.out_.outFilename_);
	outPairs.out_.transferOverwriteOpts(inOutOpts.out_);

	auto outPairsUnmappedMate = SeqIOOptions::genPairedOut(bib::files::prependFileBasename(inOutOpts.out_.outFilename_, "mateUnmapped_"));
	outPairsUnmappedMate.out_.transferOverwriteOpts(inOutOpts.out_);

	auto thrownAwayUnammpedMateOpts = SeqIOOptions::genFastqOut(bib::files::prependFileBasename(inOutOpts.out_.outFilename_, "thrownAwayMate_"));
	thrownAwayUnammpedMateOpts.out_.transferOverwriteOpts(inOutOpts.out_);

	//inverse here is if when the pairs are put into the orientation of the target region of interest they end up being not being in the same orientation or regular inverse as well
	auto outPairsInverse = SeqIOOptions::genPairedOut(bib::files::prependFileBasename(inOutOpts.out_.outFilename_, "inverse_"));
	outPairsInverse.out_.transferOverwriteOpts(inOutOpts.out_);

	auto outUnpaired = SeqIOOptions::genFastqOut(inOutOpts.out_.outFilename_);
	outUnpaired.out_.transferOverwriteOpts(inOutOpts.out_);
	//pair writers
	SeqOutput pairWriter(outPairs);
	SeqOutput mateUnmappedPairWriter(outPairsUnmappedMate);

	SeqOutput thrownAwayUnammpedMateWriter(thrownAwayUnammpedMateOpts);
	SeqOutput inversePairWriter(outPairsInverse);
	//non paired writer
	SeqOutput writer(outUnpaired);

	BamExtractor::ExtractedFilesOpts ret;

	ret.inPairs_ =                  SeqIOOptions::genPairedIn(outPairs.getPriamryOutName(), outPairs.getSecondaryOutName());
	ret.inPairsMateUnmapped_ =      SeqIOOptions::genPairedIn(outPairsUnmappedMate.getPriamryOutName(), outPairsUnmappedMate.getSecondaryOutName());
	ret.inThrownAwayUnmappedMate_ = SeqIOOptions::genFastqIn (thrownAwayUnammpedMateOpts.getPriamryOutName());
	ret.inInverse_ =                SeqIOOptions::genPairedIn(outPairsInverse.getPriamryOutName(), outPairsInverse.getSecondaryOutName());
	//no disconcordant reads as this is aiming to grab those reads
	ret.inUnpaired_ =               SeqIOOptions::genFastqIn (outUnpaired.getPriamryOutName());

	auto writeMateFilteredOff = [&ret,&originalOrientation,&writer,&bWriter](const BamTools::BamAlignment & bAln, const GenomicRegion & region){
		//unpaired read
		++ret.pairFilteredOff_;
		if (originalOrientation) {
			writer.openWrite(bamAlnToSeqInfo(bAln));
		} else {
			seqInfo outSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
					SangerQualOffset);
			//put in the orientation of the output region
			if(region.reverseSrand_){
				outSeq.reverseComplementRead(false, true);
			}
			writer.openWrite(outSeq);
		}
		bWriter.SaveAlignment(bAln);
	};

	auto writeInversePair = [&ret,&originalOrientation,&inversePairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.inverse_;
		if(originalOrientation){
			if (bAln.IsFirstMate()) {
				inversePairWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(searchAln)));
			} else {
				inversePairWriter.openWrite(PairedRead(bamAlnToSeqInfo(searchAln), bamAlnToSeqInfo(bAln)));
			}
		}else{
			if (bAln.IsFirstMate()) {
				inversePairWriter.openWrite(PairedRead(bAlnSeq, searchSeq));
			} else {
				inversePairWriter.openWrite(PairedRead(searchSeq, bAlnSeq));
			}
		}
		bWriter.SaveAlignment(bAln);
		bWriter.SaveAlignment(searchAln);
	};
	auto writeDiscordantPair = [&ret,&originalOrientation,&pairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.discordant_;
		if(originalOrientation){
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
		bWriter.SaveAlignment(bAln);
		bWriter.SaveAlignment(searchAln);
	};
	auto writeRegPair = [&ret,&originalOrientation,&pairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.pairedReads_;
		if(originalOrientation){
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
		bWriter.SaveAlignment(bAln);
		bWriter.SaveAlignment(searchAln);
	};
	auto writeMateUnmappedPair = [&ret,&originalOrientation,&mateUnmappedPairWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq,
			const BamTools::BamAlignment & searchAln, const seqInfo & searchSeq){
		++ret.pairedReadsMateUnmapped_;
		if(originalOrientation){
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
		bWriter.SaveAlignment(bAln);
		bWriter.SaveAlignment(searchAln);
	};
	auto writeThrowAwayUnmappedMate = [&ret,&originalOrientation,&writer,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		++ret.pairedReadsMateUnmapped_;
		if(originalOrientation){
			writer.openWrite(bamAlnToSeqInfo(bAln));
		}else{
			writer.openWrite(bAlnSeq);
		}
		bWriter.SaveAlignment(bAln);
	};

	auto writeTheThrownAwayUnmappedMate = [&originalOrientation,&thrownAwayUnammpedMateWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		if(originalOrientation){
			thrownAwayUnammpedMateWriter.openWrite(bamAlnToSeqInfo(bAln));
		}else{
			thrownAwayUnammpedMateWriter.openWrite(bAlnSeq);
		}
		bWriter.SaveAlignment(bAln);
	};

	auto writeTheThrownAwayMate = [&originalOrientation,&thrownAwayUnammpedMateWriter,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		if(originalOrientation){
			thrownAwayUnammpedMateWriter.openWrite(bamAlnToSeqInfo(bAln));
		}else{
			thrownAwayUnammpedMateWriter.openWrite(bAlnSeq);
		}
		bWriter.SaveAlignment(bAln);
	};

	auto writeUnmappedMateFilteredPair= [&ret,&originalOrientation,&writer,&bWriter](const BamTools::BamAlignment & bAln, const seqInfo & bAlnSeq){
		++ret.pairFilteredOffUnmapped_;
		if(originalOrientation){
			writer.openWrite(bamAlnToSeqInfo(bAln));
		}else{
			writer.openWrite(bAlnSeq);
		}
		bWriter.SaveAlignment(bAln);
	};

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
//				std::cout << bib::json::toJson(bAln) << std::endl;
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



					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
							SangerQualOffset);
					seqInfo searchSeq(search->Name, search->QueryBases,
							search->Qualities, SangerQualOffset);

					if (bAln.IsMapped()) {
						bAlnIn = region.getPercInRegion(bAln, refData) >= percInRegion;
						if (region.reverseSrand_) {
							bAlnSeq.reverseComplementRead(false, true);
						}
					}

					if(search->IsMapped()){
						searchIn = searchRegion->getPercInRegion(*search, refData) >= percInRegion;
						if(searchRegion->reverseSrand_){
							searchSeq.reverseComplementRead(false, true);
						}
					}

					if(bAln.IsMapped() && search->IsMapped()){
						if(bAln.RefID == search->RefID && std::abs(bAln.InsertSize) < insertLengthCutOff_){
							//concordant mapping to the current region, reorient
							if (searchIn && bAlnIn) {
								//check for inverse mapping
								if (bAln.IsReverseStrand() == search->IsReverseStrand()) {
									writeInversePair(bAln, bAlnSeq, *search, searchSeq);
								} else {
									writeRegPair(bAln, bAlnSeq, *search, searchSeq);
								}
							}else if(searchIn){
								writeMateFilteredOff(*search, *searchRegion);
								writeTheThrownAwayMate(bAln, bAlnSeq);
							}else if(bAlnIn){
								writeMateFilteredOff(bAln, region);
								writeTheThrownAwayMate(*search, searchSeq);
							}
						}else{
							if (searchIn && bAlnIn) {
								//if these checks end up being the same it means the seq is now in the original orientation
								//if false then they are now in the reverse complement of what they use to be
								bool bAlnCheck = region.reverseSrand_ == bAln.IsReverseStrand();
								bool searchCheck = searchRegion->reverseSrand_ == search->IsReverseStrand();
								//if the checks equal each other that means the mates are now in the opposite orientation from each other and therefore are inverse mapping
								if(bAlnCheck == searchCheck){
									writeInversePair(bAln, bAlnSeq, *search, searchSeq);
								}else{
									writeDiscordantPair(bAln, bAlnSeq, *search, searchSeq);
								}
							}else if(searchIn){
								writeMateFilteredOff(*search, *searchRegion);
								writeTheThrownAwayMate(bAln, bAlnSeq);
							}else if(bAlnIn){
								writeMateFilteredOff(bAln, region);
								writeTheThrownAwayMate(*search, searchSeq);
							}
						}
					}else if(bAln.IsMapped()){
						if (throwAwayUnmappedMate) {
							if(bAlnIn){
								writeThrowAwayUnmappedMate(bAln, bAlnSeq);
								writeTheThrownAwayUnmappedMate(*search, searchSeq);
							}
						} else {
							//first check to see if un mapped mate is likely falling within the region of interest
							/**@todo should incorporate insert size if possible */
							size_t posibleMatePosition = bAln.Position;
							size_t possibleMatePostionEnd = 0;
							bool reverseStrand = bAln.IsReverseStrand();
							if(bAln.IsReverseStrand()){
								possibleMatePostionEnd = bAln.Position > bAln.QueryBases.size() ? bAln.Position - bAln.QueryBases.size() : 0;
								posibleMatePosition = possibleMatePostionEnd > search->QueryBases.size() ? possibleMatePostionEnd - search->QueryBases.size() : 0;
							}else{
								posibleMatePosition += bAln.QueryBases.size();
								possibleMatePostionEnd = posibleMatePosition + search->QueryBases.size();
							}
							GenomicRegion possibleMateReg(search->Name, refData[bAln.RefID].RefName, posibleMatePosition, possibleMatePostionEnd, reverseStrand);
							double mateBases = search->QueryBases.size();
							if(mateBases > 0 && possibleMateReg.getOverlapLen(region)/mateBases >= percInRegion){
							//if(posibleMatePosition >= region.start_ && posibleMatePosition < region.end_){

								//since to make it here at least one of the mates had to align, just check which one did
								if(!region.reverseSrand_){
									searchSeq.reverseComplementRead(false, true);
								}
								if(bAlnIn){
									writeMateUnmappedPair(bAln, bAlnSeq, *search, searchSeq);
								}else{
									writeUnmappedMateFilteredPair(*search, searchSeq);
								}
							}else{
								if(bAlnIn){
									writeThrowAwayUnmappedMate(bAln, bAlnSeq);
									writeTheThrownAwayUnmappedMate(*search, searchSeq);
								}
							}
						}
					}else if(search->IsMapped()){
						if (throwAwayUnmappedMate) {
							if(searchIn){
								writeThrowAwayUnmappedMate(*search, searchSeq);
								writeTheThrownAwayUnmappedMate(bAln, bAlnSeq);
							}
						} else {
							/**@todo should incorporate insert size if possible */
							//first check to see if un mapped mate is likely falling within the region of interest
							size_t posibleMatePosition = search->Position;
							size_t possibleMatePostionEnd = 0;
							bool reverseStrand = search->IsReverseStrand();
							if(search->IsReverseStrand()){
								possibleMatePostionEnd = search->Position > search->QueryBases.size() ? search->Position - search->QueryBases.size() : 0;
								posibleMatePosition = possibleMatePostionEnd > bAln.QueryBases.size() ? possibleMatePostionEnd - bAln.QueryBases.size() : 0;
							}else{
								posibleMatePosition += search->QueryBases.size();
								possibleMatePostionEnd = posibleMatePosition + bAln.QueryBases.size();
							}
							GenomicRegion possibleMateReg(bAln.Name, refData[search->RefID].RefName,posibleMatePosition, possibleMatePostionEnd, reverseStrand);
							double mateBases = bAln.QueryBases.size();

							//if(posibleMatePosition >= searchRegion->start_ && posibleMatePosition < searchRegion->end_){
							if(mateBases > 0 && possibleMateReg.getOverlapLen(*searchRegion)/mateBases >= percInRegion){
								//since to make it here at least one of the mates had to align, just check which one did
								if(!searchRegion->reverseSrand_){
									bAlnSeq.reverseComplementRead(false, true);
								}
								if(searchIn){
									writeMateUnmappedPair(bAln, bAlnSeq, *search, searchSeq);
								}else{
									writeUnmappedMateFilteredPair(bAln, bAlnSeq);
								}
							} else {
								if (searchIn) {
									writeThrowAwayUnmappedMate(*search, searchSeq);
									writeTheThrownAwayUnmappedMate(bAln, bAlnSeq);
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
				if(region.getPercInRegion(bAln, refData) >= percInRegion){
					//unpaired read
					++ret.unpaiedReads_;
					if (originalOrientation) {
						writer.openWrite(bamAlnToSeqInfo(bAln));
					} else {
						seqInfo outSeq(bAln.Name, bAln.QueryBases, bAln.Qualities,
								SangerQualOffset);
						//put in the orientation of the output region
						if(region.reverseSrand_){
							outSeq.reverseComplementRead(false, true);
						}
						writer.openWrite(outSeq);
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
		OutOptions orphansSiblingsLocationFileOpt(bib::files::make_path(inOutOpts.out_.outFilename_.parent_path(), "orphansSiblingsStats.tab.txt"));
		orphansSiblingsLocationFileOpt.openFile(orphansSiblingsLocationFile);
		orphansSiblingsLocationFile << "OrphanName\tdistance\tclosetRegions\tRefId\tPosition\tIsMapped\tMateRefId\tMatePosition\tIsMateMapped" << "\n";
	}
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		if(tryToFindOrphansMate){
			//find orphans' mates if possible;
			BamTools::BamReader bReaderMateFinder;
			bReaderMateFinder.Open(inOutOpts.firstName_.string());
			checkBamOpenThrow(bReaderMateFinder, inOutOpts.firstName_.string());
			loadBamIndexThrow(bReaderMateFinder);
			while (bReaderMateFinder.GetNextAlignment(bAln)) {
				if (bAln.IsPaired()) {
					if (alnCache.has(bAln.Name)) {
						auto search = alnCache.get(bAln.Name);
						if (bAln.IsFirstMate() != search->IsFirstMate()) {
							writeTheThrownAwayMate(bAln,
									seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
											SangerQualOffset));
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
								<< "\t" << (bestRegions.empty() ? "none" : bib::conToStr(bestRegions, ","))
								<< "\t" << (search->RefID < 0 ? std::string("*") : refData[search->RefID].RefName)
								<< "\t" << search->Position
								<< "\t" << bib::boolToStr(search->IsMapped())
								<< "\t" << (search->MateRefID < 0 ? std::string("*") : refData[search->MateRefID].RefName)
								<< "\t" << search->MatePosition
								<< "\t" << bib::boolToStr(search->IsMateMapped())
								<< "\n";
					} else {
						orphansSiblingsLocationFile << search->Name
								<< "\t" << "*"
								<< "\t" << "*"
								<< "\t" << (search->RefID < 0 ? std::string("*") : refData[search->RefID].RefName)
								<< "\t" << search->Position
								<< "\t" << bib::boolToStr(search->IsMapped())
								<< "\t" << (search->MateRefID < 0 ? std::string("*") : refData[search->MateRefID].RefName)
								<< "\t" << search->MatePosition
								<< "\t" << bib::boolToStr(search->IsMateMapped())
								<< "\n";
					}
				}
			}
			if(!search->IsMapped()){
				if(search->IsMapped()){
					++ret.orphansUnmapped_;
				}else{
					++ret.unpairedUnMapped_;
				}
				alnCache.remove(name);
				continue;
			}

			auto searchRegion = alnCache.getRegion(name);
			bool searchIn = searchRegion->getPercInRegion(*search, refData)>= percInRegion;
			if (searchIn) {
				if (search->IsPaired()) {
					++ret.orphans_;
				} else {
					++ret.unpaiedReads_;
				}
				if (originalOrientation) {
					writer.openWrite(bamAlnToSeqInfo(*search));
				} else {
					seqInfo searchSeq(search->Name, search->QueryBases, search->Qualities, SangerQualOffset);
					if (searchRegion->reverseSrand_) {
						searchSeq.reverseComplementRead(false, true);
					}
					writer.openWrite(searchSeq);
				}
			}
			alnCache.remove(name);
		}
	}
	return ret;
}

}  // namespace bibseq
