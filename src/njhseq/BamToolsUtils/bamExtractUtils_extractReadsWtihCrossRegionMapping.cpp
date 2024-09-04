//
// Created by Nicholas Hathaway on 9/3/24.
//


#include "bamExtractUtils.hpp"
#include "njhseq/BamToolsUtils/BamAlnsCache.hpp"
#include "njhseq/BamToolsUtils/BamAlnsCacheWithRegion.hpp"
#include "njhseq/IO/SeqIO.h"
#include "njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp"




namespace njhseq {


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

	auto filteredInverseSinglesOpts = SeqIOOptions::genFastqOut(njh::files::prependFileBasename(outOpts.outFilename_, "filteredInverseSingles_"));
	filteredInverseSinglesOpts.out_.transferOverwriteOpts(outOpts);

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
	SeqOutput filteredInverseSinglesWriter(filteredInverseSinglesOpts);


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
	ret.inSoftClipFilteredSingles_= SeqIOOptions::genFastqIn (filteredSoftClipSinglesOpts.getPriamryOutName());


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
			if((len(outSeq) > region.getLen() || region.getLen() < 150 || extractPars.trimToRegion_) ){
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

			if((len(outSeq) > region.getLen() || region.getLen() < 150 || extractPars.trimToRegion_) ){
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
				if(extractPars.renameSingles_){
					bAlnSeq.name_ = njh::pasteAsStr(bAlnSeq.name_, "-", ret.unpairedFailedSoftClip_ + ret.unpaiedReads_);
				}
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

					bool bAlnPassSoftClipAmount = getSoftClipAmount(bAln) < extractPars.softClipHardCutOff_ && getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) < extractPars.softClipPercentageCutOff_;
					bool searchPassSoftClipAmount = getSoftClipAmount(*search) < extractPars.softClipHardCutOff_ && getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_;

					// std::cout << "extractPars.softClipPercentageCutOff_: " << extractPars.softClipPercentageCutOff_ << std::endl;
					// std::cout << "getSoftClipAmount(bAln): " << getSoftClipAmount(bAln) << std::endl;
					// std::cout << "getSoftClipAmount(*search): " << getSoftClipAmount(*search) << std::endl;
					// std::cout << "bAlnPassSoftClipAmount: " << njh::colorBool(bAlnPassSoftClipAmount) << std::endl;
					// std::cout << "searchPassSoftClipAmount: " << njh::colorBool(searchPassSoftClipAmount) << std::endl;
					// std::cout << region.genBedRecordCore().toDelimStrWithExtra() << std::endl;
					// std::cout << "bAln.QueryBases.size()   : " << bAln.QueryBases.size() << std::endl;
					// std::cout << "search->QueryBases.size(): " << search->QueryBases.size() << std::endl;


					seqInfo bAlnSeq(bAln.Name, bAln.QueryBases, bAln.Qualities, SangerQualOffset);

					seqInfo searchSeq(search->Name, search->QueryBases, search->Qualities, SangerQualOffset);

					if (bAln.IsMapped()) {
						bAlnIn = region.getPercInRegion(bAln, refData) >= extractPars.percInRegion_;
					}
					if(search->IsMapped()){
						searchIn = searchRegion->getPercInRegion(*search, refData) >= extractPars.percInRegion_;
					}
					bool balnAllFilters = bAlnIn && bAlnPassAlnSize;
					bool searchAllFilters = searchIn && searchPassAlnSize;

					if((len(bAlnSeq) > region.getLen() || region.getLen() < 150 || extractPars.trimToRegion_) && bAlnIn){

						seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
						GenomicRegion balnRegion(bAln, refData);

						uint32_t startRelative = 0;
						if(balnRegion.start_ < region.start_) {
							startRelative = region.start_ - balnRegion.start_;
						}
//						uint32_t startRelative = region.start_ - balnRegion.start_;
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

					if((len(searchSeq) > searchRegion->getLen() || searchRegion->getLen() < 150 || extractPars.trimToRegion_) && searchIn){

						seqInfo querySeq = bamAlnToSeqInfo(*search, true);

						GenomicRegion balnRegion(*search, refData);
						//std::cout << njh::json::toJson(*search) << std::endl;
//						std::cout << genCigarStr(*search) << std::endl;
//
//						std::cout << "searchRegion->genBedRecordCore().toDelimStrWithExtra(): " << std::endl;
//						std::cout << searchRegion->genBedRecordCore().toDelimStrWithExtra() << std::endl;
//						std::cout << balnRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//						std::cout << "balnRegion.start_   : " << balnRegion.start_ << std::endl;
//						std::cout << "searchRegion->start_: " << searchRegion->start_ << std::endl;
//						std::cout << "balnRegion.start_ >= searchRegion->start_: " << njh::colorBool(balnRegion.start_ >= searchRegion->start_) << std::endl;
						uint32_t startRelative = 0;
						if(balnRegion.start_ < searchRegion->start_) {
							startRelative = searchRegion->start_ - balnRegion.start_;
						}
						uint32_t endRelative = searchRegion->end_ - balnRegion.start_;

						seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
						auto alnInfo = bamAlnToAlnInfoLocal(*search);
						alignCalc::rearrangeLocal(holderSeq.seq_,  querySeq.seq_, '-'	, alnInfo.begin()->second);
						alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0	, alnInfo.begin()->second);

//						std::cout << "searchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
//						holderSeq.outPutSeqAnsi(std::cout);
//						searchSeq.outPutSeqAnsi(std::cout);
//						querySeq.outPutSeqAnsi(std::cout);
//						std::cout << "startRelative: " << startRelative << std::endl;
//						std::cout << "endRelative: " << endRelative << std::endl;
//						std::cout << "searchRegion->start_: " << searchRegion->start_ << std::endl;
//						std::cout << "balnRegion.start_: " << balnRegion.start_ << std::endl;

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

//						if(balnRegion.start_ > searchRegion->start_){
//							exit(1);
//						}
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
//
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
					if(getSoftClipAmount(bAln) > extractPars.softClipHardCutOff_ || getSoftClipAmount(bAln)/static_cast<double>(bAln.QueryBases.size()) > extractPars.softClipPercentageCutOff_){
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

								if((len(outSeq) > region.getLen() || region.getLen() < 150|| extractPars.trimToRegion_) ){
									seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
									GenomicRegion balnRegion(bAln, refData);
//									uint32_t startRelative = region.start_ - balnRegion.start_;
									uint32_t startRelative = 0;
									if (balnRegion.start_ < region.start_) {
										startRelative = region.start_ - balnRegion.start_;
									}
//						uint32_t startRelative = region.start_ - balnRegion.start_;
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
								if(extractPars.renameSingles_){
									outSeq.name_ = njh::pasteAsStr(outSeq.name_, "-", ret.unpairedFailedSoftClip_ + ret.unpaiedReads_);
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
		// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		if (extractPars.tryToFindOrphansMate_) {
			// std::cout << __FILE__ << " " << __LINE__ << std::endl;
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

//
//			std::cout << search->Name << std::endl;
//			std::cout << "\tgetAlnLen(*search): " << getAlnLen(*search) << std::endl;
//			std::cout << "\tnjh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_): " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//			std::cout << "\tsearchIn: " << njh::colorBool(searchIn) << std::endl;
//			std::cout << "\tsearchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
			if (searchIn &&
					getAlnLen(*search) >= extractPars.minAlnMapSize_) {
				// std::cout << __FILE__ << " " << __LINE__ << std::endl;
				// std::cout << "getSoftClipAmount(*search): " << getSoftClipAmount(*search) << std::endl;
				// std::cout << "\tgetSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_: " << njh::colorBool(getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_) << std::endl;

				if(getSoftClipAmount(*search) < extractPars.softClipHardCutOff_ && getSoftClipAmount(*search)/static_cast<double>(search->QueryBases.size()) < extractPars.softClipPercentageCutOff_){
					// std::cout << __FILE__ << " " << __LINE__ << std::endl;
					// std::cout << "extractPars.removeInverseOrphans_: " << njh::colorBool(extractPars.removeInverseOrphans_) << std::endl;
					// std::cout << "search->IsReverseStrand() == search->IsMateReverseStrand(): " << njh::colorBool(search->IsReverseStrand() == search->IsMateReverseStrand()) << std::endl;
					if(extractPars.removeInverseOrphans_ && search->IsReverseStrand() == search->IsMateReverseStrand()) {
						//write filtered orphan
						++ret.orphansFilteredInverse_;
						if(extractPars.writeAll_){
							if(subSamplingFunction()){
								//i'm not sure this sub sampling function should be here.....
								filteredInverseSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
							}
							bWriter.SaveAlignment(*search);
						}
					} else {
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
							// std::cout << __FILE__ << " " << __LINE__ << std::endl;
							if (search->IsPaired()) {
								++ret.orphans_;
							} else {
								++ret.unpaiedReads_;
							}
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
	//
	//						std::cout << search->Name << std::endl;
	//						std::cout << "\tgetAlnLen(*search): " << getAlnLen(*search) << std::endl;
	//						std::cout << "\tnjh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_): " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
	//						std::cout << "\tsearchIn: " << njh::colorBool(searchIn) << std::endl;
	//						std::cout << "\tsearchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
	//
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
//
//				std::cout << search->Name << std::endl;
//				std::cout << "\tgetAlnLen(*search): " << getAlnLen(*search) << std::endl;
//				std::cout << "\tnjh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_): " << njh::colorBool(getAlnLen(*search) >= extractPars.minAlnMapSize_) << std::endl;
//				std::cout << "\tsearchIn: " << njh::colorBool(searchIn) << std::endl;
//				std::cout << "\tsearchRegion->getPercInRegion(*search, refData): " << searchRegion->getPercInRegion(*search, refData) << std::endl;
//
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


//

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

} // namespace njhseq