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
/*
 * reading.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: nick
 */

#include "reading.hpp"
#include "njhseq/objects/BioDataObject/BioDataFileIO.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"

namespace njhseq {

std::vector<std::shared_ptr<GFFCore>> getGFFs(const bfs::path & filename) {
	std::vector<std::shared_ptr<GFFCore>> ret;
	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(filename))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gff = reader.readNextRecord();
	while (nullptr != gff) {
		ret.emplace_back(std::move(gff));
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gff = reader.readNextRecord();
		++count;
	}
	return ret;
}

std::vector<std::shared_ptr<Bed6RecordCore>> getBeds(
		const bfs::path & filename) {
	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(filename))};
	reader.openIn();
	std::vector<std::shared_ptr<Bed6RecordCore>> beds;
	std::shared_ptr<Bed6RecordCore> bed = reader.readNextRecord();
	while (nullptr != bed) {
		beds.emplace_back(std::move(bed));
		bed = reader.readNextRecord();
	}
	return beds;
}

std::vector<std::shared_ptr<Bed3RecordCore>> getBed3s(
		const bfs::path & filename) {
	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(filename))};
	reader.openIn();
	std::vector<std::shared_ptr<Bed3RecordCore>> beds;
	std::shared_ptr<Bed3RecordCore> bed = reader.readNextRecord();
	while (nullptr != bed) {
		beds.emplace_back(std::move(bed));
		bed = reader.readNextRecord();
	}
	return beds;
}

std::vector<std::shared_ptr<RefSeqGeneRecord>> getRefSeqGeneRecords(
		const bfs::path & filename) {
	BioDataFileIO<RefSeqGeneRecord> reader{IoOptions(InOptions(filename))};
	reader.openIn();
	std::vector<std::shared_ptr<RefSeqGeneRecord>> refSeqs;
	std::shared_ptr<RefSeqGeneRecord> refSeq = reader.readNextRecord();
	while (nullptr != refSeq) {
		refSeqs.emplace_back(std::move(refSeq));
		refSeq = reader.readNextRecord();
	}
	return refSeqs;
}


std::vector<std::shared_ptr<RepeatMaskerRecord>> getRepeatMaskerRecords(const bfs::path & filename){
	BioDataFileIO<RepeatMaskerRecord> reader{IoOptions(InOptions(filename))};
	reader.openIn();
	std::vector<std::shared_ptr<RepeatMaskerRecord>> repeatRecords;
	std::shared_ptr<RepeatMaskerRecord> repeatRecord = reader.readNextRecord();
	while (nullptr != repeatRecord) {
		repeatRecords.emplace_back(std::move(repeatRecord));
		repeatRecord = reader.readNextRecord();
	}
	return repeatRecords;
}

std::vector<std::shared_ptr<TandemRepeatFinderRecord>> getTandemRepeatFinderRecords(const bfs::path & filename){
	BioDataFileIO<TandemRepeatFinderRecord> reader{IoOptions(InOptions(filename))};
	reader.openIn();
	std::vector<std::shared_ptr<TandemRepeatFinderRecord>> repeatRecords;
	std::string currentSeqName = "";
	std::string line = "";
	while (njh::files::crossPlatGetline(*reader.inFile_, line)) {
		if (njh::beginsWith(line, "Sequence")) {
			auto toks = njh::tokenizeString(line, "whitespace");
			if (toks.size() < 2) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error in processing line " << line
						<< ", was expecting at least two values separated by white sapce but found "
						<< toks.size() << " instead " << "\n";
				throw std::runtime_error { ss.str() };
			}
			currentSeqName = toks[1];
		}
		//skip lines that don't start with a digit or is a blank line, it seems to be the only indicator that the file is on a data line, this format is stupid
		if("" == line || allWhiteSpaceStr(line) || !::isdigit(line.front())){
			continue;
		}
		std::shared_ptr<TandemRepeatFinderRecord> record = std::make_shared<TandemRepeatFinderRecord>(line);
		record->setSeqName(currentSeqName);
		repeatRecords.emplace_back(record);
	}
	return repeatRecords;
}


void checkPositionSortedBedThrow(const bfs::path & bedFnp,
		const std::string & funcName) {
	Bed3RecordCore bedRecord;
	Bed3RecordCore bedRecordNext;
	BioDataFileIO<Bed3RecordCore> bedReader { IoOptions(InOptions(bedFnp)) };
	bedReader.openIn();
	bedReader.readNextRecord(bedRecord);
	uint32_t record = 2;
	VecStr alreadySeenChroms;
	while (bedReader.readNextRecord(bedRecordNext)) {
		if (bedRecord.chrom_ == bedRecordNext.chrom_) {
			if (bedRecord.chromStart_ > bedRecordNext.chromStart_) {
				std::stringstream ss;
				ss << funcName << ", error record " << record - 1 << " starts after "
						<< record << "\n";
				ss << "Record " << record - 1 << ": " << bedRecord.toDelimStr() << "\n";
				ss << "Record " << record << ": " << bedRecordNext.toDelimStr() << "\n";
				throw std::runtime_error { ss.str() };
			}
		} else {
			alreadySeenChroms.emplace_back(bedRecord.chrom_);
			if (njh::in(bedRecordNext.chrom_, alreadySeenChroms)) {
				std::stringstream ss;
				ss << funcName << ", error record " << record << " has chrom "
						<< bedRecordNext.chrom_
						<< ", which was already seen earlier in file, should be sorted by chrom as well as positione"
						<< "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		bedRecord = bedRecordNext;
		++record;
	}
}



intersectBedLocsWtihGffRecordsPars::intersectBedLocsWtihGffRecordsPars(){}
intersectBedLocsWtihGffRecordsPars::intersectBedLocsWtihGffRecordsPars(const bfs::path & gffFnp) :
		gffFnp_(gffFnp) {
}
intersectBedLocsWtihGffRecordsPars::intersectBedLocsWtihGffRecordsPars(const bfs::path & gffFnp,
		const VecStr & extraAttributes, const VecStr & selectFeatures) :
		gffFnp_(gffFnp), extraAttributes_(extraAttributes), selectFeatures_(
				selectFeatures) {
}


Json::Value intersectBedLocsWtihGffRecordsPars::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));

	ret["gffFnp_"] = njh::json::toJson(gffFnp_);
	ret["extraAttributes_"] = njh::json::toJson(extraAttributes_);
	ret["selectFeatures_"] = njh::json::toJson(selectFeatures_);

	return ret;
}




void checkPositionSortedNoOverlapsBedThrow(const bfs::path & bedFnp,
		const std::string & funcName){
	Bed3RecordCore bedRecord;
	Bed3RecordCore bedRecordNext;
	BioDataFileIO<Bed3RecordCore> bedReader { IoOptions(InOptions(bedFnp)) };
	bedReader.openIn();
	bedReader.readNextRecord(bedRecord);
	uint32_t record = 2;
	VecStr alreadySeenChroms;
	while (bedReader.readNextRecord(bedRecordNext)) {
		if (bedRecord.chrom_ == bedRecordNext.chrom_) {
			if (bedRecord.chromStart_ > bedRecordNext.chromStart_) {
				std::stringstream ss;
				ss << funcName << ", error record " << record - 1 << " starts after "
						<< record << "\n";
				ss << "Record " << record - 1 << ": " << bedRecord.toDelimStr() << "\n";
				ss << "Record " << record << ": " << bedRecordNext.toDelimStr() << "\n";
				throw std::runtime_error { ss.str() };
			} else if (GenomicRegion(bedRecordNext).overlaps(GenomicRegion(bedRecord))){
				std::stringstream ss;
				ss << funcName << ", error record " << record - 1 << " overlaps with "
						<< record << "\n";
				ss << "Record " << record - 1 << ": " << bedRecord.toDelimStr() << "\n";
				ss << "Record " << record << ": " << bedRecordNext.toDelimStr() << "\n";
				throw std::runtime_error { ss.str() };
			}
		} else {
			alreadySeenChroms.emplace_back(bedRecord.chrom_);
			if (njh::in(bedRecordNext.chrom_, alreadySeenChroms)) {
				std::stringstream ss;
				ss << funcName << ", error record " << record << " has chrom "
						<< bedRecordNext.chrom_
						<< ", which was already seen earlier in file, should be sorted by chrom as well as positione"
						<< "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		bedRecord = bedRecordNext;
		++record;
	}
}

//Json::Value intersectBedLocsWtihGffRecords(
//		const std::vector<std::shared_ptr<Bed3RecordCore>> & beds,
//		const intersectBedLocsWtihGffRecordsPars & pars){
//
//	Json::Value ret;
//	for (const auto & inputRegion : beds) {
//		inputRegion->extraFields_.emplace_back("");
//	}
//	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed3RecordCore>>> bedsByChrome;
//	BioDataFileIO<GFFCore> reader { IoOptions(InOptions(pars.gffFnp_)) };
//	reader.openIn();
//	uint32_t count = 0;
//	std::string line = "";
//	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
//	for (const auto & b : beds) {
//		bedsByChrome[b->chrom_].emplace_back(b);
//	}
//	while (nullptr != gRecord) {
//		if (pars.selectFeatures_.empty() || njh::in(gRecord->type_, pars.selectFeatures_)) {
//			auto gRegion = GenomicRegion(*gRecord);
//			for (const auto & inputRegion : bedsByChrome[gRegion.chrom_]) {
//				if (GenomicRegion(*inputRegion).overlaps(gRegion)) {
//					if("" != inputRegion->extraFields_.back()){
//						inputRegion->extraFields_.back().append(",");
//					}
//					inputRegion->extraFields_.back().append("[");
//					inputRegion->extraFields_.back().append(
//							"ID=" + gRecord->getAttr("ID") + ";");
//					if(pars.selectFeatures_.empty() || 1 != pars.selectFeatures_.size()){
//						inputRegion->extraFields_.back().append("feature=" + gRecord->type_ + ";");
//					}
//					if (!ret.isMember(gRecord->getAttr("ID"))) {
//						ret[gRecord->getAttr("ID")] = gRecord->toJson();
//					}
//					for (const auto & attr : pars.extraAttributes_) {
//						if (gRecord->hasAttr(attr)) {
//							inputRegion->extraFields_.back().append(
//									attr + "=" + gRecord->getAttr(attr) + ";");
//						} else {
//							inputRegion->extraFields_.back().append(
//									attr + "=" + "NA" + ";");
//						}
//					}
//					inputRegion->extraFields_.back().append("]");
//				}
//			}
//		}
//		bool end = false;
//		while ('#' == reader.inFile_->peek()) {
//			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
//				end = true;
//				break;
//			}
//			njh::files::crossPlatGetline(*reader.inFile_, line);
//		}
//		if (end) {
//			break;
//		}
//		gRecord = reader.readNextRecord();
//		++count;
//	}
//	return ret;
//}

}  // namespace njhseq

