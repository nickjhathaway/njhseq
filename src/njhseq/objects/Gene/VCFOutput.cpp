//
// Created by Nicholas Hathaway on 1/30/24.
//

#include "VCFOutput.hpp"

namespace njhseq {

uint32_t VCFOutput::VCFRecord::getNumberOfAlleles() const {
	return 1 + alts_.size();
}





Json::Value VCFOutput::InfoEntry::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["id_"] = njh::json::toJson(id_);
	ret["number_"] = njh::json::toJson(number_);
	ret["type_"] = njh::json::toJson(type_);
	ret["description_"] = njh::json::toJson(description_);
	ret["source_"] = njh::json::toJson(source_);
	ret["version_"] = njh::json::toJson(version_);

	return ret;
}

Json::Value VCFOutput::FormatEntry::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["id_"] = njh::json::toJson(id_);
	ret["number_"] = njh::json::toJson(number_);
	ret["type_"] = njh::json::toJson(type_);
	ret["description_"] = njh::json::toJson(description_);
	return ret;
}

Json::Value VCFOutput::FilterEntry::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["id_"] = njh::json::toJson(id_);
	ret["description_"] = njh::json::toJson(description_);
	return ret;
}

Json::Value VCFOutput::ContigEntry::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["id_"] = njh::json::toJson(id_);
	ret["length_"] = njh::json::toJson(length_);
	ret["assembly_"] = njh::json::toJson(assembly_);
	ret["md5_"] = njh::json::toJson(md5_);
	ret["species_"] = njh::json::toJson(species_);
	ret["otherKeysValues_"] = njh::json::toJson(otherKeysValues_);

	return ret;
}

GenomicRegion VCFOutput::VCFRecord::genRegion() const {
	uint32_t start = pos_ - 1;
	uint32_t end = pos_;
	if (ref_.size() == 1 && !std::all_of(alts_.begin(), alts_.end(), [](const std::string& alt) {
		return alt.size() == 1;
	})) {
		//insertion, will give the region right before and right after the insertion
		end += 1;
	} else if (ref_.size() > 1) {
		//should only be deletions
		start += 1; //increase start to go actual deleted base
		end = start - 1 + ref_.size();
	}
	std::string uid = njh::pasteAsStr(chrom_, "-", start, "-", end);
	if(id_ != ".") {
		uid = id_;
	}
	GenomicRegion ret(uid,chrom_,start, end, false);
	ret.meta_ = info_;
	ret.meta_.addMeta("ref", ref_);
	ret.meta_.addMeta("alts", njh::conToStr(alts_, ","));
	ret.meta_.addMeta("qual", qual_);
	ret.meta_.addMeta("filter", filter_);

	return ret;
}

Json::Value VCFOutput::VCFRecord::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["chrom_"] = njh::json::toJson(chrom_);
	ret["pos_"] = njh::json::toJson(pos_);
	ret["id_"] = njh::json::toJson(id_);
	ret["ref_"] = njh::json::toJson(ref_);
	ret["alts_"] = njh::json::toJson(alts_);
	ret["qual_"] = njh::json::toJson(qual_);
	ret["filter_"] = njh::json::toJson(filter_);
	ret["info_"] = njh::json::toJson(info_);
	ret["sampleFormatInfos_"] = njh::json::toJson(sampleFormatInfos_);
	return ret;
}


void VCFOutput::sortRecords() {
	njh::sort(records_, [](const VCFRecord & r1, const VCFRecord & r2) {
		if(r1.chrom_ == r2.chrom_) {
			if(r1.pos_ == r2.pos_) {
				return r1.ref_ < r2.ref_;
			} else {
				return r1.pos_ < r2.pos_;
			}
		} else {
			return r1.chrom_ < r2.chrom_;
		}
	});
}

void VCFOutput::writeOutHeaderFieldsOtherThanFormat(std::ostream & vcfOut) const {
	//write out contigs
	for (const auto&contigKey: contigEntries_) {
		const auto & contig = contigKey.second;
		vcfOut <<"##contig=<"
		<< "ID=" << contig.id_ << ","
		<< "length=" << contig.length_;
		if(!contig.md5_.empty()) {
			vcfOut << "," << "md5=" << contig.md5_;
		}
		if(!contig.assembly_.empty()) {
			if(contig.assembly_.front() == '"' && contig.assembly_.back() == '"') {
				vcfOut << "," << "assembly=" << contig.assembly_ << "";
			}else {
				vcfOut << "," << "assembly=\"" << contig.assembly_ << "\"";
			}
		}
		if(!contig.species_.empty()) {
			if(contig.species_.front() == '"' && contig.species_.back() == '"') {
				vcfOut << "," << "species=" << contig.species_ << "";
			} else {
				vcfOut << "," << "species=\"" << contig.species_ << "\"";
			}
		}
		for(const auto & others : contig.otherKeysValues_) {
			if((std::string::npos != others.second.find(',') || njh::strHasWhitesapce(others.second)) && !(others.second.front() == '"' && others.second.back() == '"') ) {
				vcfOut << "," << others.first << "=" << "\"" << others.second << "\"";
			} else {
				vcfOut << "," << others.first << "=" << others.second;
			}
		}
		vcfOut << ">" << std::endl;
	}
	//write out infos
	for (const auto&infoKey: infoEntries_) {
		const auto & info  = infoKey.second;
		vcfOut <<"##INFO=<"
		<< "ID=" << info.id_
		<< ","<< "Number=" << info.number_
		<< "," << "Type=" << info.type_;
		if(info.description_.front() == '"' && info.description_.back() == '"') {
			vcfOut << "," << "Description=" << info.description_ << "";
		} else {
			vcfOut << "," << "Description=\"" << info.description_ << "\"";
		}
		if(!info.source_.empty()) {
			if(info.source_.front() == '"' && info.source_.back() == '"') {
				vcfOut << "," << "Source=" << info.source_ << "";
			} else {
				vcfOut << "," << "Source=\"" << info.source_ << "\"";
			}
		}
		if(!info.version_.empty()) {
			if(info.version_.front() == '"' && info.version_.back() == '"') {
				vcfOut << "," << "Version=" << info.version_ << "";
			} else {
				vcfOut << "," << "Version=\"" << info.version_ << "\"";
			}
		}
		vcfOut << ">"
		<< std::endl;
	}
	//write out filters
	for (const auto& info: filterEntries_) {
		vcfOut << "##FILTER=<"
				<< "ID=" << info.id_;
		if (info.description_.front() == '"' && info.description_.back() == '"') {
			vcfOut << "," << "Description=" << info.description_ << "";
		} else {
			vcfOut << "," << "Description=\"" << info.description_ << "\"";
		}
		vcfOut << ">" << std::endl;
	}

	//write out other metas
	for(const auto & otherMeta : otherHeaderMetaFields_) {
		vcfOut << "##" << otherMeta.first << "=<";
		bool first = true;
		if(otherMeta.second.containsMeta("ID")) {
			first = false;
			vcfOut << "ID" << "=" << otherMeta.second.getMeta("ID");
		}
		for(const auto & valuePairs : otherMeta.second.meta_) {
			if(valuePairs.first == "ID") {
				continue;
			}
			if(!first) {
				vcfOut << ",";
			}
			if((std::string::npos != valuePairs.second.find(',') || njh::strHasWhitesapce(valuePairs.second)) && !(valuePairs.second.front() == '"' && valuePairs.second.back() == '"') ) {
				vcfOut << valuePairs.first << "=" << "\"" << valuePairs.second << "\"";
			} else {
				vcfOut << valuePairs.first << "=" << valuePairs.second;
			}
			first = false;
		}
		vcfOut << ">" << std::endl;
	}

	//write out key=value pairs
	for(const auto & other : otherHeaderValuePairs_) {
		vcfOut << "##" << other.first << "=" << other.second << std::endl;
	}
}


void VCFOutput::writeOutFixedOnly(std::ostream&vcfOut, const std::vector<GenomicRegion> & selectRegions) const {
	vcfOut << "##fileformat=" << vcfFormatVersion_ << std::endl;
	writeOutHeaderFieldsOtherThanFormat(vcfOut);
	//write out
	//vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
	vcfOut << njh::conToStr(getSubVector(headerNonSampleFields_,0, 8), "\t") << std::endl;
	for(const auto & rec : records_) {
		// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(!selectRegions.empty()) {
			// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			bool overlaps = false;
			for(const auto & reg : selectRegions) {
				if(reg.overlaps(rec.genRegion())) {
					overlaps = true;
					break;
				}
			}
			// std::cout << "overlaps: " << njh::colorBool(overlaps) << std::endl;
			// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if(!overlaps) {
				// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				continue;
			}
			// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		vcfOut << rec.chrom_
		<< "\t" << rec.pos_
		<< "\t" << rec.id_
		<< "\t" << rec.ref_
		<< "\t" << njh::conToStr(rec.alts_, ",")
		<< "\t" << rec.qual_
		<< "\t" << rec.filter_;
		std::string infoOut;
		for (const auto & infoKey: infoEntries_) {
			const auto & info = infoKey.second;
			if(!infoOut.empty()) {
				infoOut +=";";
			}
			infoOut += info.id_ + "=" + rec.info_.getMeta(info.id_);
		}
		vcfOut << "\t" << infoOut;
		vcfOut << std::endl;
	}
}


void VCFOutput::writeOutFixedAndSampleMeta(std::ostream& vcfOut, const std::vector<GenomicRegion>& selectRegions) const {
	vcfOut << "##fileformat=" << vcfFormatVersion_ << std::endl;
	writeOutHeaderFieldsOtherThanFormat(vcfOut);
	//write out formats
	for (const auto & infoKey: formatEntries_) {
		const auto & info = infoKey.second;
		vcfOut <<"##FORMAT=<"
		<< "ID=" << info.id_
		<< ","<< "Number=" << info.number_
		<< "," << "Type=" << info.type_;
		if(info.description_.front() == '"' && info.description_.back() == '"') {
			vcfOut << "," << "Description=" << info.description_ << "";
		} else {
			vcfOut << "," << "Description=\"" << info.description_ << "\"";
		}
		vcfOut << ">" << std::endl;
	}
	//check samples

	//vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	vcfOut << njh::conToStr(headerNonSampleFields_, "\t") << "\t" << njh::conToStr(samples_, "\t") << std::endl;
	if(!records_.empty()) {
		const auto firstSetOfSamples = njh::vecToSet(getVectorOfMapKeys(records_.front().sampleFormatInfos_));
		const auto firstSetOfSamplesVec = VecStr(firstSetOfSamples.begin(), firstSetOfSamples.end());
		std::set<std::string> headerSamplesSet(samples_.begin(), samples_.end());
		if(firstSetOfSamples != headerSamplesSet) {
			std::vector<std::string> uniqueTo1;
			std::vector<std::string> uniqueTo2;
			std::vector<std::string> inBoth;

			njh::decompose_sets(
				headerSamplesSet.begin(), headerSamplesSet.end(),
				firstSetOfSamples.begin(), firstSetOfSamples.end(),
				std::back_inserter(uniqueTo1),
				std::back_inserter(uniqueTo2),
				std::back_inserter(inBoth));
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "samples in records don't match the header samples"  << "\n";
			ss << "header samples: " << njh::conToStr(headerSamplesSet, "\t") << "\n";
			ss << "record samples: " << njh::conToStr(firstSetOfSamplesVec, "\t") << "\n";
			ss << "samples only in header: " << njh::conToStr(uniqueTo1, "\t") << "\n";
			ss << "samples only in record: " << njh::conToStr(uniqueTo2, "\t") << "\n";
			throw std::runtime_error{ss.str()};
		}
		for(const auto & rec : records_) {

			auto currentSetOfSamples = njh::vecToSet(getVectorOfMapKeys(rec.sampleFormatInfos_));
			if(firstSetOfSamples != currentSetOfSamples) {
				std::vector<std::string> uniqueTo1;
				std::vector<std::string> uniqueTo2;
				std::vector<std::string> inBoth;

				njh::decompose_sets(firstSetOfSamples.begin(), firstSetOfSamples.end(),
					currentSetOfSamples.begin(), currentSetOfSamples.end(),
					std::back_inserter(uniqueTo1),
					std::back_inserter(uniqueTo2),
					std::back_inserter(inBoth));
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "samples different for record: " << rec.chrom_ << " " << rec.pos_ << " " << rec.ref_ << "\n";
				ss << "first set of samples: " << njh::conToStr(firstSetOfSamples, ",") << "\n";
				ss << "current set of samples: " << njh::conToStr(currentSetOfSamples, ",") << "\n";
				ss << "samples only in first set  : " << njh::conToStr(uniqueTo1, "\t") << "\n";
				ss << "samples only in current set: " << njh::conToStr(uniqueTo2, "\t") << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
	}
	std::string formatOut;
	for (const auto & formatKey: formatEntries_) {
		const auto & format = formatKey.second;
		if(!formatOut.empty()) {
			formatOut +=":";
		}
		formatOut += format.id_;
	}
	if(!records_.empty()) {
		for(const auto & rec : records_) {
			if(!selectRegions.empty()) {
				bool overlaps = false;
				for(const auto & reg : selectRegions) {
					if(reg.overlaps(rec.genRegion())) {
						overlaps = true;
						break;
					}
				}
				if(!overlaps) {
					continue;
				}
			}
			vcfOut << rec.chrom_
			<< "\t" << rec.pos_
			<< "\t" << rec.id_
			<< "\t" << rec.ref_
			<< "\t" << njh::conToStr(rec.alts_, ",")
			<< "\t" << rec.qual_
			<< "\t" << rec.filter_;
			std::string infoOut;
			for (const auto & infoKey: infoEntries_) {
				const auto & info = infoKey.second;
				if(!infoOut.empty()) {
					infoOut +=";";
				}
				infoOut += info.id_ + "=" + rec.info_.getMeta(info.id_);
			}
			vcfOut << "\t" << infoOut;
			vcfOut << "\t" << formatOut;
			for(const auto & sampleName : samples_) {
				const auto & sample = rec.sampleFormatInfos_.at(sampleName);
				std::string formatOutForSample;
				for (const auto & formatKey: formatEntries_) {
					const auto & format = formatKey.second;
					if(!formatOutForSample.empty()) {
						formatOutForSample +=":";
					}
					formatOutForSample += sample.getMeta(format.id_);
				}
				vcfOut << "\t" << formatOutForSample;
			}
			vcfOut << std::endl;
		}
	}
}
size_t VCFOutput::expectedColumnNumber() const {
	return samples_.size() + headerNonSampleFields_.size();
}

Json::Value VCFOutput::headerToJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["otherHeaderValuePairs_"] = njh::json::toJson(otherHeaderValuePairs_);
	ret["otherHeaderMetaFields_"] = njh::json::toJson(otherHeaderMetaFields_);

	ret["contigEntries_"] = njh::json::toJson(contigEntries_);
	ret["filterEntries_"] = njh::json::toJson(filterEntries_);
	ret["formatEntries_"] = njh::json::toJson(formatEntries_);
	ret["infoEntries_"] = njh::json::toJson(infoEntries_);
	ret["vcfFormatVersion_"] = njh::json::toJson(vcfFormatVersion_);

	ret["headerNonSampleFields_"] = njh::json::toJson(headerNonSampleFields_);
	ret["samples_"] = njh::json::toJson(samples_);


	return ret;
}

void VCFOutput::addInBlnaksForAnyMissingSamples(const std::set<std::string>& samples) {
	std::set<std::string> missingSamples;
	for(auto & record : records_) {
		for(const auto & samp : samples) {
			if(njh::notIn(samp, record.sampleFormatInfos_)) {
				MetaDataInName emptyMeta;
				for(const auto & format : formatEntries_) {
					uint32_t numberOfVals = 1;
					if(njh::strAllDigits(format.second.number_)) {
						numberOfVals = njh::StrToNumConverter::stoToNum<uint32_t>(format.second.number_);
					}else if(format.second.number_ == "R") {
						numberOfVals = record.getNumberOfAlleles();
					}else if(format.second.number_ == "A") {
						numberOfVals = record.alts_.size();
					}
					emptyMeta.addMeta(format.first, njh::conToStr(VecStr{numberOfVals, "."}, ",") );
				}
				record.sampleFormatInfos_.emplace(samp, emptyMeta);
				missingSamples.emplace(samp);
			}
		}
	}
	njh::addVecToSet(samples_, missingSamples);
	samples_ = VecStr(missingSamples.begin(), missingSamples.end());
}




VCFOutput::VCFRecord VCFOutput::processRecordLineForFixedData(const std::string & line) const {
	VCFRecord rec;
	auto toks = tokenizeString(line, "\t");
	if(toks.size() != expectedColumnNumber()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "expected " << expectedColumnNumber() <<" not " << toks.size() << "\n";
		ss << "error for line: " << line << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(toks.size() < 8) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "line should be at least 8 columns, not: " << toks.size() << "\n";
		ss << "error for line: " << line << "\n";
		throw std::runtime_error{ss.str()};
	}
	rec.chrom_ = toks[0];
	rec.pos_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
	rec.id_ = toks[2];
	rec.ref_ = toks[3];
	rec.alts_ = tokenizeString(toks[4], ",");
	rec.qual_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);
	rec.filter_ = toks[6];

	//info field
	VecStr warningsInfoField;

	if(std::string::npos != toks[7].find(':')) {
		warningsInfoField.emplace_back("info field can't have :");
	}
	if(njh::strHasWhitesapce( toks[7])) {
		warningsInfoField.emplace_back("info field can't whitespace");
	}
	if(!warningsInfoField.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "encoutner the following errors when processing info field for line: " << "\n";
		ss << "errors: " << njh::conToStr(warningsInfoField, ",") << "\n";
		ss << line << "\n";
		throw std::runtime_error{ss.str()};
	}

	auto infoToks = tokenizeString(toks[7], ";");
	for(const auto & infoTok : infoToks) {
		const auto equalSignPos = infoTok.find("=");
		if(equalSignPos == std::string::npos || equalSignPos == 0 || equalSignPos +1 >= infoTok.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "info toks should have an equal sign separating values" << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto key = infoTok.substr(0, equalSignPos);
		auto val = infoTok.substr(equalSignPos + 1);

		if(!njh::in(key, infoEntries_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no info entry to define " << key  << " options are " << njh::conToStr(njh::getVecOfMapKeys(infoEntries_)) << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto valCount = 1 + countOccurences(val, ",");
		if(infoEntries_.at(key).number_ == "A") {
			if(valCount != rec.alts_.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << rec.alts_.size() << " but has " << valCount << " instead " << "\n";
				ss << "key: " << key << "\n";
				ss << "val: " << val << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error{ss.str()};
			}
		} else if (infoEntries_.at(key).number_ == "R") {
			if(valCount != rec.alts_.size() + 1) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << rec.alts_.size() + 1 << " but has " << valCount << " instead " << "\n";
				ss << "key: " << key << "\n";
				ss << "val: " << val << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error{ss.str()};
			}
		} else if (njh::strAllDigits(infoEntries_.at(key).number_)) {
			if(njh::pasteAsStr(valCount) != infoEntries_.at(key).number_) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << infoEntries_.at(key).number_ << " but has " << valCount << " instead " << "\n";
				ss << "key: " << key << "\n";
				ss << "val: " << val << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		rec.info_.addMeta(key, val);
	}
	return rec;
}

VCFOutput::VCFRecord VCFOutput::processRecordLineForFixedDataAndSampleMetaData(const std::string & line) const {


	auto rec = processRecordLineForFixedData(line);
	auto toks = tokenizeString(line, "\t");
	//safety checks already done above
	//format
	if(toks.size() > 8) {
		auto formatToks = tokenizeString(toks[8], ":");
		VecStr missingFormat;
		for(const auto & f : formatToks) {
			if(!njh::in(f, formatEntries_)) {
				missingFormat.emplace_back(f);
			}
		}
		if(!missingFormat.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " missing format info on the following formats in this line: " << njh::conToStr(missingFormat, ",") << ", options: " << njh::conToStr(njh::getVecOfMapKeys(formatEntries_))<< "\n";
			ss << "line: " << line << "\n";
			throw std::runtime_error{ss.str()};
		}

		if(toks.size() > 9) {
			//process samples
			for(const auto pos : iter::range(9UL, toks.size())) {
				auto sampleToks = tokenizeString(toks[pos], ":");
				if(sampleToks.size() != formatToks.size()) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "sample info size, " <<  sampleToks.size() << ", doesn't match the expected number of " << formatToks.size() << "\n";
					ss << "sample info: " << toks[pos] << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto sampleName = samples_[pos - 9];
				MetaDataInName sampleInfo;
				for(const auto & e : iter::enumerate(sampleToks)) {
					auto valCount = 1 + countOccurences(e.element, ",");
					const auto & key = formatToks[e.index];
					const auto & val = e.element;
					if(formatEntries_.at(key).number_ == "A") {
						if(valCount != rec.alts_.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << rec.alts_.size() << " but has " << valCount << " instead " << "\n";
							ss << "key: " << key << "\n";
							ss << "val: " << val << "\n";
							ss << "line: " << line << "\n";
							throw std::runtime_error{ss.str()};
						}
					} else if (formatEntries_.at(key).number_ == "R") {
						if(valCount != rec.alts_.size() + 1) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << rec.alts_.size() + 1 << " but has " << valCount << " instead " << "\n";
							ss << "key: " << key << "\n";
							ss << "val: " << val << "\n";
							ss << "line: " << line << "\n";
							throw std::runtime_error{ss.str()};
						}
					} else if (njh::strAllDigits(formatEntries_.at(key).number_)) {
						if(njh::pasteAsStr(valCount) != formatEntries_.at(key).number_) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << formatEntries_.at(key).number_ << " but has " << valCount << " instead " << "\n";
							ss << "key: " << key << "\n";
							ss << "val: " << val << "\n";
							ss << "line: " << line << "\n";
							throw std::runtime_error{ss.str()};
						}
					}
					sampleInfo.addMeta(formatToks[e.index], e.element);
				}
				rec.sampleFormatInfos_.emplace(sampleName, sampleInfo);
			}
		}
	}

	return rec;
}


void VCFOutput::addInRecordsFixedDataFromFile(std::istream & in) {
	std::string line;
	// uint32_t count = 0;
	while(njh::files::crossPlatGetline(in, line)) {
		if(line.front() != '#') {
			// std::cout << count++ << std::endl;
			records_.emplace_back(processRecordLineForFixedData(line));
		}
	}
}

void VCFOutput::addInRecordsFromFile(std::istream & in) {
	std::string line;
	// uint32_t count = 0;
	while(njh::files::crossPlatGetline(in, line)) {
		if(line.front() != '#') {
			// std::cout << count++ << std::endl;
			records_.emplace_back(processRecordLineForFixedDataAndSampleMetaData(line));
		}
	}
}



VCFOutput VCFOutput::readInHeader(const bfs::path & fnp) {
	njh::files::checkExistenceThrow(fnp, __PRETTY_FUNCTION__);
	if(njh::files::isFileEmpty(fnp)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << fnp << " is empty" << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto firstLine = njh::files::getFirstLine(fnp);
	std::string fileFormatCheck = "##fileformat=VCFv";
	if(!njh::beginsWith(firstLine, fileFormatCheck)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << fnp << " should start with " << fileFormatCheck << "\n";
		throw std::runtime_error{ss.str()};
	} else if(firstLine.size() <= fileFormatCheck.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "should have more than just: " << firstLine << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto version = firstLine.substr(firstLine.find_first_not_of("##fileformat="));
	VCFOutput ret;
	ret.vcfFormatVersion_ = version;
	InputStream input(fnp);
	std::string line;
	while(njh::files::crossPlatGetline(input, line)) {
		if(!njh::beginsWith(line, "##")) {
			if(njh::beginsWith(line, "#CHROM")) {
				//header file
				auto toks = tokenizeString(line, "\t");
				if(toks.size() < 8) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << " header fields must at least be 8 fields, not: " << toks.size() << ", error in line: "  << "\n";
					ss << line << "\n";
					throw std::runtime_error{ss.str()};
				}
				//strict checking
				VecStr warnings;
				if("#CHROM" != toks[0]) {
					warnings.emplace_back("Field 1 must be #CHROM");
				}
				if("POS" != toks[1]) {
					warnings.emplace_back("Field 2 must be POS");
				}
				if("ID" != toks[2]) {
					warnings.emplace_back("Field 3 must be ID");
				}
				if("REF" != toks[3]) {
					warnings.emplace_back("Field 4 must be REF");
				}
				if("ALT" != toks[4]) {
					warnings.emplace_back("Field 5 must be ALT");
				}
				if("QUAL" != toks[5]) {
					warnings.emplace_back("Field 6 must be QUAL");
				}
				if("FILTER" != toks[6]) {
					warnings.emplace_back("Field 7 must be FILTER");
				}
				if("INFO" != toks[7]) {
					warnings.emplace_back("Field 8 must be INFO");
				}
				if(toks.size() >8) {
					if("FORMAT" != toks[8]) {
						warnings.emplace_back("Field 9 must be FORMAT");
					}
					ret.headerNonSampleFields_ = getSubVector(toks, 0, 9);
				} else {
					ret.headerNonSampleFields_ = getSubVector(toks, 0, 8);
				}

				if(toks.size() >9) {
					ret.samples_ = getSubVector(toks, 9);
					std::unordered_map<std::string, uint32_t> counts;
					for(const auto & sample : ret.samples_) {
						++counts[sample];
					}
					for(const auto & c : counts) {
						if(c.second > 1) {
							warnings.emplace_back(njh::pasteAsStr("can't duplicate sample names, sample ", c.first, "were found with counts: ", c.second));
						}
					}
				}
				if(!warnings.empty()) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "error, found the following warnings when processing header line: " << "\n";
					ss << line << "\n";
					ss << njh::conToStr(warnings, "\n") << "\n";
					throw std::runtime_error{ss.str()};
				}
			}
			break;
		}
		if(!njh::beginsWith(line, "##fileformat=VCFv")) {
			if(std::string::npos == line.find('=')) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "every line in header should have at least one =, error for line: "  << "\n";
				ss << line << "\n";
				throw std::runtime_error{ss.str()};
			}
			if(std::string::npos == line.find("=<")) {
				//not a meta filed, just a key=value pair
				//split at first found equal sign
				auto pos = line.find_first_of('=');
				if(pos + 1 == line.size()) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << " = shouldn't be at the very end of the line, error for line: " << "\n";
					ss << line << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::string key = line.substr(2, pos -2);
				std::string value = line.substr(pos + 1);
				ret.otherHeaderValuePairs_.emplace(key, value);
			} else {
				auto pos =  line.find("=<");
				if(2 == pos) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "=< shouldn't come right after the ##, error for line: " << "\n";
					ss << line << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(line.back() != '>') {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "in metafield headers line should always end with >, error for line: " << "\n";
					ss << line << "\n";
					throw std::runtime_error{ss.str()};
				}
				const auto metaField = line.substr(2, pos - 2);
				auto restStart = pos + 2;
				auto restEnd = line.size() - 1;
				auto rest = line.substr(restStart, restEnd - restStart);

				//process for quotation marks
				auto numberOfQuotes = countOccurences(rest, "\"");
				if(numberOfQuotes % 2 != 0) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "there should be an even number of \", error for line: "  << "\n";
					ss << "line: " << line << "\n";
					ss << "processed_portion: " << rest << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(numberOfQuotes > 0 ) {
					auto allCommaPositions = findOccurences(rest, ",");
					std::vector<size_t> commaInQuotesPositions;
					auto currentQuote = rest.find_first_of('"');
					bool inbetweenQuote = true;
					auto nextQuote = rest.find_first_of('"', currentQuote + 1);
					while(std::string::npos != nextQuote) {
						if(inbetweenQuote) {
							//process if inbetween quotes
							for(const auto commaPos : allCommaPositions) {
								if(commaPos > currentQuote && commaPos < nextQuote) {
									commaInQuotesPositions.emplace_back(commaPos);
								}
							}
						}
						inbetweenQuote = !inbetweenQuote; //toggle inbetween quote;
						currentQuote = nextQuote;
						if(nextQuote + 1 >= rest.size()) {
							nextQuote = std::string::npos;
						} else {
							nextQuote = rest.find_first_of('"', nextQuote + 1);
						}
					}
					for(auto commaPos : iter::reversed(commaInQuotesPositions)) {
						rest.replace(commaPos, 1, std::string("COMMA_IN_BETWEEN_QUOTES"));
					}
				}
				auto toks = njh::tokenizeString(rest, ",");
				std::unordered_map<std::string, std::string> valuePairs;
				for(const auto & tok : toks) {
					auto equalPos = tok.find_first_of('=');
					auto key = tok.substr(0, equalPos	);
					auto value = tok.substr(equalPos + 1);
					value = njh::replaceString(value, "COMMA_IN_BETWEEN_QUOTES", ",");
					if(njh::in(key, valuePairs)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "error, already have key: " << key << ", error for line: "  << "\n";
						ss << line << "\n";
						throw std::runtime_error{ss.str()};
					}
					valuePairs[key] = value;
				}
				MetaDataInName currentMeta;
				for(const auto & keyVal : valuePairs) {
					if(currentMeta.containsMeta(keyVal.first)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "error, already have field: " << keyVal.first << " for meta field: " << metaField <<", error in line: " << "\n";
						ss << line << "\n";
						throw std::runtime_error{ss.str()};
					}
					currentMeta.addMeta(keyVal.first, keyVal.second);
				}
				if(!currentMeta.containsMeta("ID")) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "all metafields should have at least ID field, error for meta: " << metaField << ", error on line: " << "\n";
					ss << line << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto checkForRequiredField = [&line,&currentMeta](const VecStr & requiredFields, const std::string & metafield, const std::string & funcName) {
					VecStr missingFields;
					for(const auto & f : requiredFields) {
						if(!currentMeta.containsMeta(f)) {
							missingFields.emplace_back(f);
						}
					}
					if(!missingFields.empty()) {
						std::stringstream ss;
						ss << funcName << ", error " << "error in processing metafield " << metafield << "was missing the following required fields:" << njh::conToStr(missingFields, ",") << ", only found the following: " << njh::conToStr(njh::getVecOfMapKeys(currentMeta.meta_), ",")  << " error for line:"<< "\n";
						ss << line << "\n";
						throw std::runtime_error{ss.str()};
					}
				};

				if(metaField == "INFO") {
					VecStr requiredMetaFields{"ID", "Number", "Type", "Description"};
					checkForRequiredField(requiredMetaFields, metaField, __PRETTY_FUNCTION__);
					InfoEntry info_entry(currentMeta.getMeta("ID"), currentMeta.getMeta("Number"), currentMeta.getMeta("Type"), currentMeta.getMeta("Description"));
					if(currentMeta.containsMeta("Source")) {
						info_entry.source_ = currentMeta.getMeta("Source");
					}
					if(currentMeta.containsMeta("Version")) {
						info_entry.version_ = currentMeta.getMeta("Version");
					}
					if(njh::in(info_entry.id_, ret.infoEntries_)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << " already have entry on INFO: " << info_entry.id_ << "\n";
						throw std::runtime_error{ss.str()};
					}
					ret.infoEntries_.emplace(info_entry.id_, info_entry);
				} else if (metaField == "FILTER") {
					VecStr requiredMetaFields{"ID", "Description"};
					checkForRequiredField(requiredMetaFields, metaField, __PRETTY_FUNCTION__);
					ret.filterEntries_.emplace_back(currentMeta.getMeta("ID"), currentMeta.getMeta("Description"));
				} else if (metaField == "FORMAT") {
					VecStr requiredMetaFields{"ID", "Number", "Type", "Description"};
					checkForRequiredField(requiredMetaFields, metaField, __PRETTY_FUNCTION__);
					FormatEntry format_entry(currentMeta.getMeta("ID"), currentMeta.getMeta("Number"), currentMeta.getMeta("Type"), currentMeta.getMeta("Description"));
					if(njh::in(format_entry.id_, ret.formatEntries_)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << " already have entry on FORMAT: " << format_entry.id_ << "\n";
						throw std::runtime_error{ss.str()};
					}
					ret.formatEntries_.emplace(format_entry.id_, format_entry);
				} else if (metaField == "contig") {
					VecStr requiredMetaFields{"ID", "length"};
					checkForRequiredField(requiredMetaFields, metaField, __PRETTY_FUNCTION__);
					ContigEntry contig_entry(currentMeta.getMeta("ID"), currentMeta.getMeta<uint32_t>("length"));
					if(currentMeta.containsMeta("assembly")) {
						contig_entry.assembly_ = currentMeta.getMeta("assembly");
					}
					if(currentMeta.containsMeta("md5")) {
						contig_entry.md5_ = currentMeta.getMeta("md5");
					}
					if(currentMeta.containsMeta("species")) {
						contig_entry.species_ = currentMeta.getMeta("species");
					}
					VecStr knownFields{"ID", "length", "assembly", "md5", "species"};
					for(const auto & otherMetas : currentMeta.meta_) {
						if(!njh::in(otherMetas.first, knownFields)) {
							if(njh::in(otherMetas.first, contig_entry.otherKeysValues_)) {
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "already have meta for contig meta field, " << otherMetas.first << ", error on line:" << "\n";
								ss << line << "\n";
								throw std::runtime_error{ss.str()};
							}
							contig_entry.otherKeysValues_[otherMetas.first] = otherMetas.second;
						}
					}
					if(njh::in(contig_entry.id_, ret.contigEntries_)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << " already have entry on contig: " << contig_entry.id_ << "\n";
						throw std::runtime_error{ss.str()};
					}
					ret.contigEntries_.emplace(contig_entry.id_, contig_entry);
				} else {
					ret.otherHeaderMetaFields_.emplace(metaField, currentMeta);
				}
			}
		}
	}


	return ret;
}



VCFOutput VCFOutput::comnbineVCFs(const std::vector<bfs::path> & vcfsFnps,
		const std::set<std::string> & sampleNamesSet,
		bool doNotRescueVariantCallsAccrossTargets){
	const auto& firstVcfFnp = vcfsFnps.front();
	auto firstVcf = VCFOutput::readInHeader(firstVcfFnp);
	{
		InputStream in(firstVcfFnp);
		firstVcf.addInRecordsFromFile(in);
	}
	//add in blanks for any missing samples for any of the records
	firstVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
	if(firstVcf.records_.empty()) {
		firstVcf.samples_ = VecStr{sampleNamesSet.begin(), sampleNamesSet.end()};
	}
	for(const auto & gvcfFnp : vcfsFnps) {
		if(firstVcfFnp == gvcfFnp) {
			continue;
		}
		auto currentGvcf = VCFOutput::readInHeader(gvcfFnp);

		{
			InputStream in(gvcfFnp);
			currentGvcf.addInRecordsFromFile(in);
			// firstVcf.addInRecordsFromFile(in);
		}
		//add in header
		//given we know what generated these vcf files there's not need to check the FORMAT, FILTER, or INFO tags

		//add in new contigs if any
		for(const auto & contig : currentGvcf.contigEntries_) {
			if(!njh::in(contig.first, firstVcf.contigEntries_)) {
				firstVcf.contigEntries_.emplace(contig);
			} else if(!(firstVcf.contigEntries_.at(contig.first) == contig.second)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "adding contig: " << contig.first << " but already have an entry for this but it doesn't match" << "\n";
				ss << "contig in master header: " << njh::json::writeAsOneLine(njh::json::toJson(firstVcf.contigEntries_.at(contig.first))) << "\n";
				ss << "adding contig: " << njh::json::writeAsOneLine(njh::json::toJson(contig.second)) << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		currentGvcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		firstVcf.records_.reserve(firstVcf.records_.size() + currentGvcf.records_.size());
		for(auto & record : currentGvcf.records_) {
			firstVcf.records_.emplace_back(std::move(record));
		}
		//std::move(currentPvcf.records_.begin(), currentPvcf.records_.back(). std::back_inserter(firstPVcf.records_));
	}
	firstVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
	firstVcf.sortRecords();
	//eliminate duplicate variant calls
	{

		if(firstVcf.records_.size() > 1) {
			std::set<uint32_t> positionsToErase;
			std::set<uint32_t> currentPositionsToComp;
			auto compSamePositions = [&currentPositionsToComp,&firstVcf,&positionsToErase,
				&doNotRescueVariantCallsAccrossTargets]() {
				if(currentPositionsToComp.size() > 1) {
					uint32_t bestPos = std::numeric_limits<uint32_t>::max();
					uint32_t bestNS = 0;
					for(const auto & checkPos : currentPositionsToComp) {
						auto currentNS = firstVcf.records_[checkPos].info_.getMeta<uint32_t>("NS");
						if(currentNS > bestNS) {
							bestNS = currentNS;
							bestPos = checkPos;
						}
					}
					for(const auto & checkPos : currentPositionsToComp) {
						if(bestPos != checkPos) {
							positionsToErase.emplace(checkPos);
						}
					}
					if(!doNotRescueVariantCallsAccrossTargets) {

						std::set<std::string> blankSamples;
						std::regex blankDataPattern("\\.(,\\.)*");
						for(const auto & sampInfo : firstVcf.records_[bestPos].sampleFormatInfos_) {
							//check to see if sample has all blank meta fields (e.g. . .,. .,.,. etc)
							if(std::all_of(sampInfo.second.meta_.begin(), sampInfo.second.meta_.end(),
								[&blankDataPattern](const auto & meta) {
									return std::regex_match(meta.second, blankDataPattern);
								})) {
								blankSamples.emplace(sampInfo.first);
							}
						}

						for(const auto & samp : blankSamples) {
							std::vector<uint32_t> rowPositionsWithNonBlankSamples;
							for(const auto & checkPos : currentPositionsToComp) {
								if (bestPos != checkPos) {
									//check to see if this blank sample has data for the overlapping variant call
									if (!std::all_of(firstVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.begin(),
									                firstVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.end(),
									                [&blankDataPattern](const auto& meta) {
										                return std::regex_match(meta.second, blankDataPattern);
									                }) &&
									                "." != firstVcf.records_[checkPos].sampleFormatInfos_[samp].getMeta("DP") &&
																	"0" != firstVcf.records_[checkPos].sampleFormatInfos_[samp].getMeta("DP")) {

										rowPositionsWithNonBlankSamples.emplace_back(checkPos);
									}
								}
							}

							if(!rowPositionsWithNonBlankSamples.empty()) {
								//replace with the row with the best read depth
								uint32_t bestReadDepth = 0;
								uint32_t bestRowPositionWithData = std::numeric_limits<uint32_t>::max();
								for(const auto currentPos : rowPositionsWithNonBlankSamples) {

									auto readDepth = firstVcf.records_[currentPos].sampleFormatInfos_[samp].getMeta<uint32_t>("DP");
									if(readDepth > bestReadDepth) {
										bestReadDepth = readDepth;
										bestRowPositionWithData = currentPos;
									}
								}
								firstVcf.records_[bestPos].sampleFormatInfos_[samp] = firstVcf.records_[bestRowPositionWithData].sampleFormatInfos_[samp];
								//add to bestRefPos the NS, AN, AC
								//// NS will increase by 1 (though have to investigate that it's possible a sample is blank from the first initial call becauase it might already be counted)
								firstVcf.records_[bestPos].info_.addMeta("NS",firstVcf.records_[bestPos].info_.getMeta<uint32_t>("NS") + 1,true);
																	//have to check to see if the same variants are present
								if(firstVcf.records_[bestRowPositionWithData].alts_ != firstVcf.records_[bestPos].alts_) {
									VecStr altsNotInBestPos;
									for(const auto & alt : firstVcf.records_[bestRowPositionWithData].alts_) {
										if(njh::notIn(alt, firstVcf.records_[bestPos].alts_)) {
											altsNotInBestPos.emplace_back(alt);
										}
									}
									VecStr altsNotInBestPosWithData;
									for(const auto & alt : firstVcf.records_[bestPos].alts_) {
										if(njh::notIn(alt, firstVcf.records_[bestRowPositionWithData].alts_)) {
											altsNotInBestPosWithData.emplace_back(alt);
										}
									}
									if(!altsNotInBestPosWithData.empty() && altsNotInBestPos.empty()){
									   //replacement has to add the missing alt alleles
										auto AD = firstVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD");
										auto AD_toks = tokenizeString(AD, ",");
										auto AF = firstVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AF");
										auto AF_toks = tokenizeString(AF, ",");

										std::string AD_replacement = AD_toks.front();
										std::string AF_replacement = AF_toks.front();
										std::unordered_map<std::string, uint32_t> altsInBestPosWithDataKey;
										for(const auto & idx : iter::enumerate(firstVcf.records_[bestRowPositionWithData].alts_)) {
											altsInBestPosWithDataKey[idx.element] = idx.index;
										}
										for(const auto & alt : firstVcf.records_[bestPos].alts_) {
											if(!AD_replacement.empty()) {
												AD_replacement += ",";
												AF_replacement += ",";
											}
											if(njh::in(alt, altsInBestPosWithDataKey)) {
												AD_replacement += AD_toks[altsInBestPosWithDataKey[alt] + 1];
												AF_replacement += AF_toks[altsInBestPosWithDataKey[alt] + 1];
											} else {
												AD_replacement += "0";
												AF_replacement += "0";
											}
										}
										firstVcf.records_[bestPos].sampleFormatInfos_[samp].addMeta("AD", AD_replacement, true);
										firstVcf.records_[bestPos].sampleFormatInfos_[samp].addMeta("AF", AF_replacement, true);
									} else if(!altsNotInBestPos.empty()){
										//all other samples have to add the missing alts from the replacement
										//if the replacement is also missing the original alts this will handle that as well
										auto newAltsSet = njh::vecToSet(concatVecs(firstVcf.records_[bestRowPositionWithData].alts_, firstVcf.records_[bestPos].alts_));
										VecStr newAlts(newAltsSet.begin(), newAltsSet.end());
										njh::sort(newAlts);
										std::unordered_map<std::string, uint32_t> altsInBestPosKey;
										for(const auto & idx : iter::enumerate(firstVcf.records_[bestPos].alts_)) {
											altsInBestPosKey[idx.element] = idx.index;
										}
										for(const auto & currentSample : firstVcf.samples_) {
											if(samp == currentSample) {
												continue;
											}
											//replacement has to add the missing alt alleles
											auto AD = firstVcf.records_[bestPos].sampleFormatInfos_[currentSample].getMeta("AD");
											auto AD_toks = tokenizeString(AD, ",");
											auto AF = firstVcf.records_[bestPos].sampleFormatInfos_[currentSample].getMeta("AF");
											auto AF_toks = tokenizeString(AF, ",");

											std::string AD_replacement = AD_toks.front();
											std::string AF_replacement = AF_toks.front();
											for(const auto & alt : newAlts) {
												if(!AD_replacement.empty()) {
													AD_replacement += ",";
													AF_replacement += ",";
												}
												if(njh::in(alt, altsInBestPosKey)) {
													AD_replacement += AD_toks[altsInBestPosKey[alt] + 1];
													AF_replacement += AF_toks[altsInBestPosKey[alt] + 1];
												} else {
													AD_replacement += "0";
													AF_replacement += "0";
												}
											}
											firstVcf.records_[bestPos].sampleFormatInfos_[currentSample].addMeta("AD", AD_replacement, true);
											firstVcf.records_[bestPos].sampleFormatInfos_[currentSample].addMeta("AF", AF_replacement, true);
										}

										//
										if(!altsNotInBestPosWithData.empty()) {
											//replacement has to add the missing alt alleles
											auto AD = firstVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD");
											auto AD_toks = tokenizeString(AD, ",");
											auto AF = firstVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AF");
											auto AF_toks = tokenizeString(AF, ",");

											std::string AD_replacement = AD_toks.front();
											std::string AF_replacement = AF_toks.front();
											std::unordered_map<std::string, uint32_t> altsInBestPosWithDataKey;
											for(const auto & idx : iter::enumerate(firstVcf.records_[bestRowPositionWithData].alts_)) {
												altsInBestPosWithDataKey[idx.element] = idx.index;
											}
											for(const auto & alt : newAlts) {
												if(!AD_replacement.empty()) {
													AD_replacement += ",";
													AF_replacement += ",";
												}
												if(njh::in(alt, altsInBestPosWithDataKey)) {
													AD_replacement += AD_toks[altsInBestPosWithDataKey[alt] + 1];
													AF_replacement += AF_toks[altsInBestPosWithDataKey[alt] + 1];
												} else {
													AD_replacement += "0";
													AF_replacement += "0";
												}
											}
											firstVcf.records_[bestPos].sampleFormatInfos_[samp].addMeta("AD", AD_replacement, true);
											firstVcf.records_[bestPos].sampleFormatInfos_[samp].addMeta("AF", AF_replacement, true);
										}

										//replace current alts
										//replace AC and AF with zeros so that underneath they get properly modified
										auto AC = firstVcf.records_[bestPos].info_.getMeta("AC");
										auto AC_toks = tokenizeString(AC, ",");
										auto AF = firstVcf.records_[bestPos].info_.getMeta("AF");
										auto AF_toks = tokenizeString(AF, ",");
										auto SC = firstVcf.records_[bestPos].info_.getMeta("SC");
										auto SC_toks = tokenizeString(SC, ",");
										auto PREV = firstVcf.records_[bestPos].info_.getMeta("PREV");
										auto PREV_toks = tokenizeString(PREV, ",");
										std::string AC_replacement;
										std::string AF_replacement;
										std::string SC_replacement;
										std::string PREV_replacement;
										for(const auto & alt : newAlts) {
											if(!AC_replacement.empty()) {
												AC_replacement += ",";
												AF_replacement += ",";
												SC_replacement += ",";
												PREV_replacement += ",";
											}
											if(njh::in(alt, altsInBestPosKey)) {
												AC_replacement += AC_toks[altsInBestPosKey[alt]];
												AF_replacement += AF_toks[altsInBestPosKey[alt]];
												SC_replacement += SC_toks[altsInBestPosKey[alt]];
												PREV_replacement += PREV_toks[altsInBestPosKey[alt]];
											} else {
												AC_replacement += "0";
												AF_replacement += "0";
												SC_replacement += "0";
												PREV_replacement += "0";
											}
										}
										firstVcf.records_[bestPos].info_.addMeta("AC", AC_replacement, true);
										firstVcf.records_[bestPos].info_.addMeta("AF", AF_replacement, true);
										firstVcf.records_[bestPos].info_.addMeta("SC", SC_replacement, true);
										firstVcf.records_[bestPos].info_.addMeta("PREV", PREV_replacement, true);
										firstVcf.records_[bestPos].alts_ = newAlts;
									}
								}
								//// AN will increase by the number non-zero allele calls for the sample
								auto ADs = tokenizeString(firstVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD"), ",");
								auto ANcount = std::count_if(ADs.begin(), ADs.end(), [](const std::string & str){ return "0" != str;});
								firstVcf.records_[bestPos].info_.addMeta("AN",firstVcf.records_[bestPos].info_.getMeta<uint32_t>("AN") + ANcount,true);
								//// AC will increase the counts of allele that arent' 0 (and since we are adding samples with data there shouldn't be .)
								//// subsequently AF will also have to re-calculated
								////// will have to reconstruct the AC/AF for the best ref pos

								////// AC
								auto AC = firstVcf.records_[bestPos].info_.getMeta("AC");
								auto AC_toks = tokenizeString(AC, ",");
								//the ADs go ref, variant1, variant2 etc
								if(AC_toks.size() + 1 != ADs.size()) {
									std::stringstream ss;
									ss << __PRETTY_FUNCTION__ << ", error " << "AC_toks.size() + 1 should equal ADs.size()" << "\n";
									ss << "AC_toks.size() + 1: " << AC_toks.size() + 1 << ", " << "ADs.size(): " << ADs.size() << "\n";
									ss << "AC_toks: " << njh::conToStr(AC_toks, ",") << "; " << "ADs: " << njh::conToStr(ADs, ",") << "\n";
									ss << "firstVcf.records_[bestPos].alts_                : " << njh::conToStr(firstVcf.records_[bestPos].alts_, ",") << "\n";
									ss << "firstVcf.records_[bestRowPositionWithData].alts_: " << njh::conToStr(firstVcf.records_[bestRowPositionWithData].alts_, ",") << "\n";
									throw std::runtime_error{ss.str()};
								}
								//ADs_tok includes the reference as it's first tok so skip it and adjust index by 1
								for (const auto& ADs_tok: iter::enumerate(ADs)) {
									//skip reference
									if (ADs_tok.index != 0) {
										//if ADs_tok.element does not equal 0 then add 1
										if (ADs_tok.element != "0") {
											AC_toks[ADs_tok.index - 1] = estd::to_string(
												njh::StrToNumConverter::stoToNum<uint32_t>(AC_toks[ADs_tok.index - 1]) + 1);
										}
									}
								}
								//replace the updated AC
								firstVcf.records_[bestPos].info_.addMeta("AC", njh::conToStr(AC_toks, ","), true);
								//////// AF
								// re-calculate the AFs based on the new AC and AN
								double AN = firstVcf.records_[bestPos].info_.getMeta<uint32_t>("AN");
								VecStr AFs;
								for(const auto & AC_tok : AC_toks) {
									AFs.emplace_back(estd::to_string(njh::StrToNumConverter::stoToNum<uint32_t>(AC_tok)/AN));
								}
								//replace the updated AF
								firstVcf.records_[bestPos].info_.addMeta("AF", njh::conToStr(AFs, ","), true);


								////// SC
								auto SC = firstVcf.records_[bestPos].info_.getMeta("SC");
								auto SC_toks = tokenizeString(SC, ",");
								//the ADs go ref, variant1, variant2 etc
								if(SC_toks.size() + 1 != ADs.size()) {
									std::stringstream ss;
									ss << __PRETTY_FUNCTION__ << ", error " << "SC_toks.size() + 1 should equal ADs.size()" << "\n";
									ss << "SC_toks.size() + 1: " << SC_toks.size() + 1 << ", " << "ADs.size(): " << ADs.size() << "\n";
									ss << "SC_toks: " << njh::conToStr(SC_toks, ",") << "; " << "ADs: " << njh::conToStr(ADs, ",") << "\n";
									ss << "firstVcf.records_[bestPos].alts_                : " << njh::conToStr(firstVcf.records_[bestPos].alts_, ",") << "\n";
									ss << "firstVcf.records_[bestRowPositionWithData].alts_: " << njh::conToStr(firstVcf.records_[bestRowPositionWithData].alts_, ",") << "\n";
									throw std::runtime_error{ss.str()};
								}
								//ADs_tok includes the reference as it's first tok so skip it and adjust index by 1
								for (const auto& ADs_tok: iter::enumerate(ADs)) {
									//skip reference
									if (ADs_tok.index != 0) {
										//if ADs_tok.element does not equal 0 then add 1
										if (ADs_tok.element != "0") {
											SC_toks[ADs_tok.index - 1] = estd::to_string(
												njh::StrToNumConverter::stoToNum<uint32_t>(SC_toks[ADs_tok.index - 1]) + 1);
										}
									}
								}
								//////// PREV
								/// re-calculate PREV based on the new SC and NS
								double NS = firstVcf.records_[bestPos].info_.getMeta<uint32_t>("NS");
								VecStr PREVs;
								for(const auto & SC_tok : SC_toks) {
									PREVs.emplace_back(estd::to_string(njh::StrToNumConverter::stoToNum<uint32_t>(SC_tok)/NS));
								}
								//replace the updated PREV
								firstVcf.records_[bestPos].info_.addMeta("PREV", njh::conToStr(PREVs, ","), true);
							}
						}
					}
				}
				currentPositionsToComp.clear();
			};
			for(uint32_t pos = 1; pos < firstVcf.records_.size(); ++pos) {
				if(firstVcf.records_[pos].pos_ == firstVcf.records_[pos-1].pos_ &&
					firstVcf.records_[pos].ref_ == firstVcf.records_[pos-1].ref_) {
					currentPositionsToComp.emplace(pos);
					currentPositionsToComp.emplace(pos-1);
					} else {
						compSamePositions();
					}
			}
			compSamePositions();
			// std::cout << njh::conToStr(positionsToErase, ",") << std::endl;
			for(const auto posToErase : iter::reversed(positionsToErase)) {
				firstVcf.records_.erase(firstVcf.records_.begin() + posToErase);
			}
		}
	}
  return firstVcf;
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
}
} // namespace njhseq
