//
// Created by Nicholas Hathaway on 1/30/24.
//

#include "VCFOutput.hpp"

namespace njhseq {

namespace {
/**
 * \brief whether a header value is already wrapped in double quotes
 */
bool isQuoted(const std::string & val) {
	return val.size() >= 2 && '"' == val.front() && '"' == val.back();
}

/**
 * \brief generate the INFO column for a record in the header's INFO order, keys the record doesn't have are left off and "." is used when there are none
 */
std::string genInfoColumn(const std::map<std::string, VCFOutput::InfoEntry> & infoEntries, const MetaDataInName & info) {
	std::string infoOut;
	for (const auto & infoKey: infoEntries) {
		const auto & id = infoKey.second.id_;
		const auto search = info.meta_.find(id);
		if (info.meta_.end() == search) {
			continue;
		}
		if(!infoOut.empty()) {
			infoOut +=";";
		}
		if (infoKey.second.type_ == "Flag") {
			infoOut += id;
		} else {
			infoOut += id + "=" + search->second;
		}
	}
	if(infoOut.empty()) {
		infoOut = ".";
	}
	return infoOut;
}

/**
 * \brief the nth ':' separated value of a sample's values
 */
std::string_view nthSampleField(const std::string & values, const size_t n) {
	size_t start = 0;
	for (size_t i = 0; i < n; ++i) {
		start = values.find(':', start);
		if (std::string::npos == start) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "sample values: " << values << " don't have a value at position " << n << "\n";
			throw std::runtime_error{ss.str()};
		}
		++start;
	}
	const auto end = values.find(':', start);
	return std::string_view(values).substr(start, std::string::npos == end ? std::string::npos : end - start);
}
}  // namespace

bool VCFOutput::isCommaSepSingleChars(const std::string & val, const std::string & allowed) {
	if (val.empty() || 0 == val.size() % 2) {
		return false;
	}
	for (size_t pos = 0; pos < val.size(); ++pos) {
		if (0 == pos % 2) {
			if (std::string::npos == allowed.find(val[pos])) {
				return false;
			}
		} else if (',' != val[pos]) {
			return false;
		}
	}
	return true;
}

uint32_t VCFOutput::VCFRecord::getNumberOfAlleles() const {
	return 1 + alts_.size();
}


void VCFOutput::VCFRecord::autoAddTYPEField() {
	//going to need a better determination of TYPE to be able to determine complex types better
	std::string TYPE;

	for(const auto & alt : alts_) {
		if(!TYPE.empty()) {
			TYPE +=",";
		}
		if(ref_.size() == 1 && alt.size() == 1) {
			TYPE += "snp";
		} else if(ref_.size() > 1 && alt.size() == ref_.size()) {
			uint32_t snpCount = 0;
			for(const auto pos : iter::range(ref_.size())) {
				if(ref_[pos] != alt[pos]) {
						++snpCount;
				}
			}
			if(snpCount == 1) {
				TYPE += "snp";
			} else {
				TYPE += "mnp";
			}
			// //if all bases don't equal each other than it's a multiple snp TYPE, if some equal each other then it's a 'complex' TYPE
			// if(std::all_of(ref_.begin(), ref_.end(), [&alt](char c, size_t i = 0) mutable {
			// 				return c != alt[i++];
			// 		})) {
			// 	TYPE += "mnp";
			// } else {
			// 	TYPE += "complex";
			// }
		} else if(ref_.size() > alt.size()) {
			TYPE += "del";
		} else if (ref_.size() < alt.size()) {
			TYPE += "ins";
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " did not determine TYPE" << "\n";
			ss << "ref_: " << ref_ << "\n";
			ss << "alt: " << alt << "\n";
			throw std::runtime_error{ss.str()};
		}
	}

	info_.addMeta("TYPE", TYPE, true);

}


size_t VCFOutput::VCFRecord::formatKeyIndex(const std::string & key) const {
	const auto search = std::find(formatKeys_.begin(), formatKeys_.end(), key);
	return formatKeys_.end() == search ? std::string::npos : static_cast<size_t>(search - formatKeys_.begin());
}

std::vector<std::string_view> VCFOutput::VCFRecord::splitSampleValues(const std::string & values) {
	std::vector<std::string_view> ret;
	size_t start = 0;
	while (true) {
		const auto end = values.find(':', start);
		if (std::string::npos == end) {
			ret.emplace_back(values.data() + start, values.size() - start);
			break;
		}
		ret.emplace_back(values.data() + start, end - start);
		start = end + 1;
	}
	return ret;
}

std::string VCFOutput::VCFRecord::getSampleField(const size_t sampleIndex, const std::string & key) const {
	const auto keyIndex = formatKeyIndex(key);
	if (std::string::npos == keyIndex) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no FORMAT field " << key << ", current fields are: " << njh::conToStr(formatKeys_, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return std::string(nthSampleField(sampleValues_.at(sampleIndex), keyIndex));
}

MetaDataInName VCFOutput::VCFRecord::getSampleMeta(const size_t sampleIndex) const {
	MetaDataInName ret;
	if (formatKeys_.empty()) {
		return ret;
	}
	const auto vals = splitSampleValues(sampleValues_.at(sampleIndex));
	if (vals.size() != formatKeys_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "sample values: " << sampleValues_.at(sampleIndex) << " should have " << formatKeys_.size() << " values for " << njh::conToStr(formatKeys_, ":") << "\n";
		throw std::runtime_error{ss.str()};
	}
	for (const auto & e : iter::enumerate(formatKeys_)) {
		ret.addMeta(e.element, std::string(vals[e.index]));
	}
	return ret;
}

void VCFOutput::VCFRecord::setSampleMeta(const size_t sampleIndex, const MetaDataInName & meta) {
	//keys new to the record are added in sorted order
	auto metaKeys = njh::getVecOfMapKeys(meta.meta_);
	njh::sort(metaKeys);
	for (const auto & key : metaKeys) {
		if (std::string::npos == formatKeyIndex(key)) {
			setFormatValues(key, VecStr(sampleValues_.size(), "."));
		}
	}
	if (sampleIndex >= sampleValues_.size()) {
		sampleValues_.resize(sampleIndex + 1);
	}
	std::string values;
	for (size_t keyIndex = 0; keyIndex < formatKeys_.size(); ++keyIndex) {
		if (keyIndex > 0) {
			values += ":";
		}
		const auto search = meta.meta_.find(formatKeys_[keyIndex]);
		values += meta.meta_.end() == search ? std::string(".") : search->second;
	}
	sampleValues_[sampleIndex] = std::move(values);
}

void VCFOutput::VCFRecord::setFormatValues(const std::string & key, const VecStr & values, const bool atFront) {
	if (values.size() != sampleValues_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "need a value for each of the " << sampleValues_.size() << " samples, not " << values.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	const auto keyIndex = formatKeyIndex(key);
	if (std::string::npos == keyIndex) {
		const bool hadKeys = !formatKeys_.empty();
		if (atFront) {
			formatKeys_.insert(formatKeys_.begin(), key);
		} else {
			formatKeys_.emplace_back(key);
		}
		for (size_t sampleIndex = 0; sampleIndex < sampleValues_.size(); ++sampleIndex) {
			auto & sampleValues = sampleValues_[sampleIndex];
			if (!hadKeys) {
				sampleValues = values[sampleIndex];
			} else if (atFront) {
				sampleValues.insert(0, values[sampleIndex] + ":");
			} else {
				sampleValues.append(":" + values[sampleIndex]);
			}
		}
	} else {
		for (size_t sampleIndex = 0; sampleIndex < sampleValues_.size(); ++sampleIndex) {
			auto & sampleValues = sampleValues_[sampleIndex];
			const auto field = nthSampleField(sampleValues, keyIndex);
			sampleValues.replace(field.data() - sampleValues.data(), field.size(), values[sampleIndex]);
		}
	}
}

void VCFOutput::VCFRecord::addFieldDefaultValue(const std::string & field, const std::string & defaultValue, const bool replace) {
	VecStr values;
	values.reserve(sampleValues_.size());
	for (const auto & sampleValues : sampleValues_) {
		if (!replace && std::string::npos != formatKeyIndex(field)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "attempting to add FORMAT field " << field << " that's already present, use replace = true to replace" << "\n";
			throw std::runtime_error{ss.str()};
		}
		//check if all other fields are empty, then use an empty value
		bool allBlanks = true;
		if (!formatKeys_.empty()) {
			for (const auto & val : splitSampleValues(sampleValues)) {
				if (!isCommaSepSingleChars(std::string(val), ".")) {
					allBlanks = false;
					break;
				}
			}
		}
		values.emplace_back(allBlanks ? std::string(".") : defaultValue);
	}
	setFormatValues(field, values);
}

void VCFOutput::VCFRecord::autoAdd_AN_AC_AF_InfoFields() {
	std::unordered_map<std::string, uint32_t> alleleCounts;
	uint32_t AN = 0;
	const auto GTIndex = formatKeyIndex("GT");
	for (const auto & sampleValues : sampleValues_) {
		if (std::string::npos == GTIndex) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field GT in field if going to add total AN, AC, and AF fields " << "\n";
			ss << "Current fields are : " << njh::conToStr(formatKeys_, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		const auto GT = nthSampleField(sampleValues, GTIndex);
		if ("." != GT) {
			auto GT_toks = tokenizeString(std::string(GT), "/");
			for (const auto& tok: GT_toks) {
				++alleleCounts[tok];
				++AN;
			}
		}
	}
	auto alleles = getVectorOfMapKeys(alleleCounts);
	std::vector<uint32_t> allelesNumeric;
	allelesNumeric.reserve(alleles.size());
	for (auto& allele: alleles) {
		allelesNumeric.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(allele));
	}
	if (vectorMaximum(allelesNumeric) > alts_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " genotype number: " <<
				vectorMaximum(allelesNumeric) << " can't be more than the alts_.size(): " << alts_.size() << " + 1 " << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<double> afs;
	std::string ACs;
	//the allele counts include the reference so alts.size() will be equal to the genotype number, e.g. for alts of size of 2, there will be 0, 1, 2 (will skip 0)
	for (const auto pos: iter::range(1UL, alts_.size() + 1)) {
		if (!ACs.empty()) {
			ACs += ",";
		}
		//if missing from the counts will default to zero
		ACs += estd::to_string(alleleCounts[estd::to_string(pos)]);
		//no called genotypes gives an AF of 0
		afs.emplace_back(0 == AN ? 0.0 : alleleCounts[estd::to_string(pos)] / static_cast<double>(AN));
	}
	info_.addMeta("AN", AN, true);
	info_.addMeta("AC", ACs, true);
	info_.addMeta("AF", njh::conToStr(afs, ","), true);
}


void VCFOutput::VCFRecord::autoAddWeightedAFRealField(const std::vector<size_t> & sampleOrder) {
	std::vector<double> weighted_real_afs(alts_.size(), 0); //total read per each alt observation
	double reference_sum = 0;
	const auto DPIndex = formatKeyIndex("DP");
	const auto ADIndex = formatKeyIndex("AD");
	auto addSample = [&](const std::string & sampleValues) {
		if (std::string::npos == DPIndex || std::string::npos == ADIndex) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field DP and AD in field if going to add weighted AF_REAL " << "\n";
			ss << "Current fields are : " << njh::conToStr(formatKeys_, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		const auto DPVal = nthSampleField(sampleValues, DPIndex);
		if("." != DPVal) {
			double DP = njh::StrToNumConverter::stoToNum<uint32_t>(std::string(DPVal));
			if(0 == DP) {
				return;
			}
			auto AD_toks = tokenizeString(std::string(nthSampleField(sampleValues, ADIndex)), ",");
			if(AD_toks.size() != 1 + alts_.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "AD_toks.size(): " << AD_toks.size() << " does not equal the alts_ and ref size: " << 1 + alts_.size() << "\n";
				throw std::runtime_error{ss.str()};
			}
			for(const auto & e : iter::enumerate(AD_toks)) {
				if (0 != e.index) {
					weighted_real_afs[e.index - 1] += njh::StrToNumConverter::stoToNum<uint32_t>(AD_toks[e.index])/DP;
				} else {
					reference_sum += njh::StrToNumConverter::stoToNum<uint32_t>(AD_toks[e.index])/DP;
				}
			}
		}
	};
	if (sampleOrder.empty()) {
		for (const auto & sampleValues : sampleValues_) {
			addSample(sampleValues);
		}
	} else {
		for (const auto sampleIndex : sampleOrder) {
			addSample(sampleValues_.at(sampleIndex));
		}
	}
	auto weighted_real_afs_sum = vectorSum(weighted_real_afs) + reference_sum;
	//no samples with depth leaves the AFs at 0
	if(weighted_real_afs_sum > 0) {
		for (auto & af : weighted_real_afs) {
			af /= weighted_real_afs_sum;
		}
	}
	info_.addMeta("AF_REAL", njh::conToStr(weighted_real_afs, ","), true);
}

void VCFOutput::VCFRecord::autoAddTotalDP_RO_AO_InfoFields() {
	uint32_t totalDP = 0;
	uint32_t totalRO = 0; //total reads for reference
	std::vector<uint32_t> totalAOs(alts_.size(), 0); //total read per each alt observation
	const auto DPIndex = formatKeyIndex("DP");
	const auto ADIndex = formatKeyIndex("AD");
	for (const auto & sampleValues : sampleValues_) {
		if (std::string::npos == DPIndex || std::string::npos == ADIndex) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field DP and AD in field if going to add total DP, RO, and AO fields " << "\n";
			ss << "Current fields are : " << njh::conToStr(formatKeys_, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		const auto DPVal = nthSampleField(sampleValues, DPIndex);
		if("." != DPVal) {
			totalDP += njh::StrToNumConverter::stoToNum<uint32_t>(std::string(DPVal));
			auto AD_toks = tokenizeString(std::string(nthSampleField(sampleValues, ADIndex)), ",");
			if(AD_toks.size() != 1 + alts_.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "AD_toks.size(): " << AD_toks.size() << " does not equal the alts_ and ref size: " << 1 + alts_.size() << "\n";
				throw std::runtime_error{ss.str()};
			}
			totalRO += njh::StrToNumConverter::stoToNum<uint32_t>(AD_toks[0]);
			for(const auto & e : iter::enumerate(AD_toks)) {
				if(e.index != 0) {
					totalAOs[e.index - 1] += njh::StrToNumConverter::stoToNum<uint32_t>(AD_toks[e.index]);
				}
			}
		}
	}

	info_.addMeta("DP", totalDP, true);
	info_.addMeta("RO", totalRO, true);
	info_.addMeta("AO", njh::conToStr(totalAOs, ","), true);

}



void VCFOutput::VCFRecord::addGTField(uint32_t ploidy) {
	const auto ADIndex = formatKeyIndex("AD");
	VecStr GTs;
	GTs.reserve(sampleValues_.size());
	for (const auto & sampleValues : sampleValues_) {
		if (std::string::npos == ADIndex) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field AD in field if going to add GT field " << "\n";
			ss << "Current fields are : " << njh::conToStr(formatKeys_, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		const std::string raw_AD(nthSampleField(sampleValues, ADIndex));
		//all depths are missing or zero
		if (isCommaSepSingleChars(raw_AD, ".0")) {
			GTs.emplace_back(".");
		} else {
			std::vector<uint32_t> depths(1 + alts_.size());
			auto ADs = tokenizeString(raw_AD, ",");
			std::transform(ADs.begin(), ADs.end(), depths.begin(), [](const auto &p) {
				return njh::StrToNumConverter::stoToNum<uint32_t>(p);
			});
			std::vector<uint32_t> rank(ADs.size());
			njh::iota(rank, 0U);
			njh::sort(rank, [&depths](const auto &pos1, const auto &pos2) {
				return depths[pos1] > depths[pos2];
			});

			//count alleles with depth
			uint32_t allelesWithDepth = std::count_if(depths.begin(), depths.end(), [](const auto d) { return d > 0; });

			std::vector<uint32_t> gts;
			if (ploidy <= allelesWithDepth) {
				for (uint32_t pos = 0; pos < ploidy; ++pos) {
					gts.emplace_back(rank[pos]);
				}
			} else {
				//fill what you can at first
				for (const auto pos: iter::range(rank.size())) {
					if(depths[rank[pos]] > 0) {
						gts.emplace_back(rank[pos]);
					}
				}
				auto sum = vectorSum(depths);
				std::vector<double> freqs(ADs.size());
				std::transform(depths.begin(), depths.end(), freqs.begin(), [&sum](const auto &p) { return p / sum; });
				uint32_t diff = ploidy - allelesWithDepth;
				for (const auto pos: iter::range(rank.size())) {
					auto add = std::min(diff, static_cast<uint32_t>(std::round(freqs[rank[pos]] * diff)));
					addOtherVec(gts, std::vector(add, rank[pos]));
					diff -= add;
					if (diff == 0) {
						break;
					}
				}
				if (diff != 0) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
							" diff should be zero, not:  " << diff << "\n";
					throw std::runtime_error{ss.str()};
				}
			}

			njh::sort(gts);
			std::string GT;
			for (const auto gt : gts) {
				if (!GT.empty()) {
					GT += "/";
				}
				GT += std::to_string(gt);
			}
			GTs.emplace_back(GT);
		}
	}
	if (!sampleValues_.empty()) {
		setFormatValues("GT", GTs, true);
	}
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
	uint32_t end = start + ref_.size();
	//a deletion keeps the anchor base and the end of the ref
	const bool pureDeletion = ref_.size() > 1 && std::all_of(alts_.begin(), alts_.end(), [this](const std::string& alt) {
		return !alt.empty() && alt.size() < ref_.size() && alt.front() == ref_.front() && njh::endsWith(ref_, alt.substr(1));
	});
	if (ref_.size() == 1 && !std::all_of(alts_.begin(), alts_.end(), [](const std::string& alt) {
		return alt.size() == 1;
	})) {
		//insertion, will give the region right before and right after the insertion
		end += 1;
	} else if (pureDeletion) {
		//only the deleted bases, skipping the anchor base
		start += 1;
	}
	//anything else (snps, mnps, complex records) is the full ref span
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
	ret["formatKeys_"] = njh::json::toJson(formatKeys_);
	ret["sampleValues_"] = njh::json::toJson(sampleValues_);
	return ret;
}

void VCFOutput::allAddGTFields(uint32_t ploidy) {
	if (njh::notIn("GT", formatEntries_)) {
		formatEntries_.emplace("GT", FormatEntry(
			                       "GT", "1", "String",
			                       "Genotype"
		                       ));
	}
	njh::for_each(
		records_, [&ploidy](auto &rec) { rec.addGTField(ploidy); });
}


void VCFOutput::allAutoAddDPFields() {
	if(njh::notIn("DP", infoEntries_)) {
		infoEntries_["DP"] =
			InfoEntry("DP", "1", "Integer", "Total read depth at the locus");
	}
	if(njh::notIn("RO", infoEntries_)) {
		infoEntries_["RO"] =
			InfoEntry("RO", "1", "Integer", "Read Count of full observations of the reference haplotype.");
	}
	if(njh::notIn("AO", infoEntries_)) {
		infoEntries_["AO"] =
			InfoEntry("AO", "A", "Integer", "Read Count of full observations of this alternate haplotype.");
	}
	njh::for_each(
	records_, [](auto &rec) { rec.autoAddTotalDP_RO_AO_InfoFields(); });
}



void VCFOutput::allAutoAdd_AN_AC_AF_InfoFields() {
	if (njh::notIn("AN", infoEntries_)) {
		infoEntries_.emplace("AN", InfoEntry(
			                     "AN", "1", "Integer",
			                     "Total Allele Depth dependent on ploidy, sum of AC with rest of depth being ref")
		);
	}
	if (njh::notIn("AC", infoEntries_)) {
		infoEntries_.emplace("AC", InfoEntry(
			                     "AC", "A", "Integer", "Allele Count dependent on set ploidy"
		                     ));
	}
	if (njh::notIn("AF", infoEntries_)) {
		infoEntries_.emplace("AF", InfoEntry(
			                     "AF", "A", "Float", "Allele Frequency dependent on ploidy, calulated AC/AN"
		                     ));
	}

	njh::for_each(
		records_, [](auto& rec) { rec.autoAdd_AN_AC_AF_InfoFields(); });

}

void VCFOutput::allAutoAddWeightedAFRealField() {
	if (njh::notIn("AF_REAL", infoEntries_)) {
		infoEntries_.emplace("AF_REAL", InfoEntry(
													 "AF_REAL", "A", "Float", "Allele Frequency not dependent on ploidy, calculated as AC/AN weighted by within sample frequencies"
												 ));
	}
	//the frequencies are summed in sample name order, only needed when samples_ isn't sorted
	std::vector<size_t> sampleOrder;
	if (!std::is_sorted(samples_.begin(), samples_.end())) {
		sampleOrder.resize(samples_.size());
		njh::iota(sampleOrder, static_cast<size_t>(0));
		njh::sort(sampleOrder, [this](const size_t s1, const size_t s2) { return samples_[s1] < samples_[s2]; });
	}
	njh::for_each(
		records_, [&sampleOrder](auto& rec) { rec.autoAddWeightedAFRealField(sampleOrder); });

}






void VCFOutput::allAutoAddTYPEFields() {
	if(njh::notIn("TYPE", infoEntries_)) {
		infoEntries_["TYPE"] =
			InfoEntry("TYPE", "A", "String", "The type of allele, either snp, mnp, ins, del, or complex");
	}
	njh::for_each(
records_, [](auto &rec) { rec.autoAddTYPEField(); });
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
	for (const auto& contigKey: contigEntries_) {
		const auto& contig = contigKey.second;
		vcfOut << "##contig=<"
				<< "ID=" << contig.id_ << ","
				<< "length=" << contig.length_;
		if (!contig.md5_.empty()) {
			vcfOut << "," << "md5=" << contig.md5_;
		}
		if(!contig.assembly_.empty()) {
			if(isQuoted(contig.assembly_)) {
				vcfOut << "," << "assembly=" << contig.assembly_ << "";
			}else {
				vcfOut << "," << "assembly=\"" << contig.assembly_ << "\"";
			}
		}
		if(!contig.species_.empty()) {
			if(isQuoted(contig.species_)) {
				vcfOut << "," << "species=" << contig.species_ << "";
			} else {
				vcfOut << "," << "species=\"" << contig.species_ << "\"";
			}
		}
		for(const auto & others : contig.otherKeysValues_) {
			if((std::string::npos != others.second.find(',') || njh::strHasWhitesapce(others.second)) && !(isQuoted(others.second)) ) {
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
		if(isQuoted(info.description_)) {
			vcfOut << "," << "Description=" << info.description_ << "";
		} else {
			vcfOut << "," << "Description=\"" << info.description_ << "\"";
		}
		if(!info.source_.empty()) {
			if(isQuoted(info.source_)) {
				vcfOut << "," << "Source=" << info.source_ << "";
			} else {
				vcfOut << "," << "Source=\"" << info.source_ << "\"";
			}
		}
		if(!info.version_.empty()) {
			if(isQuoted(info.version_)) {
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
		if (isQuoted(info.description_)) {
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
			if((std::string::npos != valuePairs.second.find(',') || njh::strHasWhitesapce(valuePairs.second)) && !(isQuoted(valuePairs.second)) ) {
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

void VCFOutput::changeContigNames(const std::unordered_map<std::string, std::string> & name_key) {
	VecStr missing_name;
	std::map<std::string, ContigEntry> replacement_map;
	for (auto & contig : contigEntries_) {
		if (njh::notIn(contig.second.id_, name_key)) {
			missing_name.emplace_back(contig.first);
		} else {
			auto replacement_name =  njh::mapAt(name_key, contig.second.id_);
			auto replacement_contig = contig.second;
			replacement_contig.id_ =  replacement_name;
			if (njh::in(replacement_name, replacement_map)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have replacement name: " << replacement_name << "\n";
				throw std::runtime_error{ss.str()};
			}
			replacement_map.emplace(replacement_name, replacement_contig);
		}
	}
	if (!missing_name.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "missing the following contig names from renmaing key: " << njh::conToStr(missing_name, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
	contigEntries_ = replacement_map;
	//rename in records
	for (auto & record: records_) {
		record.chrom_ = name_key.at(record.chrom_);
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
			const auto recRegion = rec.genRegion();
			for(const auto & reg : selectRegions) {
				if(reg.overlaps(recRegion)) {
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
		<< "\t" << (rec.qual_ == std::numeric_limits<uint32_t>::max() ? "." : estd::to_string(rec.qual_))
		<< "\t" << rec.filter_;
		vcfOut << "\t" << genInfoColumn(infoEntries_, rec.info_);
		vcfOut << std::endl;
	}
}


void VCFOutput::writeOutFixedAndSampleMeta(std::ostream& vcfOut, const std::vector<GenomicRegion>& selectRegions) const {
	vcfOut << "##fileformat=" << vcfFormatVersion_ << std::endl;
	writeOutHeaderFieldsOtherThanFormat(vcfOut);
	//write out formats
	VecStr formatOutputOrder;
	//force GT to be first field
	if(njh::in(std::string("GT"), formatEntries_)) {
		formatOutputOrder.emplace_back("GT");
	}
	for (const auto & infoKey: formatEntries_) {
		if(infoKey.first != "GT") {
			formatOutputOrder.emplace_back(infoKey.first);
		}
	}
	for (const auto & infoKey: formatOutputOrder) {
		const auto & info = formatEntries_.at(infoKey);
		vcfOut <<"##FORMAT=<"
		<< "ID=" << info.id_
		<< ","<< "Number=" << info.number_
		<< "," << "Type=" << info.type_;
		if(isQuoted(info.description_)) {
			vcfOut << "," << "Description=" << info.description_ << "";
		} else {
			vcfOut << "," << "Description=\"" << info.description_ << "\"";
		}
		vcfOut << ">" << std::endl;
	}
	//check samples

	//vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	vcfOut << njh::conToStr(headerNonSampleFields_, "\t") << "\t" << njh::conToStr(samples_, "\t") << std::endl;
	for(const auto & rec : records_) {
		if(rec.sampleValues_.size() != samples_.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "samples in records don't match the header samples" << "\n";
			ss << "record " << rec.chrom_ << " " << rec.pos_ << " " << rec.ref_ << " has " << rec.sampleValues_.size() << " samples, the header has " << samples_.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	std::string formatOut;
	for (const auto & infoKey: formatOutputOrder) {
		const auto & format = formatEntries_.at(infoKey);
		if(!formatOut.empty()) {
			formatOut +=":";
		}
		formatOut += format.id_;
	}
	for(const auto & rec : records_) {
		if(!selectRegions.empty()) {
			bool overlaps = false;
			const auto recRegion = rec.genRegion();
			for(const auto & reg : selectRegions) {
				if(reg.overlaps(recRegion)) {
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
		<< "\t" << (rec.qual_ == std::numeric_limits<uint32_t>::max() ? "." : estd::to_string(rec.qual_))
		<< "\t" << rec.filter_;
		vcfOut << "\t" << genInfoColumn(infoEntries_, rec.info_);
		vcfOut << "\t" << formatOut;
		if(rec.formatKeys_ == formatOutputOrder) {
			for(const auto & sampleValues : rec.sampleValues_) {
				vcfOut << "\t" << sampleValues;
			}
		} else {
			//put the values in the header's FORMAT order, a FORMAT field the record doesn't have is written as missing, e.g. a FORMAT field only some of the combined files had
			std::vector<size_t> keyIndexes;
			keyIndexes.reserve(formatOutputOrder.size());
			for (const auto & key : formatOutputOrder) {
				keyIndexes.emplace_back(rec.formatKeyIndex(key));
			}
			for(const auto & sampleValues : rec.sampleValues_) {
				const auto vals = rec.formatKeys_.empty() ? std::vector<std::string_view>{} : VCFRecord::splitSampleValues(sampleValues);
				vcfOut << "\t";
				for (const auto & e : iter::enumerate(keyIndexes)) {
					if (0 != e.index) {
						vcfOut << ":";
					}
					if (std::string::npos == e.element || e.element >= vals.size()) {
						vcfOut << ".";
					} else {
						vcfOut << vals[e.element];
					}
				}
			}
		}
		vcfOut << std::endl;
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

void VCFOutput::addInBlanksForAnyMissingSamples(const std::set<std::string>& samples) {
	//the final samples are the header samples and the requested samples, sorted
	std::set<std::string> allSamples(samples.begin(), samples.end());
	njh::addVecToSet(samples_, allSamples);
	alignSamples(VecStr(allSamples.begin(), allSamples.end()));
}

std::unordered_map<std::string, size_t> VCFOutput::getSampleIndexes() const {
	std::unordered_map<std::string, size_t> ret;
	ret.reserve(samples_.size());
	for (const auto & e : iter::enumerate(samples_)) {
		ret.emplace(e.element, e.index);
	}
	return ret;
}

void VCFOutput::alignSamples(const VecStr & newSamples) {
	std::unordered_map<std::string, size_t> newIndexes;
	newIndexes.reserve(newSamples.size());
	for (const auto & e : iter::enumerate(newSamples)) {
		newIndexes.emplace(e.element, e.index);
	}
	//where each current sample goes
	std::vector<size_t> toNewIndex;
	toNewIndex.reserve(samples_.size());
	for (const auto & samp : samples_) {
		const auto search = newIndexes.find(samp);
		if (newIndexes.end() == search) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "sample " << samp << " is missing from the new samples" << "\n";
			throw std::runtime_error{ss.str()};
		}
		toNewIndex.emplace_back(search->second);
	}
	const bool sameSamples = newSamples == samples_;
	for (auto & record : records_) {
		if (record.sampleValues_.empty() && !samples_.empty() && !record.formatKeys_.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "record " << record.chrom_ << " " << record.pos_ << " " << record.ref_ << " has FORMAT keys but no sample values" << "\n";
			throw std::runtime_error{ss.str()};
		}
		const bool hasSamples = !record.sampleValues_.empty();
		if (hasSamples && record.sampleValues_.size() != samples_.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "record " << record.chrom_ << " " << record.pos_ << " " << record.ref_ << " has " << record.sampleValues_.size() << " samples but there are " << samples_.size() << " samples" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if (sameSamples && hasSamples) {
			continue;
		}
		//blank samples get a value for every FORMAT entry, so the record gets any FORMAT entry it doesn't have, with missing values for its current samples
		for (const auto & format : formatEntries_) {
			if (std::string::npos == record.formatKeyIndex(format.first)) {
				record.setFormatValues(format.first, VecStr(record.sampleValues_.size(), "."));
			}
		}
		//the blank entry only depends on the record's number of alleles so it's made once per record
		std::string blank;
		for (const auto & e : iter::enumerate(record.formatKeys_)) {
			if (0 != e.index) {
				blank += ":";
			}
			uint32_t numberOfVals = 1;
			const auto format = formatEntries_.find(e.element);
			if (formatEntries_.end() != format) {
				if(njh::strAllDigits(format->second.number_)) {
					numberOfVals = njh::StrToNumConverter::stoToNum<uint32_t>(format->second.number_);
				}else if(format->second.number_ == "R") {
					numberOfVals = record.getNumberOfAlleles();
				}else if(format->second.number_ == "A") {
					numberOfVals = record.alts_.size();
				}
			}
			blank += njh::conToStr(VecStr{numberOfVals, "."}, ",");
		}
		VecStr aligned(newSamples.size(), blank);
		if (hasSamples) {
			for (const auto & e : iter::enumerate(toNewIndex)) {
				aligned[e.element] = std::move(record.sampleValues_[e.index]);
			}
		}
		record.sampleValues_ = std::move(aligned);
	}
	samples_ = newSamples;
}


VCFOutput::VCFRecord VCFOutput::processRecordLineForFixedData(const std::string & line) const {
	return processRecordLineForFixedData(line, tokenizeString(line, "\t"));
}

VCFOutput::VCFRecord VCFOutput::processRecordLineForFixedData(const std::string & line, const VecStr & toks) const {
	VCFRecord rec;
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
	rec.qual_ = toks[5] == "." ? std::numeric_limits<uint32_t>::max() : njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);
	rec.filter_ = toks[6];

	//info field
	VecStr warningsInfoField;

	// if(std::string::npos != toks[7].find(':')) {
	// 	warningsInfoField.emplace_back("info field can't have :");
	// }
	if(njh::strHasWhitesapce( toks[7])) {
		warningsInfoField.emplace_back("info field can't whitespace");
	}
	if(!warningsInfoField.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "encountered the following errors when processing info field for line: " << "\n";
		ss << "errors: " << njh::conToStr(warningsInfoField, ",") << "\n";
		ss << line << "\n";
		throw std::runtime_error{ss.str()};
	}

	//"." is an empty INFO column
	const auto infoToks = "." == toks[7] ? VecStr{} : tokenizeString(toks[7], ";");
	for(const auto & infoTok : infoToks) {
		std::string key;
		std::string val;
		uint32_t valCount = 0;
		if (njh::in(infoTok, infoEntries_) && infoEntries_.at(infoTok).number_ == "0" && infoEntries_.at(infoTok).type_ == "Flag") {
			key = infoTok;
			val = "true";
		} else {
			const auto equalSignPos = infoTok.find("=");
			if(equalSignPos == std::string::npos || equalSignPos == 0 || equalSignPos +1 >= infoTok.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "info toks should have an equal sign separating values" << "\n";
				ss << "infoTok: " << infoTok << "\n";
				ss << "info field: " << toks[7] << "\n";
				throw std::runtime_error{ss.str()};
			}
			key = infoTok.substr(0, equalSignPos);
			val = infoTok.substr(equalSignPos + 1);
			valCount = 1 + countOccurences(val, ",");
		}

		if(!njh::in(key, infoEntries_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no info entry to define " << key  << " options are " << njh::conToStr(njh::getVecOfMapKeys(infoEntries_)) << "\n";
			throw std::runtime_error{ss.str()};
		}
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


	auto toks = tokenizeString(line, "\t");
	auto rec = processRecordLineForFixedData(line, toks);
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
			if (toks.size() - 9 != samples_.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "line has " << toks.size() - 9 << " sample columns but the header has " << samples_.size() << " samples" << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error{ss.str()};
			}
			if (njh::vecToSet(formatToks).size() != formatToks.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "FORMAT has a repeated key: " << toks[8] << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error{ss.str()};
			}
			//the Number for each FORMAT key in this line
			std::vector<const std::string *> formatNumbers;
			std::vector<bool> formatNumberIsDigits;
			formatNumbers.reserve(formatToks.size());
			formatNumberIsDigits.reserve(formatToks.size());
			for(const auto & f : formatToks) {
				formatNumbers.emplace_back(&formatEntries_.at(f).number_);
				formatNumberIsDigits.emplace_back(njh::strAllDigits(formatEntries_.at(f).number_));
			}
			//check the samples
			for(const auto pos : iter::range(9UL, toks.size())) {
				const auto sampleToks = VCFRecord::splitSampleValues(toks[pos]);
				if(sampleToks.size() != formatToks.size()) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "sample info size, " <<  sampleToks.size() << ", doesn't match the expected number of " << formatToks.size() << "\n";
					ss << "sample info: " << toks[pos] << "\n";
					throw std::runtime_error{ss.str()};
				}
				for(const auto & e : iter::enumerate(sampleToks)) {
					const auto valCount = 1 + static_cast<uint32_t>(std::count(e.element.begin(), e.element.end(), ','));
					const auto & key = formatToks[e.index];
					const auto & val = e.element;
					const auto & number = *formatNumbers[e.index];
					if(number == "A") {
						if(valCount != rec.alts_.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << rec.alts_.size() << " but has " << valCount << " instead " << "\n";
							ss << "key: " << key << "\n";
							ss << "val: " << val << "\n";
							ss << "line: " << line << "\n";
							throw std::runtime_error{ss.str()};
						}
					} else if (number == "R") {
						if(valCount != rec.alts_.size() + 1) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << rec.alts_.size() + 1 << " but has " << valCount << " instead " << "\n";
							ss << "key: " << key << "\n";
							ss << "val: " << val << "\n";
							ss << "line: " << line << "\n";
							throw std::runtime_error{ss.str()};
						}
					} else if (formatNumberIsDigits[e.index]) {
						if(std::to_string(valCount) != number) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "info entry " << key << " should have " << number << " but has " << valCount << " instead " << "\n";
							ss << "key: " << key << "\n";
							ss << "val: " << val << "\n";
							ss << "line: " << line << "\n";
							throw std::runtime_error{ss.str()};
						}
					}
				}
			}
			//the sample columns are kept as is
			rec.formatKeys_ = std::move(formatToks);
			rec.sampleValues_.reserve(toks.size() - 9);
			std::move(toks.begin() + 9, toks.end(), std::back_inserter(rec.sampleValues_));
		}
	}

	return rec;
}


void VCFOutput::addInRecordsFixedDataFromFile(std::istream & in) {
	std::string line;
	// uint32_t count = 0;
	while(njh::files::crossPlatGetline(in, line)) {
		if(!line.empty() && line.front() != '#') {
			// std::cout << count++ << std::endl;
			records_.emplace_back(processRecordLineForFixedData(line));
		}
	}
}

void VCFOutput::addInRecordsFromFile(std::istream & in) {
	std::string line;
	// uint32_t count = 0;
	while(njh::files::crossPlatGetline(in, line)) {
		if(!line.empty() && line.front() != '#') {
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
	auto version = firstLine.substr(std::string("##fileformat=").size());
	VCFOutput ret;
	ret.vcfFormatVersion_ = version;
	InputStream input(fnp);
	std::string line;
	bool foundHeaderLine = false;
	while(njh::files::crossPlatGetline(input, line)) {
		if(!njh::beginsWith(line, "##")) {
			if(njh::beginsWith(line, "#CHROM")) {
				foundHeaderLine = true;
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
					if(std::string::npos == equalPos) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "meta field entry " << tok << " has no =, values with commas need to be in quotes, error for line: " << "\n";
						ss << line << "\n";
						throw std::runtime_error{ss.str()};
					}
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
	if(!foundHeaderLine) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << fnp << " has no #CHROM header line after the ## meta lines" << "\n";
		throw std::runtime_error{ss.str()};
	}

	return ret;
}


VCFOutput VCFOutput::combineVCFs(const std::vector<bfs::path> &vcfsFnps,
																	const combineVCFsPars &pars) {
	std::set<std::string> sampleNamesSet;
	for (const auto & fnp : vcfsFnps) {
		auto vcfHeader = VCFOutput::readInHeader(fnp);
		njh::addVecToSet(vcfHeader.samples_, sampleNamesSet);
	}
	return combineVCFs(vcfsFnps, sampleNamesSet, pars);
}

namespace {

/**
 * \brief whether a VCF value is missing, e.g. "." or ".,.,."
 */
bool isMissingVCFValue(const std::string_view & val) {
	return !val.empty() && std::all_of(val.begin(), val.end(), [](const char c) { return '.' == c || ',' == c; });
}

/**
 * \brief whether a sample has read data for a record, i.e. has a DP that is set and non-zero
 */
bool sampleHasData(const VCFOutput::VCFRecord & rec, const size_t sample) {
	const auto DPIndex = rec.formatKeyIndex("DP");
	if (std::string::npos == DPIndex || sample >= rec.sampleValues_.size()) {
		return false;
	}
	const auto DP = nthSampleField(rec.sampleValues_[sample], DPIndex);
	return !isMissingVCFValue(DP) && "0" != DP;
}

uint32_t getSampleDP(const VCFOutput::VCFRecord & rec, const size_t sample) {
	return sampleHasData(rec, sample) ? njh::StrToNumConverter::stoToNum<uint32_t>(rec.getSampleField(sample, "DP")) : 0;
}

/**
 * \brief get the AD field (ref, alt1, alt2, etc.) as counts, missing values are treated as 0
 */
std::vector<uint32_t> getADCounts(const VCFOutput::VCFRecord & rec, const std::string & AD) {
	const auto toks = tokenizeString(AD, ",");
	if (toks.size() != rec.getNumberOfAlleles()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "AD: " << AD << " should have " << rec.getNumberOfAlleles()
			 << " values for " << rec.chrom_ << ":" << rec.pos_ << " REF: " << rec.ref_ << " ALT: " << njh::conToStr(rec.alts_, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<uint32_t> ret;
	ret.reserve(toks.size());
	for (const auto & tok : toks) {
		ret.emplace_back(isMissingVCFValue(tok) ? 0 : njh::StrToNumConverter::stoToNum<uint32_t>(tok));
	}
	return ret;
}

/**
 * \brief re-order a per allele value (Number=A or Number=R) from one set of alts to another, alts not in fromAlts get fill
 */
std::string remapAlleleValues(const std::string & value,
                              const VecStr & fromAlts,
                              const VecStr & toAlts,
                              const bool includesRef,
                              const std::string & fill) {
	const size_t offset = includesRef ? 1 : 0;
	if ("." == value) {
		return njh::conToStr(VecStr(toAlts.size() + offset, "."), ",");
	}
	const auto toks = tokenizeString(value, ",");
	if (toks.size() != fromAlts.size() + offset) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "value: " << value << " should have " << fromAlts.size() + offset << " values, not " << toks.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, size_t> fromIndex;
	for (const auto & alt : iter::enumerate(fromAlts)) {
		fromIndex[alt.element] = alt.index;
	}
	VecStr ret;
	ret.reserve(toAlts.size() + offset);
	if (includesRef) {
		ret.emplace_back(toks.front());
	}
	for (const auto & alt : toAlts) {
		const auto idx = fromIndex.find(alt);
		ret.emplace_back(fromIndex.end() == idx ? fill : toks[idx->second + offset]);
	}
	return njh::conToStr(ret, ",");
}

/**
 * \brief change a record's alts to newAlts (which must contain all current alts) and re-order every Number=A/R INFO and FORMAT field to match
 */
void setRecordAlts(VCFOutput::VCFRecord & rec, const VecStr & newAlts, const VCFOutput & header) {
	if (rec.alts_ == newAlts) {
		return;
	}
	auto isPerAllele = [](const std::string & number) { return "A" == number || "R" == number; };
	auto numericFill = [](const std::string & type) { return std::string("Integer" == type || "Float" == type ? "0" : "."); };
	std::string currentField;
	try {
		for (const auto & info : header.infoEntries_) {
			if (isPerAllele(info.second.number_) && rec.info_.containsMeta(info.first)) {
				currentField = "INFO/" + info.first;
				rec.info_.addMeta(info.first,
				                  remapAlleleValues(rec.info_.getMeta(info.first), rec.alts_, newAlts, "R" == info.second.number_, numericFill(info.second.type_)),
				                  true);
			}
		}
		//the per allele FORMAT fields the record has
		std::vector<std::pair<size_t, const VCFOutput::FormatEntry *>> perAlleleFormats;
		for (const auto & format : header.formatEntries_) {
			const auto keyIndex = rec.formatKeyIndex(format.first);
			if (isPerAllele(format.second.number_) && std::string::npos != keyIndex) {
				perAlleleFormats.emplace_back(keyIndex, &format.second);
			}
		}
		if (!perAlleleFormats.empty()) {
			for (size_t sample = 0; sample < rec.sampleValues_.size(); ++sample) {
				//samples without data get missing values for the new alleles rather than 0
				const bool hasData = sampleHasData(rec, sample);
				VecStr vals;
				for (const auto & val : VCFOutput::VCFRecord::splitSampleValues(rec.sampleValues_[sample])) {
					vals.emplace_back(val);
				}
				for (const auto & format : perAlleleFormats) {
					currentField = "FORMAT/" + format.second->id_ + " for sample " + (sample < header.samples_.size() ? header.samples_[sample] : estd::to_string(sample));
					vals.at(format.first) = remapAlleleValues(vals.at(format.first), rec.alts_, newAlts, "R" == format.second->number_,
					                                          hasData ? numericFill(format.second->type_) : std::string("."));
				}
				rec.sampleValues_[sample] = njh::conToStr(vals, ":");
			}
		}
	} catch (const std::exception & e) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "re-ordering " << currentField << " for " << rec.chrom_ << ":" << rec.pos_
			 << " REF: " << rec.ref_ << " ALT: " << njh::conToStr(rec.alts_, ",") << " to ALT: " << njh::conToStr(newAlts, ",") << "\n";
		ss << e.what();
		throw std::runtime_error{ss.str()};
	}
	rec.alts_ = newAlts;
}

/**
 * \brief merge records at the same CHROM/POS/REF (e.g. called by overlapping targets) into the record with the highest NS
 *
 * By default a sample with no data in the kept record is rescued from the overlapping record with the highest depth.
 * With pars.combinedOverlappingCallsAcrossTargets the sample's depths from every overlapping record are summed instead.
 *
 * \param vcf the vcf holding the records
 * \param group the positions of the records in vcf.records_ to merge
 * \param pars the combining parameters
 * \param toErase the records that are merged away are marked true here
 */
void mergeOverlappingRecords(VCFOutput & vcf,
                             const std::vector<uint32_t> & group,
                             const VCFOutput::combineVCFsPars & pars,
                             std::vector<bool> & toErase) {
	auto getNS = [&vcf](const uint32_t pos) {
		return vcf.records_[pos].info_.containsMeta("NS") ? vcf.records_[pos].info_.getMeta<uint32_t>("NS") : 0U;
	};
	/**@todo need to adjust for when forcing a positioning with <*> and one overlapping things calls a variant and the other one does not */
	uint32_t bestPos = group.front();
	for (const auto pos : group) {
		if (getNS(pos) > getNS(bestPos)) {
			bestPos = pos;
		}
	}
	auto & best = vcf.records_[bestPos];

	//mark the rest for removal and add their targets to the kept record
	std::set<std::string> targets;
	for (const auto pos : group) {
		if (pos != bestPos) {
			toErase[pos] = true;
		}
		if (vcf.records_[pos].info_.containsMeta("TARGET")) {
			njh::addVecToSet(tokenizeString(vcf.records_[pos].info_.getMeta("TARGET"), "::"), targets);
		}
	}
	if (best.info_.containsMeta("TARGET")) {
		best.info_.addMeta("TARGET", njh::conToStr(targets, "::"), true);
	}

	if (pars.doNotRescueVariantCallsAcrossTargets && !pars.combinedOverlappingCallsAcrossTargets) {
		return;
	}

	//determine which records each sample's data will come from
	std::map<size_t, std::vector<uint32_t>> sourcesPerSample;
	for (size_t samp = 0; samp < best.sampleValues_.size(); ++samp) {
		std::vector<uint32_t> rowsWithData;
		for (const auto pos : group) {
			if (sampleHasData(vcf.records_[pos], samp)) {
				rowsWithData.emplace_back(pos);
			}
		}
		if (rowsWithData.empty()) {
			continue;
		}
		if (pars.combinedOverlappingCallsAcrossTargets) {
			if (rowsWithData.size() > 1 || rowsWithData.front() != bestPos) {
				sourcesPerSample[samp] = rowsWithData;
			}
		} else if (!sampleHasData(best, samp)) {
			//rescue from the record with the highest depth
			uint32_t bestRow = rowsWithData.front();
			for (const auto pos : rowsWithData) {
				if (getSampleDP(vcf.records_[pos], samp) > getSampleDP(vcf.records_[bestRow], samp)) {
					bestRow = pos;
				}
			}
			sourcesPerSample[samp] = {bestRow};
		}
	}
	if (sourcesPerSample.empty()) {
		return;
	}

	//put the kept record and every contributing record on the same alts before merging any sample data,
	//alts only in contributing records are appended in sorted order
	std::set<uint32_t> contributingRows;
	std::set<std::string> extraAlts;
	for (const auto & sources : sourcesPerSample) {
		for (const auto pos : sources.second) {
			contributingRows.emplace(pos);
			for (const auto & alt : vcf.records_[pos].alts_) {
				if (njh::notIn(alt, best.alts_)) {
					extraAlts.emplace(alt);
				}
			}
		}
	}
	auto newAlts = best.alts_;
	newAlts.insert(newAlts.end(), extraAlts.begin(), extraAlts.end());
	setRecordAlts(best, newAlts, vcf);
	for (const auto pos : contributingRows) {
		setRecordAlts(vcf.records_[pos], newAlts, vcf);
	}

	//merge the sample data, AN_REAL and AC_REAL are microhaplotype counts that can't be recalculated from the sample
	//data so increase them by 1 for each allele a sample has that it didn't have in the kept record
	uint32_t AN_REAL = best.info_.containsMeta("AN_REAL") ? best.info_.getMeta<uint32_t>("AN_REAL") : 0U;
	std::vector<uint32_t> AC_REAL(newAlts.size(), 0);
	if (best.info_.containsMeta("AC_REAL")) {
		const auto AC_toks = tokenizeString(best.info_.getMeta("AC_REAL"), ",");
		std::transform(AC_toks.begin(), AC_toks.end(), AC_REAL.begin(), [](const std::string & tok) {
			return njh::StrToNumConverter::stoToNum<uint32_t>(tok);
		});
	}
	for (const auto & sources : sourcesPerSample) {
		const auto & samp = sources.first;
		const auto before = sampleHasData(best, samp)
			                    ? getADCounts(best, best.getSampleField(samp, "AD"))
			                    : std::vector<uint32_t>(best.getNumberOfAlleles(), 0);
		//start from the record with the most depth so the fields that aren't summed come from the best supported call
		uint32_t baseRow = sources.second.front();
		for (const auto pos : sources.second) {
			if (getSampleDP(vcf.records_[pos], samp) > getSampleDP(vcf.records_[baseRow], samp)) {
				baseRow = pos;
			}
		}
		auto merged = vcf.records_[baseRow].getSampleMeta(samp);
		if (sources.second.size() > 1) {
			uint32_t DP = 0;
			std::vector<uint32_t> AD(best.getNumberOfAlleles(), 0);
			for (const auto pos : sources.second) {
				DP += getSampleDP(vcf.records_[pos], samp);
				const auto rowAD = getADCounts(vcf.records_[pos], vcf.records_[pos].getSampleField(samp, "AD"));
				std::transform(AD.begin(), AD.end(), rowAD.begin(), AD.begin(), std::plus<>());
			}
			std::vector<double> AF;
			AF.reserve(AD.size());
			for (const auto count : AD) {
				AF.emplace_back(count / static_cast<double>(DP));
			}
			merged.addMeta("DP", DP, true);
			merged.addMeta("AD", njh::conToStr(AD, ","), true);
			merged.addMeta("AF", njh::conToStr(AF, ","), true);
		}
		const auto after = getADCounts(best, merged.getMeta("AD"));
		for (const auto idx : iter::range(after.size())) {
			if (after[idx] > 0 && 0 == before[idx]) {
				++AN_REAL;
				if (idx > 0) {
					++AC_REAL[idx - 1];
				}
			}
		}
		best.setSampleMeta(samp, merged);
	}

	//recalculate the sample based counts from the merged sample data
	uint32_t NS = 0;
	std::vector<uint32_t> SC(newAlts.size(), 0);
	for (size_t samp = 0; samp < best.sampleValues_.size(); ++samp) {
		if (sampleHasData(best, samp)) {
			++NS;
			const auto AD = getADCounts(best, best.getSampleField(samp, "AD"));
			for (const auto idx : iter::range<size_t>(1, AD.size())) {
				if (AD[idx] > 0) {
					++SC[idx - 1];
				}
			}
		}
	}
	std::vector<double> PREV;
	std::vector<double> UNWEIGHTED_AF_REAL;
	for (const auto idx : iter::range(newAlts.size())) {
		PREV.emplace_back(NS > 0 ? SC[idx] / static_cast<double>(NS) : 0.0);
		UNWEIGHTED_AF_REAL.emplace_back(AN_REAL > 0 ? AC_REAL[idx] / static_cast<double>(AN_REAL) : 0.0);
	}
	best.info_.addMeta("NS", NS, true);
	best.info_.addMeta("SC", njh::conToStr(SC, ","), true);
	best.info_.addMeta("PREV", njh::conToStr(PREV, ","), true);
	best.info_.addMeta("AN_REAL", AN_REAL, true);
	best.info_.addMeta("AC_REAL", njh::conToStr(AC_REAL, ","), true);
	best.info_.addMeta("UNWEIGHTED_AF_REAL", njh::conToStr(UNWEIGHTED_AF_REAL, ","), true);
}

} // namespace

VCFOutput VCFOutput::combineVCFs(const std::vector<bfs::path> &vcfsFnps,
                                  const std::set<std::string> &sampleNamesSet,
                                  const combineVCFsPars &pars) {
	if (vcfsFnps.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no vcf files given to combine" << "\n";
		throw std::runtime_error{ss.str()};
	}
	//skip any file given more than once so its records aren't counted twice
	std::vector<bfs::path> uniqueVcfFnps;
	for (const auto & fnp : vcfsFnps) {
		if (njh::notIn(fnp, uniqueVcfFnps)) {
			uniqueVcfFnps.emplace_back(fnp);
		}
	}
	//combine the headers first so the samples' blank values use every FORMAT entry
	auto combined = VCFOutput::readInHeader(uniqueVcfFnps.front());
	std::vector<VCFOutput> headers;
	headers.reserve(uniqueVcfFnps.size());
	headers.emplace_back(combined);
	for (const auto idx : iter::range<size_t>(1, uniqueVcfFnps.size())) {
		auto currentVcf = VCFOutput::readInHeader(uniqueVcfFnps[idx]);
		//add in header lines the combined header doesn't have yet, an INFO or FORMAT with the same ID must match
		auto addEntries = [&uniqueVcfFnps, &idx](auto & combinedEntries, const auto & currentEntries, const std::string & entryType) {
			for (const auto & entry : currentEntries) {
				const auto search = combinedEntries.find(entry.first);
				if (combinedEntries.end() == search) {
					combinedEntries.emplace(entry);
				} else if (!(search->second == entry.second)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << entryType << " " << entry.first << " in " << uniqueVcfFnps[idx]
					   << " doesn't match the definition in the files before it" << "\n";
					throw std::runtime_error{ss.str()};
				}
			}
		};
		addEntries(combined.infoEntries_, currentVcf.infoEntries_, "INFO");
		addEntries(combined.formatEntries_, currentVcf.formatEntries_, "FORMAT");
		for (const auto & filter : currentVcf.filterEntries_) {
			if (std::none_of(combined.filterEntries_.begin(), combined.filterEntries_.end(),
			                 [&filter](const FilterEntry & other) { return other.id_ == filter.id_; })) {
				combined.filterEntries_.emplace_back(filter);
			}
		}
		for (const auto & other : currentVcf.otherHeaderMetaFields_) {
			const auto range = combined.otherHeaderMetaFields_.equal_range(other.first);
			if (std::none_of(range.first, range.second,
			                 [&other](const auto & existing) { return existing.second.meta_ == other.second.meta_; })) {
				combined.otherHeaderMetaFields_.emplace(other);
			}
		}
		for (const auto & other : currentVcf.otherHeaderValuePairs_) {
			if (!njh::in(other.first, combined.otherHeaderValuePairs_)) {
				combined.otherHeaderValuePairs_.emplace(other);
			}
		}

		//add in new contigs if any
		for (const auto & contig : currentVcf.contigEntries_) {
			if (!njh::in(contig.first, combined.contigEntries_)) {
				combined.contigEntries_.emplace(contig);
			} else if (!(combined.contigEntries_.at(contig.first) == contig.second)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "adding contig: " << contig.first << " but already have an entry for this but it doesn't match" << "\n";
				ss << "contig in master header: " << njh::json::writeAsOneLine(njh::json::toJson(combined.contigEntries_.at(contig.first))) << "\n";
				ss << "adding contig: " << njh::json::writeAsOneLine(njh::json::toJson(contig.second)) << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		headers.emplace_back(std::move(currentVcf));
	}
	//read in the records, the samples are the requested samples, the first file's samples and the samples of every file with records
	std::set<std::string> allSamples(sampleNamesSet.begin(), sampleNamesSet.end());
	njh::addVecToSet(headers.front().samples_, allSamples);
	for (auto & header : headers) {
		InputStream in(uniqueVcfFnps[&header - headers.data()]);
		header.addInRecordsFromFile(in);
		if (!header.records_.empty()) {
			njh::addVecToSet(header.samples_, allSamples);
		}
	}
	const VecStr allSamplesVec(allSamples.begin(), allSamples.end());
	combined.samples_ = allSamplesVec;
	for (auto & header : headers) {
		if (header.records_.empty()) {
			continue;
		}
		header.formatEntries_ = combined.formatEntries_;
		header.alignSamples(allSamplesVec);
		combined.records_.reserve(combined.records_.size() + header.records_.size());
		std::move(header.records_.begin(), header.records_.end(), std::back_inserter(combined.records_));
		header.records_.clear();
		header.records_.shrink_to_fit();
	}
	if (combined.records_.empty()) {
		combined.samples_ = VecStr{sampleNamesSet.begin(), sampleNamesSet.end()};
	}
	combined.sortRecords();

	//merge records at the same CHROM/POS/REF, sorting puts them next to each other
	std::vector<bool> toErase(combined.records_.size(), false);
	auto sameSite = [&combined](const uint32_t pos1, const uint32_t pos2) {
		return combined.records_[pos1].chrom_ == combined.records_[pos2].chrom_ &&
		       combined.records_[pos1].pos_ == combined.records_[pos2].pos_ &&
		       combined.records_[pos1].ref_ == combined.records_[pos2].ref_;
	};
	for (uint32_t groupStart = 0; groupStart < combined.records_.size();) {
		uint32_t groupEnd = groupStart + 1;
		while (groupEnd < combined.records_.size() && sameSite(groupStart, groupEnd)) {
			++groupEnd;
		}
		if (groupEnd - groupStart > 1) {
			std::vector<uint32_t> group(groupEnd - groupStart);
			njh::iota(group, groupStart);
			mergeOverlappingRecords(combined, group, pars, toErase);
		}
		groupStart = groupEnd;
	}
	std::vector<VCFRecord> keptRecords;
	keptRecords.reserve(combined.records_.size());
	for (const auto pos : iter::range(combined.records_.size())) {
		if (!toErase[pos]) {
			keptRecords.emplace_back(std::move(combined.records_[pos]));
		}
	}
	combined.records_ = std::move(keptRecords);

	//redetermine GT since when combining counts the GT could have changed
	combined.allAddGTFields(pars.ploidy);
	combined.allAutoAddDPFields();
	combined.allAutoAddTYPEFields();
	combined.allAutoAdd_AN_AC_AF_InfoFields();
	combined.allAutoAddWeightedAFRealField();
	combined.allAddDefaultFormatField("GQ", 40, FormatEntry("GQ", "1", "Float", "Genotype Quality"), true);
	return combined;
}
} // namespace njhseq
