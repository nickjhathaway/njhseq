//
// Created by Nicholas Hathaway on 1/30/24.
//

#include "VCFOutput.hpp"

namespace njhseq {

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


void VCFOutput::VCFRecord::autoAdd_AN_AC_AF_InfoFields() {
	std::unordered_map<std::string, uint32_t> alleleCounts;
	uint32_t AN = 0;
	std::regex blankDataPattern("\\.(,\\.)*");
	for (auto& samp: sampleFormatInfos_) {
		if (!samp.second.containsMeta("GT")) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field GT in field if going to add total AN, AC, and AF fields " << "\n";
			ss << "Current fields are : " << njh::conToStr(njh::getVecOfMapKeys(samp.second.meta_), ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		if ("." != samp.second.getMeta("GT")) {
			auto GT = samp.second.getMeta("GT");
			auto GT_toks = tokenizeString(GT, "/");
			for (const auto& tok: GT_toks) {
				++alleleCounts[tok];
				++AN;
			}
		}
	}
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	// for (const auto & allele_count : alleleCounts) {
	// 	std::cout << allele_count.first << " : " << allele_count.second << std::endl;
	// }
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto alleles = getVectorOfMapKeys(alleleCounts);
	std::vector<uint32_t> allelesNumeric;
	allelesNumeric.reserve(alleles.size());
	for (auto& allele: alleles) {
		allelesNumeric.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(allele));
	}
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if (vectorMaximum(allelesNumeric) > alts_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " genotype number: " <<
				vectorMaximum(allelesNumeric) << " can't be more than the alts_.size(): " << alts_.size() << " + 1 " << "\n";
		throw std::runtime_error{ss.str()};
	}
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::vector<double> afs;
	std::string ACs;
	//the allele counts include the reference so alts.size() will be equal to the genotype number, e.g. for alts of size of 2, there will be 0, 1, 2 (will skip 0)
	for (const auto pos: iter::range(1UL, alts_.size() + 1)) {
		if (!ACs.empty()) {
			ACs += ",";
		}
		//if missing from the counts will default to zero
		ACs += estd::to_string(alleleCounts[estd::to_string(pos)]);
		afs.emplace_back(alleleCounts[estd::to_string(pos)] / static_cast<double>(AN));
	}
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	info_.addMeta("AN", AN, true);
	info_.addMeta("AC", ACs, true);
	info_.addMeta("AF", njh::conToStr(afs, ","), true);
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
}


void VCFOutput::VCFRecord::autoAddWeightedAFRealField() {
	std::vector<double> weighted_real_afs(alts_.size(), 0); //total read per each alt observation
	std::regex blankDataPattern("\\.(,\\.)*");
	double reference_sum = 0;
	for (auto &samp: sampleFormatInfos_) {
		if (!samp.second.containsMeta("DP") || !samp.second.containsMeta("AD")) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field DP and AD in field if going to add weighted AF_REAL " << "\n";
			ss << "Current fields are : " << njh::conToStr(njh::getVecOfMapKeys(samp.second.meta_), ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		if("." != samp.second.getMeta("DP")) {
			double DP = samp.second.getMeta<uint32_t>("DP");
			auto AD_toks = tokenizeString(samp.second.getMeta("AD"), ",");
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
	}
	auto weighted_real_afs_sum = vectorSum(weighted_real_afs) + reference_sum;
	for (auto & af : weighted_real_afs) {
		af /= weighted_real_afs_sum;
	}
	info_.addMeta("AF_REAL", njh::conToStr(weighted_real_afs, ","), true);
}

void VCFOutput::VCFRecord::autoAddTotalDP_RO_AO_InfoFields() {
	uint32_t totalDP = 0;
	uint32_t totalRO = 0; //total reads for reference
	std::vector<uint32_t> totalAOs(alts_.size(), 0); //total read per each alt observation
	std::regex blankDataPattern("\\.(,\\.)*");
	for (auto &samp: sampleFormatInfos_) {
		if (!samp.second.containsMeta("DP") || !samp.second.containsMeta("AD")) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field DP and AD in field if going to add total DP, RO, and AO fields " << "\n";
			ss << "Current fields are : " << njh::conToStr(njh::getVecOfMapKeys(samp.second.meta_), ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		if("." != samp.second.getMeta("DP")) {
			totalDP += samp.second.getMeta<uint32_t>("DP");
			auto AD_toks = tokenizeString(samp.second.getMeta("AD"), ",");
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
	// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
	std::regex blankDataPattern("[\\.0](,[\\.0])*");
	for (auto &samp: sampleFormatInfos_) {
		if (!samp.second.containsMeta("AD")) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " <<
					" need to have field AD in field if going to add GT field " << "\n";
			ss << "Current fields are : " << njh::conToStr(njh::getVecOfMapKeys(samp.second.meta_), ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto raw_AD = samp.second.getMeta("AD");
		if (std::regex_match(raw_AD, blankDataPattern)) {
			samp.second.addMeta("GT", ".", true);
		} else {
			std::vector<uint32_t> depths(1 + alts_.size());
			auto ADs = tokenizeString(raw_AD, ",");
			// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
			// std::cout << "ADs: " << raw_AD << std::endl;
			std::transform(ADs.begin(), ADs.end(), depths.begin(), [](const auto &p) {
				return njh::StrToNumConverter::stoToNum<uint32_t>(p);
			});
			// std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			std::vector<uint32_t> rank(ADs.size());
			njh::iota(rank, 0U);
			njh::sort(rank, [&depths](const auto &pos1, const auto &pos2) {
				return depths[pos1] > depths[pos2];
			});

			//count alleles with depth
			uint32_t allelesWithDepth = std::count_if(depths.begin(), depths.end(), [](const auto d) { return d > 0; });
			// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
			// std::cout << "sample: " << samp.first << std::endl;
			// std::cout << "ploidy: " << ploidy << std::endl;
			// std::cout << "depths: " << njh::conToStr(depths, ",") << std::endl;
			// std::cout << "allelesWithDepth: " << allelesWithDepth << std::endl;
			// std::cout << "rank: " << njh::conToStr(rank, ",") << std::endl;
			// std::cout << "ploidy <= allelesWithDepth: " << njh::colorBool(ploidy <= allelesWithDepth) << std::endl;

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
				// std::cout << "diff: " << diff << std::endl;
				// std::cout << "gts: " << njh::conToStr(gts, ",") << std::endl;
				// std::cout << "gts.size(): " << gts.size()<< std::endl;
				for (const auto pos: iter::range(rank.size())) {
					// std::cout << "\tdiff: " << diff << std::endl;
					// std::cout << "\tfreqs[rank[pos]]: " << freqs[rank[pos]] << std::endl;
					// std::cout << "\tstd::round(freqs[rank[pos]] * diff)): " << std::round(freqs[rank[pos]] * diff) << std::endl;
					auto add = std::min(diff, static_cast<uint32_t>(std::round(freqs[rank[pos]] * diff)));
					// std::cout << "\tadd: " << add << std::endl;
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
			auto GT = njh::conToStr(gts, "/");
			samp.second.addMeta("GT", GT, true);
			// std::cout << "GT: " << GT << std::endl;
			// std::cout << std::endl;
			// if(1 == gts.size()) {
			// 	exit(1);
			// }
		}
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
	uint32_t end = pos_;
	if(std::all_of(alts_.begin(), alts_.end(), [this](const std::string& alt) {
		return ref_.size() == alt.size();
	})) {
		end = start + ref_.size();
	} else if (ref_.size() == 1 && !std::all_of(alts_.begin(), alts_.end(), [](const std::string& alt) {
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
	njh::for_each(
		records_, [](auto& rec) { rec.autoAddWeightedAFRealField(); });

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
		<< "\t" << (rec.qual_ == std::numeric_limits<uint32_t>::max() ? "." : estd::to_string(rec.qual_))
		<< "\t" << rec.filter_;
		std::string infoOut;
		for (const auto & infoKey: infoEntries_) {
			const auto & info = infoKey.second;
			if (infoKey.second.type_ == "Flag") {
				if (rec.info_.containsMeta(info.id_)) {
					if(!infoOut.empty()) {
						infoOut +=";";
					}
					infoOut += info.id_;
				}
			} else {
				if(!infoOut.empty()) {
					infoOut +=";";
				}
				infoOut += info.id_ + "=" + rec.info_.getMeta(info.id_);
			}
		}
		vcfOut << "\t" << infoOut;
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
	for (const auto & infoKey: formatOutputOrder) {
		const auto & format = formatEntries_.at(infoKey);
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
			<< "\t" << (rec.qual_ == std::numeric_limits<uint32_t>::max() ? "." : estd::to_string(rec.qual_))
			<< "\t" << rec.filter_;
			std::string infoOut;
			for (const auto & infoKey: infoEntries_) {
				const auto & info = infoKey.second;

				if (infoKey.second.type_ == "Flag") {
					if (rec.info_.containsMeta(info.id_)) {
						if(!infoOut.empty()) {
							infoOut +=";";
						}
						infoOut += info.id_;
					}
				} else {
					if(!infoOut.empty()) {
						infoOut +=";";
					}
					infoOut += info.id_ + "=" + rec.info_.getMeta(info.id_);
				}
			}
			vcfOut << "\t" << infoOut;
			vcfOut << "\t" << formatOut;
			for(const auto & sampleName : samples_) {
				const auto & sample = rec.sampleFormatInfos_.at(sampleName);
				std::string formatOutForSample;
				for (const auto & infoKey: formatOutputOrder) {
					const auto & format = formatEntries_.at(infoKey);
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

	auto infoToks = tokenizeString(toks[7], ";");
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


VCFOutput VCFOutput::comnbineVCFs(const std::vector<bfs::path> &vcfsFnps,
																	const comnbineVCFsPars &pars) {
	std::set<std::string> sampleNamesSet;
	for (const auto & fnp : vcfsFnps) {
		auto vcfHeader = VCFOutput::readInHeader(fnp);
		njh::addVecToSet(vcfHeader.samples_, sampleNamesSet);
	}
	return comnbineVCFs(vcfsFnps, sampleNamesSet, pars);
}

namespace {

/**
 * \brief whether a VCF value is missing, e.g. "." or ".,.,."
 */
bool isMissingVCFValue(const std::string & val) {
	return !val.empty() && std::all_of(val.begin(), val.end(), [](const char c) { return '.' == c || ',' == c; });
}

/**
 * \brief whether a sample has read data for a record, i.e. has a DP that is set and non-zero
 */
bool sampleHasData(const VCFOutput::VCFRecord & rec, const std::string & sample) {
	const auto sampInfo = rec.sampleFormatInfos_.find(sample);
	if (rec.sampleFormatInfos_.end() == sampInfo || !sampInfo->second.containsMeta("DP")) {
		return false;
	}
	const auto DP = sampInfo->second.getMeta("DP");
	return !isMissingVCFValue(DP) && "0" != DP;
}

uint32_t getSampleDP(const VCFOutput::VCFRecord & rec, const std::string & sample) {
	return sampleHasData(rec, sample) ? rec.sampleFormatInfos_.at(sample).getMeta<uint32_t>("DP") : 0;
}

/**
 * \brief get the AD field (ref, alt1, alt2, etc.) as counts, missing values are treated as 0
 */
std::vector<uint32_t> getADCounts(const VCFOutput::VCFRecord & rec, const MetaDataInName & sampInfo) {
	const auto toks = tokenizeString(sampInfo.getMeta("AD"), ",");
	if (toks.size() != rec.getNumberOfAlleles()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "AD: " << sampInfo.getMeta("AD") << " should have " << rec.getNumberOfAlleles()
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
		for (auto & sampInfo : rec.sampleFormatInfos_) {
			//samples without data get missing values for the new alleles rather than 0
			const bool hasData = sampleHasData(rec, sampInfo.first);
			for (const auto & format : header.formatEntries_) {
				if (isPerAllele(format.second.number_) && sampInfo.second.containsMeta(format.first)) {
					currentField = "FORMAT/" + format.first + " for sample " + sampInfo.first;
					sampInfo.second.addMeta(format.first,
					                        remapAlleleValues(sampInfo.second.getMeta(format.first), rec.alts_, newAlts, "R" == format.second.number_,
					                                          hasData ? numericFill(format.second.type_) : std::string(".")),
					                        true);
				}
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
                             const VCFOutput::comnbineVCFsPars & pars,
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
	std::map<std::string, std::vector<uint32_t>> sourcesPerSample;
	for (const auto & sampInfo : best.sampleFormatInfos_) {
		const auto & samp = sampInfo.first;
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
			                    ? getADCounts(best, best.sampleFormatInfos_.at(samp))
			                    : std::vector<uint32_t>(best.getNumberOfAlleles(), 0);
		//start from the record with the most depth so the fields that aren't summed come from the best supported call
		uint32_t baseRow = sources.second.front();
		for (const auto pos : sources.second) {
			if (getSampleDP(vcf.records_[pos], samp) > getSampleDP(vcf.records_[baseRow], samp)) {
				baseRow = pos;
			}
		}
		auto merged = vcf.records_[baseRow].sampleFormatInfos_.at(samp);
		if (sources.second.size() > 1) {
			uint32_t DP = 0;
			std::vector<uint32_t> AD(best.getNumberOfAlleles(), 0);
			for (const auto pos : sources.second) {
				DP += getSampleDP(vcf.records_[pos], samp);
				const auto rowAD = getADCounts(vcf.records_[pos], vcf.records_[pos].sampleFormatInfos_.at(samp));
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
		const auto after = getADCounts(best, merged);
		for (const auto idx : iter::range(after.size())) {
			if (after[idx] > 0 && 0 == before[idx]) {
				++AN_REAL;
				if (idx > 0) {
					++AC_REAL[idx - 1];
				}
			}
		}
		best.sampleFormatInfos_[samp] = merged;
	}

	//recalculate the sample based counts from the merged sample data
	uint32_t NS = 0;
	std::vector<uint32_t> SC(newAlts.size(), 0);
	for (const auto & sampInfo : best.sampleFormatInfos_) {
		if (sampleHasData(best, sampInfo.first)) {
			++NS;
			const auto AD = getADCounts(best, sampInfo.second);
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

VCFOutput VCFOutput::comnbineVCFs(const std::vector<bfs::path> &vcfsFnps,
                                  const std::set<std::string> &sampleNamesSet,
                                  const comnbineVCFsPars &pars) {
	if (vcfsFnps.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no vcf files given to combine" << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto readInVcf = [](const bfs::path & fnp) {
		auto ret = VCFOutput::readInHeader(fnp);
		InputStream in(fnp);
		ret.addInRecordsFromFile(in);
		return ret;
	};
	//skip any file given more than once so its records aren't counted twice
	std::vector<bfs::path> uniqueVcfFnps;
	for (const auto & fnp : vcfsFnps) {
		if (njh::notIn(fnp, uniqueVcfFnps)) {
			uniqueVcfFnps.emplace_back(fnp);
		}
	}
	auto combined = readInVcf(uniqueVcfFnps.front());
	for (const auto idx : iter::range<size_t>(1, uniqueVcfFnps.size())) {
		auto currentVcf = readInVcf(uniqueVcfFnps[idx]);
		//add in header
		//given we know what generated these vcf files there's not need to check the FORMAT, FILTER, or INFO tags

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
		combined.records_.reserve(combined.records_.size() + currentVcf.records_.size());
		std::move(currentVcf.records_.begin(), currentVcf.records_.end(), std::back_inserter(combined.records_));
	}
	//add in blanks for any missing samples for any of the records
	combined.addInBlnaksForAnyMissingSamples(sampleNamesSet);
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
