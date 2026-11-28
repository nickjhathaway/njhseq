//
// Operations on whole PMOs, pmo.hpp has the classes for the format itself
//

#include "pmoUtils.hpp"
#include "njhseq/objects/BioDataObject/BioDataFileIO.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"
#include "njhseq/utils/urlUtils.hpp"
#include "njhseq/GenomeUtils/TargetGenomeAssignment.hpp"

namespace njhseq::pmo {

PortableMicrohaplotypeObject readPmo(const bfs::path & fnp) {
	InputStream in{InOptions(fnp)};
	auto ret = PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	ret.validate();
	return ret;
}

void writePmo(const PortableMicrohaplotypeObject & pmoObj, const OutOptions & outOpts) {
	OutputStream out(outOpts);
	out << pmoObj.to_json().dump() << "\n";
}

namespace {
[[noreturn]] void throwConflict(const std::string & what) {
	std::stringstream ss;
	ss << "njhseq::pmo::MinimumPmoMerger, error, " << what << "\n";
	throw std::runtime_error{ss.str()};
}

//add the fields of other to base, the same field with different values is an error
nlohmann::json mergeFields(nlohmann::json base, const nlohmann::json & other, const std::string & name) {
	for (const auto & field : other.items()) {
		if (base.contains(field.key()) && base[field.key()] != field.value()) {
			throwConflict(name + " has conflicting values for " + field.key() + ": " + base[field.key()].dump() + " and " + field.value().dump());
		}
		base[field.key()] = field.value();
	}
	return base;
}

std::map<std::string, nlohmann::json> mergeExtras(const std::map<std::string, nlohmann::json> & base, const std::map<std::string, nlohmann::json> & other, const std::string & name) {
	return mergeFields(nlohmann::json(base), nlohmann::json(other), name).get<std::map<std::string, nlohmann::json>>();
}

//a list section as json, null when it's missing
template<typename T>
nlohmann::json sectionJson(const std::optional<std::vector<T>> & section) {
	nlohmann::json ret;
	if (section.has_value()) {
		ret = nlohmann::json::array();
		for (const auto & entry : section.value()) {
			ret.push_back(entry.to_json());
		}
	}
	return ret;
}

//add the entries of other missing from base
template<typename T>
void addMissing(std::optional<std::vector<T>> & base, const std::optional<std::vector<T>> & other) {
	if (!other.has_value()) {
		return;
	}
	if (!base.has_value()) {
		base = other;
		return;
	}
	std::set<std::string> current;
	for (const auto & entry : base.value()) {
		current.emplace(entry.to_json().dump());
	}
	for (const auto & entry : other.value()) {
		if (current.emplace(entry.to_json().dump()).second) {
			base->emplace_back(entry);
		}
	}
}
}  // namespace

MinimumPmoMerger::MinimumPmoMerger(const PmoGenerationMethod & generationMethod) {
	auto today = njh::getCurrentDate();
	merged_.pmo_header_.pmo_version_ = PMO_SCHEMA_VERSION;
	merged_.pmo_header_.creation_date_ = today.substr(0, today.find('_'));
	merged_.pmo_header_.generation_method_ = generationMethod;
}

void MinimumPmoMerger::add(const PortableMicrohaplotypeObject & pmoObj) {
	if (pmoObj.project_info_.has_value() || pmoObj.sequencing_info_.has_value()) {
		throwConflict("only minimum PMOs can be merged, this PMO has project or sequencing sections");
	}
	//the bioinformatics sections have to be the same in every PMO, their one run is the run of the merged PMO
	const nlohmann::json bioinformatics{
		{"bioinformatics_methods_info", sectionJson(pmoObj.bioinformatics_methods_info_)},
		{"bioinformatics_run_info", sectionJson(pmoObj.bioinformatics_run_info_)}};
	const uint32_t runCount = pmoObj.bioinformatics_run_info_.has_value() ? pmoObj.bioinformatics_run_info_->size() : 0;
	if (!bioinformatics_.has_value()) {
		if (runCount > 1) {
			throwConflict(njh::pasteAsStr("only PMOs with at most one bioinformatics run can be merged, this PMO has ", runCount));
		}
		bioinformatics_ = bioinformatics;
		merged_.bioinformatics_methods_info_ = pmoObj.bioinformatics_methods_info_;
		merged_.bioinformatics_run_info_ = pmoObj.bioinformatics_run_info_;
	} else if (bioinformatics_.value() != bioinformatics) {
		throwConflict("PMOs with different bioinformatics_methods_info or bioinformatics_run_info can't be merged");
	}
	auto checkRunId = [&runCount](const std::optional<uint32_t> & runId, const std::string & section) {
		if (runId.has_value() && runId.value() >= runCount) {
			throwConflict(njh::pasteAsStr(section, " refers to bioinformatics run ", runId.value(), " but there are only ", runCount, " bioinformatics_run_info"));
		}
	};
	//genomes with the same name and version become one
	std::vector<uint32_t> genomeMap;
	if (pmoObj.targeted_genomes_.has_value()) {
		for (const auto & genome : pmoObj.targeted_genomes_.value()) {
			const auto key = genome.name_ + "\t" + genome.genome_version_;
			auto search = genomeIds_.find(key);
			if (genomeIds_.end() == search) {
				if (!merged_.targeted_genomes_.has_value()) {
					merged_.targeted_genomes_ = std::vector<GenomeInfo>{};
				}
				search = genomeIds_.emplace(key, merged_.targeted_genomes_->size()).first;
				merged_.targeted_genomes_->emplace_back(genome);
			} else {
				auto & current = merged_.targeted_genomes_->at(search->second);
				current = GenomeInfo::from_json(mergeFields(current.to_json(), genome.to_json(), "genome " + genome.name_ + " " + genome.genome_version_));
			}
			genomeMap.emplace_back(search->second);
		}
	}
	auto mapGenome = [&genomeMap](GenomicLocation & loc) {
		if (loc.genome_id_ >= genomeMap.size()) {
			throwConflict(njh::pasteAsStr("a genomic location refers to genome ", loc.genome_id_, " but there are only ", genomeMap.size(), " targeted_genomes"));
		}
		loc.genome_id_ = genomeMap[loc.genome_id_];
	};
	auto mapGenomes = [&mapGenome](RepresentativeMicrohaplotype & mhap) {
		if (mhap.associated_protein_variants_.has_value()) {
			for (auto & variant : mhap.associated_protein_variants_.value()) {
				mapGenome(variant.protein_location_);
				if (variant.codon_genomic_location_.has_value()) {
					mapGenome(variant.codon_genomic_location_.value());
				}
			}
		}
		if (mhap.associated_seq_variants_.has_value()) {
			for (auto & variant : mhap.associated_seq_variants_.value()) {
				mapGenome(variant);
			}
		}
	};

	std::vector<uint32_t> targetMap;
	for (const auto & target : pmoObj.target_info_) {
		auto search = targetIds_.find(target.target_name_);
		if (targetIds_.end() == search) {
			search = targetIds_.emplace(target.target_name_, merged_.target_info_.size()).first;
			auto newTarget = target;
			if (newTarget.insert_location_.has_value()) {
				mapGenome(newTarget.insert_location_.value());
			}
			if (newTarget.markers_of_interest_.has_value()) {
				for (auto & marker : newTarget.markers_of_interest_.value()) {
					mapGenome(marker.marker_location_);
				}
			}
			merged_.target_info_.emplace_back(newTarget);
		} else {
			const auto & current = merged_.target_info_[search->second];
			if (current.forward_primer_.seq_ != target.forward_primer_.seq_ || current.reverse_primer_.seq_ != target.reverse_primer_.seq_) {
				throwConflict("target " + target.target_name_ + " has different primers in different PMOs");
			}
		}
		targetMap.emplace_back(search->second);
	}

	std::vector<uint32_t> panelMap;
	for (const auto & panel : pmoObj.panel_info_) {
		auto search = panelIds_.find(panel.panel_name_);
		if (panelIds_.end() == search) {
			search = panelIds_.emplace(panel.panel_name_, merged_.panel_info_.size()).first;
			PanelInfo newPanel;
			newPanel.panel_name_ = panel.panel_name_;
			merged_.panel_info_.emplace_back(newPanel);
		}
		auto & reactions = merged_.panel_info_[search->second].reactions_;
		for (const auto & reaction : panel.reactions_) {
			auto reactionSearch = std::find_if(reactions.begin(), reactions.end(), [&reaction](const ReactionInfo & current) {
				return current.reaction_name_ == reaction.reaction_name_;
			});
			if (reactions.end() == reactionSearch) {
				ReactionInfo newReaction;
				newReaction.reaction_name_ = reaction.reaction_name_;
				reactions.emplace_back(newReaction);
				reactionSearch = std::prev(reactions.end());
			}
			std::set<uint32_t> panelTargets(reactionSearch->panel_targets_.begin(), reactionSearch->panel_targets_.end());
			for (const auto & target : reaction.panel_targets_) {
				panelTargets.emplace(targetMap.at(target));
			}
			reactionSearch->panel_targets_ = std::vector<uint32_t>(panelTargets.begin(), panelTargets.end());
		}
		panelMap.emplace_back(search->second);
	}

	std::vector<uint32_t> specimenMap;
	for (const auto & specimen : pmoObj.specimen_info_) {
		auto search = specimenIds_.find(specimen.specimen_name_);
		if (specimenIds_.end() == search) {
			search = specimenIds_.emplace(specimen.specimen_name_, merged_.specimen_info_.size()).first;
			merged_.specimen_info_.emplace_back(specimen);
		} else {
			//specimens with just a name have nothing to add, this is most specimens when merging per target PMOs
			auto specimenJson = specimen.to_json();
			if (specimenJson.size() > 1) {
				merged_.specimen_info_[search->second] = SpecimenInfo::from_json(
						mergeFields(merged_.specimen_info_[search->second].to_json(), specimenJson, "specimen " + specimen.specimen_name_));
			}
		}
		specimenMap.emplace_back(search->second);
	}

	std::vector<uint32_t> libraryMap;
	for (const auto & library : pmoObj.library_sample_info_) {
		auto mappedLibrary = library;
		mappedLibrary.specimen_id_ = specimenMap.at(library.specimen_id_);
		mappedLibrary.panel_id_ = panelMap.at(library.panel_id_);
		auto search = libraryIds_.find(library.library_sample_name_);
		if (libraryIds_.end() == search) {
			search = libraryIds_.emplace(library.library_sample_name_, merged_.library_sample_info_.size()).first;
			merged_.library_sample_info_.emplace_back(mappedLibrary);
		} else {
			auto libraryJson = mappedLibrary.to_json();
			auto currentJson = merged_.library_sample_info_[search->second].to_json();
			if (libraryJson != currentJson) {
				merged_.library_sample_info_[search->second] = LibrarySampleInfo::from_json(
						mergeFields(currentJson, libraryJson, "library sample " + library.library_sample_name_));
			}
		}
		libraryMap.emplace_back(search->second);
	}

	//a target gets one entry in representative_microhaplotypes, with one haplotype per sequence
	std::vector<uint32_t> repPosMap;
	std::vector<std::vector<uint32_t>> mhapIdMap;
	for (const auto & representatives : pmoObj.representative_microhaplotypes_.targets_) {
		const auto targetId = targetMap.at(representatives.target_id_);
		auto repSearch = repPosByTargetId_.find(targetId);
		if (repPosByTargetId_.end() == repSearch) {
			repSearch = repPosByTargetId_.emplace(targetId, merged_.representative_microhaplotypes_.targets_.size()).first;
			RepresentativeMicrohaplotypesForTarget newTarget;
			newTarget.target_id_ = targetId;
			newTarget.mhap_location_ = representatives.mhap_location_;
			if (newTarget.mhap_location_.has_value()) {
				mapGenome(newTarget.mhap_location_.value());
			}
			merged_.representative_microhaplotypes_.targets_.emplace_back(newTarget);
		}
		repPosMap.emplace_back(repSearch->second);
		auto & mhaps = merged_.representative_microhaplotypes_.targets_[repSearch->second].microhaplotypes_;
		mhapIdMap.emplace_back();
		for (const auto & mhap : representatives.microhaplotypes_) {
			auto newMhap = mhap;
			mapGenomes(newMhap);
			auto seqSearch = mhapIdBySeq_[targetId].find(mhap.seq_);
			if (mhapIdBySeq_[targetId].end() == seqSearch) {
				seqSearch = mhapIdBySeq_[targetId].emplace(mhap.seq_, mhaps.size()).first;
				//names from separate runs can be reused for different sequences, those are left unnamed
				if (newMhap.microhaplotype_name_.has_value() && !mhapNames_[targetId].emplace(newMhap.microhaplotype_name_.value()).second) {
					newMhap.microhaplotype_name_.reset();
				}
				mhaps.emplace_back(newMhap);
			} else {
				//the same sequence from another PMO, keep variants from both
				auto & current = mhaps[seqSearch->second];
				addMissing(current.associated_protein_variants_, newMhap.associated_protein_variants_);
				addMissing(current.associated_seq_variants_, newMhap.associated_seq_variants_);
			}
			mhapIdMap.back().emplace_back(seqSearch->second);
		}
	}

	for (const auto & run : pmoObj.detected_microhaplotypes_) {
		checkRunId(run.bioinformatics_run_id_, "detected_microhaplotypes");
		for (const auto & sampleResults : run.library_samples_) {
			const auto libraryId = libraryMap.at(sampleResults.library_sample_id_);
			auto detectedSearch = detectedPosByLibraryId_.find(libraryId);
			if (detectedPosByLibraryId_.end() == detectedSearch) {
				detectedSearch = detectedPosByLibraryId_.emplace(libraryId, detected_.library_samples_.size()).first;
				DetectedMicrohaplotypesForSample newSample;
				newSample.library_sample_id_ = libraryId;
				detected_.library_samples_.emplace_back(newSample);
			}
			auto & targetResults = detected_.library_samples_[detectedSearch->second].target_results_;
			for (const auto & results : sampleResults.target_results_) {
				const auto repPos = repPosMap.at(results.mhaps_target_id_);
				if (!repPosDoneByLibraryId_[libraryId].emplace(repPos).second) {
					throwConflict("library sample " + merged_.library_sample_info_[libraryId].library_sample_name_ + " has results for target "
							+ merged_.target_info_[merged_.representative_microhaplotypes_.targets_[repPos].target_id_].target_name_ + " more than once");
				}
				DetectedMicrohaplotypesForTarget newResults;
				newResults.mhaps_target_id_ = repPos;
				for (const auto & mhap : results.mhaps_) {
					auto newMhap = mhap;
					newMhap.mhap_id_ = mhapIdMap.at(results.mhaps_target_id_).at(mhap.mhap_id_);
					newResults.mhaps_.emplace_back(newMhap);
				}
				targetResults.emplace_back(newResults);
			}
		}
	}

	//read counts are per library sample and per target_info entry, so they map like library_sample_info and target_info
	if (pmoObj.read_counts_by_stage_.has_value()) {
		if (!readCounts_.has_value()) {
			readCounts_ = ReadCountsByStage{};
		}
		for (const auto & counts : pmoObj.read_counts_by_stage_.value()) {
			checkRunId(counts.bioinformatics_run_id_, "read_counts_by_stage");
			readCounts_->extras_ = mergeExtras(readCounts_->extras_, counts.extras_, "read_counts_by_stage");
			for (const auto & libraryCounts : counts.read_counts_by_library_sample_by_stage_) {
				if (libraryCounts.library_sample_id_ >= libraryMap.size()) {
					throwConflict(njh::pasteAsStr("read counts refer to library sample ", libraryCounts.library_sample_id_, " but there are only ", libraryMap.size(), " library_sample_info"));
				}
				const auto libraryId = libraryMap[libraryCounts.library_sample_id_];
				const auto & libraryName = merged_.library_sample_info_[libraryId].library_sample_name_;
				auto countsSearch = readCountsPosByLibraryId_.find(libraryId);
				if (readCountsPosByLibraryId_.end() == countsSearch) {
					countsSearch = readCountsPosByLibraryId_.emplace(libraryId, readCounts_->read_counts_by_library_sample_by_stage_.size()).first;
					ReadCountsByStageForLibrarySample newCounts;
					newCounts.library_sample_id_ = libraryId;
					newCounts.total_raw_count_ = libraryCounts.total_raw_count_;
					readCounts_->read_counts_by_library_sample_by_stage_.emplace_back(newCounts);
				}
				auto & current = readCounts_->read_counts_by_library_sample_by_stage_[countsSearch->second];
				if (current.total_raw_count_ != libraryCounts.total_raw_count_) {
					throwConflict(njh::pasteAsStr("library sample ", libraryName, " has conflicting values for total_raw_count: ", current.total_raw_count_, " and ", libraryCounts.total_raw_count_));
				}
				current.extras_ = mergeExtras(current.extras_, libraryCounts.extras_, "read counts of library sample " + libraryName);
				if (!libraryCounts.read_counts_for_targets_.has_value()) {
					continue;
				}
				if (!current.read_counts_for_targets_.has_value()) {
					current.read_counts_for_targets_ = std::vector<ReadCountsByStageForTarget>{};
				}
				for (const auto & targetCounts : libraryCounts.read_counts_for_targets_.value()) {
					if (targetCounts.target_id_ >= targetMap.size()) {
						throwConflict(njh::pasteAsStr("read counts of library sample ", libraryName, " refer to target ", targetCounts.target_id_, " but there are only ", targetMap.size(), " target_info"));
					}
					auto newCounts = targetCounts;
					newCounts.target_id_ = targetMap[targetCounts.target_id_];
					if (!readCountTargetsDoneByLibraryId_[libraryId].emplace(newCounts.target_id_).second) {
						throwConflict("library sample " + libraryName + " has read counts for target " + merged_.target_info_[newCounts.target_id_].target_name_ + " more than once");
					}
					current.read_counts_for_targets_->emplace_back(newCounts);
				}
			}
		}
	}
}

PortableMicrohaplotypeObject MinimumPmoMerger::finish() {
	const bool hasRun = merged_.bioinformatics_run_info_.has_value() && !merged_.bioinformatics_run_info_->empty();
	if (hasRun) {
		detected_.bioinformatics_run_id_ = 0;
	}
	merged_.detected_microhaplotypes_.emplace_back(detected_);
	if (readCounts_.has_value()) {
		auto & libraryCounts = readCounts_->read_counts_by_library_sample_by_stage_;
		std::sort(libraryCounts.begin(), libraryCounts.end(), [](const ReadCountsByStageForLibrarySample & a, const ReadCountsByStageForLibrarySample & b) {
			return a.library_sample_id_ < b.library_sample_id_;
		});
		if (hasRun) {
			readCounts_->bioinformatics_run_id_ = 0;
		}
		merged_.read_counts_by_stage_ = std::vector<ReadCountsByStage>{std::move(readCounts_.value())};
	}
	merged_.validate();
	return std::move(merged_);
}

PortableMicrohaplotypeObject mergeMinimumPmos(const std::vector<PortableMicrohaplotypeObject> & pmos,
		const PmoGenerationMethod & generationMethod) {
	if (pmos.empty()) {
		throwConflict("no PMOs given to merge");
	}
	MinimumPmoMerger merger(generationMethod);
	for (const auto & pmoObj : pmos) {
		merger.add(pmoObj);
	}
	return merger.finish();
}

AddSpecimenMetaResults addSpecimenMeta(PortableMicrohaplotypeObject & pmoObj, const table & meta, const std::string & sampleCol) {
	const auto sampleColPos = meta.getColPos(sampleCol);
	std::unordered_map<std::string, uint32_t> specimenBySample;
	for (const auto specimenId : iter::range(pmoObj.specimen_info_.size())) {
		specimenBySample[pmoObj.specimen_info_[specimenId].specimen_name_] = specimenId;
	}
	for (const auto & library : pmoObj.library_sample_info_) {
		specimenBySample[library.library_sample_name_] = library.specimen_id_;
	}
	AddSpecimenMetaResults ret;
	std::set<std::string> seenSamples;
	for (const auto & row : meta.content_) {
		const auto & sample = row[sampleColPos];
		if (!seenSamples.emplace(sample).second) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, sample " << sample << " is in the meta table more than once\n";
			throw std::runtime_error{ss.str()};
		}
		auto search = specimenBySample.find(sample);
		if (specimenBySample.end() == search) {
			ret.metaSamplesNotInPmo_.emplace_back(sample);
			continue;
		}
		auto specimenJson = pmoObj.specimen_info_[search->second].to_json();
		for (const auto colPos : iter::range(meta.columnNames_.size())) {
			if (colPos == sampleColPos || njh::in(row[colPos], PMO_NA_STRINGS())) {
				continue;
			}
			//schema fields can be numbers, booleans or lists rather than text
			const auto & colName = meta.columnNames_[colPos];
			auto candidate = specimenJson;
			candidate[colName] = row[colPos];
			bool added = false;
			std::string lastError;
			std::vector<std::function<nlohmann::json()>> asValues{
				[&row, &colPos]() {return nlohmann::json(row[colPos]);},
				[&row, &colPos]() {return nlohmann::json::parse(row[colPos]);},
				[&row, &colPos]() {return nlohmann::json::array({nlohmann::json::parse(row[colPos])});}};
			for (const auto & asValue : asValues) {
				try {
					candidate[colName] = asValue();
					SpecimenInfo::from_json(candidate);
					added = true;
					break;
				} catch (const std::exception & e) {
					lastError = e.what();
				}
			}
			if (!added) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, " << row[colPos] << " isn't a valid value for specimen field " << colName << " for " << sample << ": " << lastError << "\n";
				throw std::runtime_error{ss.str()};
			}
			specimenJson = candidate;
		}
		pmoObj.specimen_info_[search->second] = SpecimenInfo::from_json(specimenJson);
		++ret.specimensUpdated_;
	}
	pmoObj.validate();
	return ret;
}

namespace {
//the genes of a gff by chromosome, as GENE_ID:GENE_NAME
std::unordered_map<std::string, std::vector<std::pair<GenomicRegion, std::string>>> readGenes(const bfs::path & gffFnp) {
	std::unordered_map<std::string, std::vector<std::pair<GenomicRegion, std::string>>> ret;
	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(gffFnp))};
	reader.openIn();
	std::string line;
	auto gRecord = reader.readNextRecord();
	while (nullptr != gRecord) {
		if ("gene" == gRecord->type_ || "pseudogene" == gRecord->type_ || njh::endsWith(gRecord->type_, "_gene")) {
			auto name = gRecord->hasAttr("Name") ? urldecode(gRecord->getAttr("Name")) : std::string("NA");
			if (name.empty()) {
				name = "NA";
			}
			ret[gRecord->seqid_].emplace_back(GenomicRegion(*gRecord), njh::pasteAsStr(urldecode(gRecord->getAttr("ID")), ":", name));
		}
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
		gRecord = reader.readNextRecord();
	}
	for (auto & chromGenes : ret) {
		std::stable_sort(chromGenes.second.begin(), chromGenes.second.end(), [](const auto & a, const auto & b) {
			return a.first.start_ < b.first.start_;
		});
	}
	return ret;
}
}  // namespace

AddTargetGenomeInfoResults addTargetGenomeInfo(PortableMicrohaplotypeObject & pmoObj, const AddTargetGenomeInfoPars & pars) {
	const auto locations = TargetGenomeLocations::read(pars.locationsDir_);
	TargetGenomeChoicePars choicePars;
	choicePars.preferredGenomes_ = pars.preferredGenomes_;
	const auto chosen = chooseTargetGenomes(locations, choicePars);
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::pair<GenomicRegion, std::string>>>> genesByGenome;
	std::unordered_map<std::string, std::unique_ptr<TwoBit::TwoBitFile>> twoBits;

	AddTargetGenomeInfoResults ret;
	for (auto & target : pmoObj.target_info_) {
		auto chosenSearch = chosen.find(locations.baseTarget(target.target_name_));
		if (chosen.end() == chosenSearch) {
			ret.targetsNotFound_.emplace_back(target.target_name_);
			continue;
		}
		const auto & hit = chosenSearch->second;

		if (!pmoObj.targeted_genomes_.has_value()) {
			pmoObj.targeted_genomes_ = std::vector<GenomeInfo>{};
		}
		auto & genomes = pmoObj.targeted_genomes_.value();
		auto genomeSearch = std::find_if(genomes.begin(), genomes.end(), [&hit](const GenomeInfo & genome) {
			return genome.name_ == hit.genome_;
		});
		if (genomes.end() == genomeSearch) {
			GenomeInfo genome;
			genome.name_ = hit.genome_;
			genome.genome_version_ = hit.genome_;
			genome.url_ = "NA";
			genome.validate();
			genomes.emplace_back(genome);
			genomeSearch = std::prev(genomes.end());
		}
		const uint32_t genomeId = std::distance(genomes.begin(), genomeSearch);
		auto genLocation = [&genomeId, &hit](uint32_t start, uint32_t end, char strand) {
			GenomicLocation loc;
			loc.genome_id_ = genomeId;
			loc.chrom_ = hit.chrom_;
			loc.start_ = start;
			loc.end_ = end;
			loc.strand_ = std::string(1, strand);
			return loc;
		};
		//mates analyzed as separate targets cover their own window of the insert
		uint32_t insertStart = hit.insertStart_;
		uint32_t insertEnd = hit.insertEnd_;
		char insertStrand = hit.strand_;
		if (target.target_name_ != chosenSearch->first) {
			const auto mateRegion = locations.innerRegion(target.target_name_, hit);
			if (nullptr != mateRegion) {
				insertStart = mateRegion->chromStart_;
				insertEnd = mateRegion->chromEnd_;
				insertStrand = mateRegion->strand_;
			}
		}
		const GenomicRegion insertRegion(target.target_name_, hit.chrom_, insertStart, insertEnd, '-' == insertStrand);
		auto insert = genLocation(insertStart, insertEnd, insertStrand);
		if (!pars.genomeDir_.empty()) {
			auto twoBitSearch = twoBits.find(hit.genome_);
			if (twoBits.end() == twoBitSearch) {
				auto twoBitFnp = njh::files::make_path(pars.genomeDir_, hit.genome_ + ".2bit");
				twoBitSearch = twoBits.emplace(hit.genome_, bfs::exists(twoBitFnp) ? std::make_unique<TwoBit::TwoBitFile>(twoBitFnp) : nullptr).first;
			}
			if (nullptr != twoBitSearch->second) {
				insert.ref_seq_ = insertRegion.extractSeq(*twoBitSearch->second).seq_;
			}
		}
		target.insert_location_ = insert;
		target.forward_primer_.location_ = genLocation(hit.fwdStart_, hit.fwdEnd_, hit.strand_);
		target.reverse_primer_.location_ = genLocation(hit.revStart_, hit.revEnd_, '+' == hit.strand_ ? '-' : '+');

		target.gene_name_.reset();
		if (!pars.gffDir_.empty()) {
			auto genesSearch = genesByGenome.find(hit.genome_);
			if (genesByGenome.end() == genesSearch) {
				auto gffFnp = njh::files::make_path(pars.gffDir_, hit.genome_ + ".gff");
				genesSearch = genesByGenome.emplace(hit.genome_, bfs::exists(gffFnp) ? readGenes(gffFnp)
						: std::unordered_map<std::string, std::vector<std::pair<GenomicRegion, std::string>>>{}).first;
			}
			VecStr geneNames;
			auto chromSearch = genesSearch->second.find(hit.chrom_);
			if (genesSearch->second.end() != chromSearch) {
				for (const auto & gene : chromSearch->second) {
					if (gene.first.overlaps(insertRegion) && !njh::in(gene.second, geneNames)) {
						geneNames.emplace_back(gene.second);
					}
				}
			}
			if (!geneNames.empty()) {
				target.gene_name_ = njh::conToStr(geneNames, ",");
				++ret.targetsWithGenes_;
			}
		}
		target.validate();
		++ret.targetsPerGenome_[hit.genome_];
	}
	pmoObj.validate();
	return ret;
}

}  // namespace njhseq::pmo
