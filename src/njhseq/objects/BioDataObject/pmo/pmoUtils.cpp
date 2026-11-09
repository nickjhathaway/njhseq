//
// Operations on whole PMOs, pmo.hpp has the classes for the format itself
//

#include "pmoUtils.hpp"

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
	if (pmoObj.project_info_.has_value() || pmoObj.sequencing_info_.has_value()
			|| pmoObj.bioinformatics_methods_info_.has_value() || pmoObj.bioinformatics_run_info_.has_value()
			|| pmoObj.read_counts_by_stage_.has_value()) {
		throwConflict("only minimum PMOs can be merged, this PMO has project, sequencing, bioinformatics or read count sections");
	}
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
}

PortableMicrohaplotypeObject MinimumPmoMerger::finish() {
	merged_.detected_microhaplotypes_.emplace_back(detected_);
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

}  // namespace njhseq::pmo
