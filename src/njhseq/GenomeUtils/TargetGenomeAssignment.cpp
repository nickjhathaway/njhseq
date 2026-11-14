//
// Choosing the one genome each target is placed in when SeekDeep genTargetInfoFromGenomes found it in several genomes
//

#include "TargetGenomeAssignment.hpp"
#include "njhseq/objects/BioDataObject/reading.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"

namespace njhseq {

namespace {
//the genomes of a locations directory by the files with suffix, empty files are skipped
std::map<std::string, bfs::path> genomeFiles(const bfs::path & locationsDir, const std::string & suffix) {
	std::map<std::string, bfs::path> ret;
	for (const auto & file : njh::files::listAllFiles(locationsDir, false, VecStr{suffix})) {
		const auto fileName = file.first.filename().string();
		if (file.second || !njh::endsWith(fileName, suffix) || 0 == bfs::file_size(file.first)) {
			continue;
		}
		ret.emplace(fileName.substr(0, fileName.size() - suffix.size()), file.first);
	}
	return ret;
}
}  // namespace

TargetGenomeLocations TargetGenomeLocations::read(const bfs::path & locationsDir) {
	TargetGenomeLocations ret;
	for (const auto & genomeFile : genomeFiles(locationsDir, "_infos.tab.txt")) {
		table infos(genomeFile.second.string(), "\t", true);
		for (const auto & row : infos) {
			auto num = [&infos, &row](const std::string & col) {
				return njh::StrToNumConverter::stoToNum<uint32_t>(row[infos.getColPos(col)]);
			};
			TargetGenomeHit hit;
			hit.genome_ = genomeFile.first;
			hit.chrom_ = row[infos.getColPos("#chrom")];
			hit.strand_ = row[infos.getColPos("strand")].front();
			hit.insertStart_ = num("insertStart");
			hit.insertEnd_ = num("insertStop");
			hit.fwdStart_ = num("fPrimerStart");
			hit.fwdEnd_ = num("fPrimerStop");
			hit.revStart_ = num("rPrimerStart");
			hit.revEnd_ = num("rPrimerStop");
			hit.primerErrors_ = num("errorsInFPrimer") + num("errorsInRPrimer");
			ret.hitsByTarget_[row[infos.getColPos("target")]].emplace_back(hit);
		}
	}
	for (const auto & genomeFile : genomeFiles(locationsDir, "_inner.bed")) {
		ret.innerRegionsByGenome_[genomeFile.first] = getBeds(genomeFile.second);
	}
	return ret;
}

std::string TargetGenomeLocations::baseTarget(const std::string & name) const {
	if (!njh::in(name, hitsByTarget_) && (njh::endsWith(name, "-R1") || njh::endsWith(name, "-R2"))) {
		return name.substr(0, name.size() - 3);
	}
	return name;
}

std::shared_ptr<Bed6RecordCore> TargetGenomeLocations::innerRegion(const std::string & name, const TargetGenomeHit & hit) const {
	auto regionsSearch = innerRegionsByGenome_.find(hit.genome_);
	if (innerRegionsByGenome_.end() == regionsSearch) {
		return nullptr;
	}
	for (const auto & region : regionsSearch->second) {
		if (region->name_ == name && region->chrom_ == hit.chrom_ && region->chromStart_ < hit.insertEnd_ && hit.insertStart_ < region->chromEnd_) {
			return region;
		}
	}
	return nullptr;
}

std::map<std::string, TargetGenomeHit> chooseTargetGenomes(const TargetGenomeLocations & locations, const TargetGenomeChoicePars & pars) {
	auto preferredRank = [&pars](const std::string & genome) {
		return static_cast<uint32_t>(std::distance(pars.preferredGenomes_.begin(),
				std::find(pars.preferredGenomes_.begin(), pars.preferredGenomes_.end(), genome)));
	};
	const auto notPreferred = static_cast<uint32_t>(pars.preferredGenomes_.size());
	std::map<std::string, TargetGenomeHit> ret;
	for (const auto & targetHits : locations.hitsByTarget_) {
		const auto & target = targetHits.first;
		std::vector<TargetGenomeHit> candidates;
		for (const auto & hit : targetHits.second) {
			if (!pars.onlyPreferred_ || preferredRank(hit.genome_) < notPreferred) {
				candidates.emplace_back(hit);
			}
		}
		if (candidates.empty()) {
			continue;
		}
		ret.emplace(target, *std::min_element(candidates.begin(), candidates.end(),
				[&preferredRank, &target](const TargetGenomeHit & a, const TargetGenomeHit & b) {
			const bool aPrefix = njh::beginsWith(target, a.genome_);
			const bool bPrefix = njh::beginsWith(target, b.genome_);
			return std::make_tuple(preferredRank(a.genome_), a.primerErrors_, !aPrefix, a.genome_, a.chrom_, a.insertStart_)
					< std::make_tuple(preferredRank(b.genome_), b.primerErrors_, !bPrefix, b.genome_, b.chrom_, b.insertStart_);
		}));
	}
	return ret;
}

}  // namespace njhseq
