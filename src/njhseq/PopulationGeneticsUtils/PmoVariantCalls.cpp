//
// Adding the variant calls of collapseAndCallVariants to a PMO
//

#include "PmoVariantCalls.hpp"
#include "njhseq/IO/SeqIO/SeqInput.hpp"

namespace njhseq {

void PmoGenomePars::setOptions(seqSetUp & setUp) {
	setUp.setOption(name_, "--pmoGenomeName", "Name of the genome variants were called against for the PMO, defaults to the genome file name without its extension", false, "PMO");
	setUp.setOption(version_, "--pmoGenomeVersion", "Version of the genome for the PMO, defaults to the genome name", false, "PMO");
	setUp.setOption(taxonIds_, "--pmoGenomeTaxonIds", "NCBI taxonomy ids of the genome for the PMO, e.g. 36329", false, "PMO");
	setUp.setOption(url_, "--pmoGenomeUrl", "Where the genome can be downloaded, for the PMO", false, "PMO");
	setUp.setOption(gffUrl_, "--pmoGffUrl", "Where the genome's annotation can be downloaded, for the PMO", false, "PMO");
}

pmo::GenomeInfo PmoGenomePars::genGenomeInfo(const bfs::path & genomeFnp) const {
	pmo::GenomeInfo ret;
	ret.name_ = name_.empty() ? bfs::path(genomeFnp).filename().replace_extension("").string() : name_;
	ret.genome_version_ = version_.empty() ? ret.name_ : version_;
	ret.taxon_id_ = taxonIds_;
	ret.url_ = url_;
	if (!gffUrl_.empty()) {
		ret.gff_url_ = gffUrl_;
	}
	ret.validate();
	return ret;
}

AddVariantCallsToPmoResults addVariantCallsToPmo(pmo::PortableMicrohaplotypeObject & pmoObj,
		const std::string & targetName,
		const CollapseAndCallVariantsPars & varCallPars,
		const TranslatorByAlignment::TranslatorByAlignmentResult & calls,
		const PmoGenomePars & genomePars) {
	auto targetSearch = std::find_if(pmoObj.target_info_.begin(), pmoObj.target_info_.end(), [&targetName](const pmo::TargetInfo & target) {
		return target.target_name_ == targetName;
	});
	if (pmoObj.target_info_.end() == targetSearch) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, target " << targetName << " isn't in the PMO\n";
		throw std::runtime_error{ss.str()};
	}
	const uint32_t targetId = std::distance(pmoObj.target_info_.begin(), targetSearch);
	auto repSearch = std::find_if(pmoObj.representative_microhaplotypes_.targets_.begin(), pmoObj.representative_microhaplotypes_.targets_.end(),
			[&targetId](const pmo::RepresentativeMicrohaplotypesForTarget & reps) {
		return reps.target_id_ == targetId;
	});
	AddVariantCallsToPmoResults ret;
	if (pmoObj.representative_microhaplotypes_.targets_.end() == repSearch) {
		return ret;
	}

	//the genome, reusing an entry with the same name and version
	const auto genome = genomePars.genGenomeInfo(varCallPars.transPars.lzPars_.genomeFnp);
	if (!pmoObj.targeted_genomes_.has_value()) {
		pmoObj.targeted_genomes_ = std::vector<pmo::GenomeInfo>{};
	}
	auto & genomes = pmoObj.targeted_genomes_.value();
	auto genomeSearch = std::find_if(genomes.begin(), genomes.end(), [&genome](const pmo::GenomeInfo & current) {
		return current.name_ == genome.name_ && current.genome_version_ == genome.genome_version_;
	});
	if (genomes.end() == genomeSearch) {
		genomes.emplace_back(genome);
		genomeSearch = std::prev(genomes.end());
	}
	const uint32_t genomeId = std::distance(genomes.begin(), genomeSearch);

	std::unordered_map<std::string, std::string> callNameBySeq;
	{
		SeqInput reader(SeqIOOptions::genFastaInGz(njh::files::make_path(varCallPars.outputDirectory, "uniqueSeqs.fasta.gz")));
		reader.openIn();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			callNameBySeq[seq.seq_] = seq.name_;
		}
	}
	std::unordered_map<std::string, std::string> geneByTranscript;
	for (const auto & gene : calls.transcriptInfosForGene_) {
		for (const auto & transcript : gene.second) {
			geneByTranscript[transcript.first] = gene.first;
		}
	}
	//the gene's Name from the gff, characters the schema doesn't allow become _
	std::unordered_map<std::string, std::string> altNameByTranscript;
	for (const auto & transcript : calls.translationInfoForTranscirpt_) {
		auto name = transcript.second->geneName_;
		if (name.empty() || njh::in(name, pmo::PMO_NA_STRINGS()) || name == transcript.second->geneID_) {
			continue;
		}
		for (auto & letter : name) {
			if (!std::isalnum(static_cast<unsigned char>(letter)) && '-' != letter && '.' != letter && '_' != letter) {
				letter = '_';
			}
		}
		altNameByTranscript[transcript.first] = name;
	}
	auto genLocation = [&genomeId](const std::string & chrom, const uint32_t start, const uint32_t end,
			const std::string & refSeq, const std::string & altSeq, const std::optional<char> strand) {
		pmo::GenomicLocation loc;
		loc.genome_id_ = genomeId;
		loc.chrom_ = chrom;
		loc.start_ = start;
		loc.end_ = end;
		loc.ref_seq_ = refSeq;
		loc.alt_seq_ = altSeq;
		if (strand.has_value()) {
			loc.strand_ = std::string(1, strand.value());
		}
		return loc;
	};
	auto codonStr = [](const TranslatorByAlignment::Codon & codon) {
		return njh::pasteAsStr(std::get<0>(codon.bases_), std::get<1>(codon.bases_), std::get<2>(codon.bases_));
	};
	//key1 = transcript, key2 = amino acid position, value = the codon's genomic location
	std::unordered_map<std::string, std::unordered_map<uint32_t, Bed6RecordCore>> codonLocations;

	for (auto & mhap : repSearch->microhaplotypes_) {
		mhap.associated_protein_variants_.reset();
		mhap.associated_seq_variants_.reset();
		auto nameSearch = callNameBySeq.find(mhap.seq_);
		const bool hasProtein = callNameBySeq.end() != nameSearch && njh::in(nameSearch->second, calls.fullAATypedWithCodonInfo_);
		const bool hasSeq = callNameBySeq.end() != nameSearch && njh::in(nameSearch->second, calls.seqAlns_);
		if (!hasProtein && !hasSeq) {
			++ret.mhapsWithoutCalls_;
			continue;
		}
		++ret.mhapsWithCalls_;
		if (hasProtein) {
			std::vector<pmo::ProteinVariant> proteinVariants;
			for (const auto & aa : calls.fullAATypedWithCodonInfo_.at(nameSearch->second)) {
				pmo::ProteinVariant variant;
				//amino acid positions are 1-based with the end included, start and end are the same for one amino acid
				variant.protein_location_ = genLocation(aa.transcriptName_, aa.zeroBasedPos_ + 1, aa.zeroBasedPos_ + 1,
						std::string(1, aa.refCod_.aa_), std::string(1, aa.cod_.aa_), std::nullopt);
				auto geneSearch = geneByTranscript.find(aa.transcriptName_);
				if (geneByTranscript.end() != geneSearch) {
					variant.gene_name_ = geneSearch->second;
				}
				auto altNameSearch = altNameByTranscript.find(aa.transcriptName_);
				if (altNameByTranscript.end() != altNameSearch) {
					variant.alternative_gene_name_ = altNameSearch->second;
				}
				auto & codonLocs = codonLocations[aa.transcriptName_];
				auto codonSearch = codonLocs.find(aa.zeroBasedPos_);
				if (codonLocs.end() == codonSearch) {
					codonSearch = codonLocs.emplace(aa.zeroBasedPos_,
							calls.translationInfoForTranscirpt_.at(aa.transcriptName_)->genBedFromAAPositions(aa.zeroBasedPos_, aa.zeroBasedPos_ + 1)).first;
				}
				const auto & codonLoc = codonSearch->second;
				variant.codon_genomic_location_ = genLocation(codonLoc.chrom_, codonLoc.chromStart_, codonLoc.chromEnd_,
						codonStr(aa.refCod_), codonStr(aa.cod_), codonLoc.strand_);
				proteinVariants.emplace_back(variant);
			}
			ret.proteinVariants_ += proteinVariants.size();
			if (!proteinVariants.empty()) {
				mhap.associated_protein_variants_ = proteinVariants;
			}
		}
		if (hasSeq) {
			std::vector<pmo::GenomicLocation> seqVariants;
			for (const auto & aln : calls.seqAlns_.at(nameSearch->second)) {
				const auto gPos = aln.gRegion_.genBed3RecordCore();
				for (const auto & m : aln.comp_.distances_.mismatches_) {
					const uint32_t pos = gPos.chromStart_ + m.second.refBasePos;
					seqVariants.emplace_back(genLocation(gPos.chrom_, pos, pos + 1,
							std::string(1, m.second.refBase), std::string(1, m.second.seqBase), '+'));
				}
				for (const auto & indel : aln.comp_.distances_.alignmentGaps_) {
					const uint32_t pos = gPos.chromStart_ + indel.second.refPos_;
					const std::string gaps(indel.second.gapedSequence_.size(), '-');
					if (indel.second.ref_) {
						//gap in the reference, an insertion
						seqVariants.emplace_back(genLocation(gPos.chrom_, pos, pos + 1, gaps, indel.second.gapedSequence_, '+'));
					} else {
						seqVariants.emplace_back(genLocation(gPos.chrom_, pos, pos + indel.second.gapedSequence_.size(), indel.second.gapedSequence_, gaps, '+'));
					}
				}
			}
			std::stable_sort(seqVariants.begin(), seqVariants.end(), [](const pmo::GenomicLocation & a, const pmo::GenomicLocation & b) {
				return std::tie(a.chrom_, a.start_, a.end_) < std::tie(b.chrom_, b.start_, b.end_);
			});
			ret.seqVariants_ += seqVariants.size();
			if (!seqVariants.empty()) {
				mhap.associated_seq_variants_ = seqVariants;
			}
		}
	}
	genomeSearch->validate();
	repSearch->validate();
	return ret;
}

}  // namespace njhseq
