/*
 * GenomicAminoAcidPositionTyper.cpp
 *
 *  Created on: Sep 29, 2018
 *      Author: nick
 */



#include "GenomicAminoAcidPositionTyper.hpp"
#include "njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo.h"
#include "njhseq/BamToolsUtils.h"


namespace njhseq {

GenomicAminoAcidPositionTyper::GeneAminoTyperInfo::GeneAminoTyperInfo(const std::string & geneId) :
		geneId_(geneId), altId_(geneId) {
}
GenomicAminoAcidPositionTyper::GeneAminoTyperInfo::GeneAminoTyperInfo(const std::string & geneId,
		const std::map<uint32_t, char> & aminos) :
		geneId_(geneId), altId_(geneId), aminos_(aminos) {
}

GenomicAminoAcidPositionTyper::GenomicAminoAcidPositionTyper(const bfs::path & proteinMutantTypingFnp,
		bool inputZeroBased) :
		proteinMutantTypingFnp_(proteinMutantTypingFnp),
		inputZeroBased_(inputZeroBased){
	proteinMutantTypingTab_ = table(proteinMutantTypingFnp, "\t", true);
	proteinMutantTypingTab_.checkForColumnsThrow(VecStr { "ID", "AAPosition" }, __PRETTY_FUNCTION__);

	for (const auto & row : proteinMutantTypingTab_) {
		uint32_t aaPos = inputZeroBased ?
				njh::StrToNumConverter::stoToNum<uint32_t>(
						row[proteinMutantTypingTab_.getColPos("AAPosition")]) :
				njh::StrToNumConverter::stoToNum<uint32_t>(
						row[proteinMutantTypingTab_.getColPos("AAPosition")]) - 1;
		if(njh::in(aaPos, aminoPositionsForTyping_[row[proteinMutantTypingTab_.getColPos("ID")]])){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in reading in "
					<< proteinMutantTypingFnp << " gene id " << row[proteinMutantTypingTab_.getColPos("ID")]
					<< " already has position: " << aaPos  << "\n";
			throw std::runtime_error { ss.str() };
		}
		aminoPositionsForTyping_[row[proteinMutantTypingTab_.getColPos("ID")]].emplace_back(
				aaPos);
		if(proteinMutantTypingTab_.containsColumn("Gene")){
			altNamesForIds_[row[proteinMutantTypingTab_.getColPos("ID")]].emplace(row[proteinMutantTypingTab_.getColPos("Gene")]);
		}
	}
	for (const auto & geneId : aminoPositionsForTyping_) {
		if (altNamesForIds_[geneId.first].size() > 1) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in reading in "
					<< proteinMutantTypingFnp << " gene id " << geneId.first
					<< " had more than one alternative ID in Gene column, found: "
					<< njh::conToStr(altNamesForIds_[geneId.first], ", ") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	for(auto & positions : aminoPositionsForTyping_){
		njh::sort(positions.second);
	}
}

GenomicAminoAcidPositionTyper::GenomicAminoAcidPositionTyper(
		const std::unordered_map<std::string, std::vector<uint32_t>> & aminoPositionsForTyping) :
		inputZeroBased_(true), aminoPositionsForTyping_(aminoPositionsForTyping) {
	proteinMutantTypingTab_ = table(VecStr{"ID", "AAPosition"});
	for(const auto & gene : aminoPositionsForTyping){
		for(const auto & pos : gene.second){
			proteinMutantTypingTab_.addRow(gene.first, pos);
		}
	}
}

std::set<std::string> GenomicAminoAcidPositionTyper::getGeneIds() const {
	auto idsVec = proteinMutantTypingTab_.getColumnLevels("ID");
	return std::set<std::string>(idsVec.begin(), idsVec.end());
}

std::map<uint32_t, char> GenomicAminoAcidPositionTyper::typeAlignment(
		const BamTools::BamAlignment & bAln,
		const GeneFromGffs & currentGene,
		const GeneSeqInfo & currentGeneInfo,
		TwoBit::TwoBitFile & tReader,
		aligner & alignerObj,
		const BamTools::RefVector & refData) {

	auto results = std::make_shared<AlignmentResults>(bAln, refData, true);
	std::map<uint32_t, char> aminoTyping;

	results->setRefSeq(tReader);
	results->setComparison(true);

	bool endsAtStopCodon = false;
	uint32_t transStart = 0;
	std::unordered_map<size_t, alnInfoLocal> balnAlnInfo = bamAlnToAlnInfoLocal(bAln);
	auto genePosInfoByGDna = currentGeneInfo.getInfosByGDNAPos();
	const auto & transcript = currentGene.mRNAs_.front();
	seqInfo balnSeq(bAln.Name);
	std::vector<GFFCore> cDNAIntersectedWith;
	for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
		if (results->gRegion_.overlaps(*cDna)) {
			cDNAIntersectedWith.emplace_back(*cDna);
		}
	}
	if (cDNAIntersectedWith.size() == 1
			&& results->gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
			&& results->gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
		balnSeq = *(results->alnSeq_);
		if (currentGene.gene_->isReverseStrand()) {
			if (genePosInfoByGDna.at(results->gRegion_.start_).cDNAPos_
					== currentGeneInfo.cDna_.seq_.size() - 1) {
				endsAtStopCodon = true;
			}
			uint32_t gPos = results->gRegion_.end_ - 1;
			auto codon = genePosInfoByGDna.at(gPos).codonPos_;
			while (0 != codon) {
				--gPos;
				codon = genePosInfoByGDna.at(gPos).codonPos_;
				++transStart;
			}
		} else {
			if (genePosInfoByGDna.at(results->gRegion_.end_ - 1).cDNAPos_
					== currentGeneInfo.cDna_.seq_.size() - 1) {
				endsAtStopCodon = true;
			}
			uint32_t gPos = results->gRegion_.start_;
			uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
			while (0 != codon) {
				++gPos;
				codon = genePosInfoByGDna.at(gPos).codonPos_;
				++transStart;
			}
		}
	} else {
		njh::sort(cDNAIntersectedWith,
				[](const GenomicRegion & reg1, const GenomicRegion & reg2) {
					if(reg1.start_ < reg2.start_) {
						return true;
					}
					return false;
				});

		if (currentGene.gene_->isReverseStrand()) {
			auto cDnaStop = cDNAIntersectedWith.back().end_;
			uint32_t gPos = std::min(cDnaStop, results->gRegion_.end_) - 1;
			auto codon = genePosInfoByGDna.at(gPos).codonPos_;
			while (0 != codon) {
				--gPos;
				codon = genePosInfoByGDna.at(gPos).codonPos_;
				++transStart;
			}
		} else {
			auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
			uint32_t gPos = std::max(cDnaStart, results->gRegion_.start_);
			uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
			while (0 != codon) {
				++gPos;
				codon = genePosInfoByGDna.at(gPos).codonPos_;
				++transStart;
			}
		}
		std::vector<uint32_t> starts;
		std::vector<uint32_t> ends;
		for (const auto & cDna : cDNAIntersectedWith) {
			auto cDnaStart = cDna.start_ - 1;
			auto detStart = std::max(cDnaStart, results->gRegion_.start_);
			auto detStop = std::min(cDna.end_, results->gRegion_.end_);
			ends.emplace_back(detStop);
			starts.emplace_back(detStart);
			detStart -= results->gRegion_.start_;
			detStop -= results->gRegion_.start_;
			auto alnStart = getAlnPosForRealPos(results->refSeqAligned_->seq_,
					detStart);
			auto alnStop = getAlnPosForRealPos(results->refSeqAligned_->seq_,
					detStop - 1);
			balnSeq.append(
					results->alnSeqAligned_->getSubRead(alnStart,
							alnStop - alnStart + 1));
		}
		uint32_t cDnaStart = *std::min_element(starts.begin(), starts.end());
		uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
		if (currentGene.gene_->isReverseStrand()) {
			if (genePosInfoByGDna.at(cDnaStart).cDNAPos_
					== currentGeneInfo.cDna_.seq_.size() - 1) {
				endsAtStopCodon = true;
			}
		} else {
			if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
					== currentGeneInfo.cDna_.seq_.size() - 1) {
				endsAtStopCodon = true;
			}
		}
		balnSeq.removeGaps();
	}
	if (currentGene.gene_->isReverseStrand()) {
		balnSeq.reverseComplementRead(false, true);
	}
	auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);

	alignerObj.alignCacheGlobal(currentGeneInfo.protein_, balnSeqTrans);
	alignerObj.profilePrimerAlignment(currentGeneInfo.protein_, balnSeqTrans);

	if (njh::in(currentGene.gene_->getIDAttr(), aminoPositionsForTyping_)) {
		uint32_t proteinAlnStart =
				alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
		uint32_t proteinAlnStop =
				alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
		uint32_t proteinStart = getRealPosForAlnPos(
				alignerObj.alignObjectA_.seqBase_.seq_, proteinAlnStart);
		uint32_t proteinStop = getRealPosForAlnPos(
				alignerObj.alignObjectA_.seqBase_.seq_, proteinAlnStop);
		for (const auto & pos : aminoPositionsForTyping_.at(currentGene.gene_->getIDAttr())) {
			if (pos < proteinStart || pos > proteinStop) {
				aminoTyping[pos] = ' ';
			} else {
				auto posAln = getAlnPosForRealPos(
						alignerObj.alignObjectA_.seqBase_.seq_, pos);
				aminoTyping[pos] = alignerObj.alignObjectB_.seqBase_.seq_[posAln];
			}
		}
	}
	return aminoTyping;
}



}  // namespace njhseq
