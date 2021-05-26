/*
 * TranslatorByAlignment.cpp
 *
 *  Created on: Aug 14, 2019
 *      Author: nicholashathaway
 */


#include "TranslatorByAlignment.hpp"


namespace njhseq {



TranslatorByAlignment::TranslatorByAlignmentPars::TranslatorByAlignmentPars(){
	lzPars_.coverage = 100;
	lzPars_.identity = 70;
}


void TranslatorByAlignment::TranslatorByAlignmentPars::setOptions(seqSetUp & setUp){
	setUp.setOption(knownAminoAcidMutationsFnp_, "--knownAminoAcidChangesFnp",
			"Known Amino Acid Changes, must have at least 2 columns, positions are 1-postion-based (first position is 1), 1)TranscriptID, 2)AAPosition ", false, "Translation Output");
	setUp.setOption(gffFnp_, "--gffFnp",
			"Gff file to intersect the final haplotypes with genes to get translations", "" != knownAminoAcidMutationsFnp_, "Translation Output");
	setUp.setOption(lzPars_.genomeFnp, "--genomeFnp",
			"Genome file so final haplotypes can be mapped to a genome", "" != gffFnp_ || "" != knownAminoAcidMutationsFnp_, "Translation Output");
}




std::tuple<char, char, char> TranslatorByAlignment::TranslateSeqRes::getCodonForAARefPos(
		uint32_t aaRefPos) const {
	if (aaRefPos < std::get<0>(firstAminoInfo_).aaPos_
			|| aaRefPos > std::get<0>(lastAminoInfo_).aaPos_) {
		return std::make_tuple('-', '-', '-');
	}

	uint32_t alnPosForLoc = getAlnPosForRealPos(refAlnTranslation_.seq_,
			aaRefPos);
	uint32_t queryCDNAPos = getRealPosForAlnPos(queryAlnTranslation_.seq_,
			alnPosForLoc) * 3;

	if (queryCDNAPos + 2 >= cDna_.seq_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "" << queryCDNAPos + 2
				<< " is greater than cDna_.seq_.size(): "
				<< cDna_.seq_.size() << "\n";
		throw std::runtime_error { ss.str() };
	}
	return std::make_tuple(cDna_.seq_[queryCDNAPos + 0], cDna_.seq_[queryCDNAPos + 1],
			cDna_.seq_[queryCDNAPos + 2]);

}




TranslatorByAlignment::VariantsInfo::VariantsInfo(const Bed3RecordCore & region, const seqInfo & refSeq) : region_(region),
		seqBase_(refSeq){

}


char TranslatorByAlignment::VariantsInfo::getBaseForGenomicRegionNoCheck(const uint32_t pos) const{
	return seqBase_.seq_[pos - region_.chromStart_];
}

char TranslatorByAlignment::VariantsInfo::getBaseForGenomicRegion(const uint32_t pos) const{
	if(pos < region_.chromStart_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "position: " << pos << " is before the start of the region: " << region_.chromStart_<< "\n";
		throw std::runtime_error{ss.str()};
	}
	uint32_t relativePos = pos -region_.chromStart_;
	if(relativePos > len(seqBase_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "position: " << pos << ", relative position: " << relativePos << " is out of range, seq size is: " << len(seqBase_)<< "\n";
		throw std::runtime_error{ss.str()};
	}
	return seqBase_.seq_[relativePos];
}
void TranslatorByAlignment::VariantsInfo::writeVCF(const OutOptions & vcfOutOpts) const {
	OutputStream vcfOut(vcfOutOpts);

	std::unordered_set<uint32_t> positionsSet;
	for(const auto & snps : snpsFinal){
		positionsSet.emplace(snps.first);
	}


	std::map<uint32_t, std::map<std::string,uint32_t>> insertionsFinalForVCF;
	std::map<uint32_t, std::map<std::string,uint32_t>> deletionsFinalForVCF;
	for(const auto & ins : insertionsFinal){
		if(0 == ins.first ){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle insertion at position 0"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		insertionsFinalForVCF[ins.first - 1] = ins.second;
	}
	for(const auto & del : deletionsFinal){
		if(0 == del.first ){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle insertion at position 0"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		deletionsFinalForVCF[del.first - 1] = del.second;
	}

	for(const auto & ins : insertionsFinalForVCF){
		positionsSet.emplace(ins.first);
	}
	for(const auto & del : deletionsFinalForVCF){
		positionsSet.emplace(del.first);
	}


	vcfOut << "##fileformat=VCFv4.0" << std::endl;
	vcfOut << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth\">" << std::endl;
	vcfOut << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
	vcfOut << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">" << std::endl;
	vcfOut << "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele Count\">" << std::endl;
	vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

	std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());

	njh::sort(positions);
	for(const auto & pos : positions){
		if (njh::in(pos, insertionsFinalForVCF) || njh::in(pos, snpsFinal)) {
			vcfOut <<  region_.chrom_
					<< "\t" << pos+ 1
					<< "\t" << "."
					<< "\t";
			std::vector<std::string> alts;
			std::vector<uint32_t> altsCounts;
			std::vector<double> altsFreqs;
			vcfOut << getBaseForGenomicRegion(pos) << "\t";
			if(njh::in(pos, snpsFinal)){
				uint32_t snpCount = 0;
				for(const auto & b : snpsFinal.at(pos)){
					snpCount+= b.second;
					alts.emplace_back(std::string(1, b.first));
					altsCounts.emplace_back(b.second);
					altsFreqs.emplace_back(b.second/static_cast<double>(depthPerPosition.at(pos)));
				}
			}
			if (njh::in(pos, insertionsFinalForVCF)) {
				for (const auto & ins : insertionsFinalForVCF[pos]) {
					alts.emplace_back(njh::pasteAsStr(getBaseForGenomicRegion(pos), ins.first));
					altsCounts.emplace_back(ins.second);
					altsFreqs.emplace_back(ins.second/static_cast<double>(depthPerPosition.at(pos)));
				}
			}

			vcfOut << njh::conToStr(alts, ",")
			<< "\t40\tPASS\t";
			vcfOut
					<< "DP=" << depthPerPosition.at(pos) << ";"
					<< "NS=" << samplesPerPosition.at(pos).size() << ";"
					<< "AC=" << njh::conToStr(altsCounts, ",") << ";"
					<< "AF=" << njh::conToStr(altsFreqs, ",")
			<< std::endl;
		}
		if (njh::in(pos, deletionsFinalForVCF)) {
			for (const auto & d : deletionsFinalForVCF[pos]) {
				vcfOut <<  region_.chrom_
						<< "\t" << pos + 1
						<< "\t" << "."
						<< "\t";
				vcfOut << getBaseForGenomicRegion(pos) << d.first
				<< "\t" << getBaseForGenomicRegion(pos) << "\t";
				vcfOut << "40\tPASS\t";
				vcfOut
						<< "DP=" << depthPerPosition.at(pos) << ";"
						<< "NS=" << samplesPerPosition.at(pos).size() << ";"
						<< "AC=" << d.second << ";"
						<< "AF=" << d.second/static_cast<double>(depthPerPosition.at(pos))
				<< std::endl;
			}
		}
	}
}


void TranslatorByAlignment::VariantsInfo::writeSNPTable(const OutOptions &snpTabOutOpts)const{


	OutputStream snpTabOut(snpTabOutOpts);


	std::unordered_set<uint32_t> positionsSet;
	for(const auto & snps : snpsFinal){
		positionsSet.emplace(snps.first);
	}


	std::map<uint32_t, std::map<std::string,uint32_t>> insertionsFinalForVCF;
	std::map<uint32_t, std::map<std::string,uint32_t>> deletionsFinalForVCF;
	for(const auto & ins : insertionsFinal){
		if(0 == ins.first ){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle insertion at position 0"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		insertionsFinalForVCF[ins.first - 1] = ins.second;
	}
	for(const auto & del : deletionsFinal){
		if(0 == del.first ){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle insertion at position 0"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		deletionsFinalForVCF[del.first - 1] = del.second;
	}

	for(const auto & ins : insertionsFinalForVCF){
		positionsSet.emplace(ins.first);
	}
	for(const auto & del : deletionsFinalForVCF){
		positionsSet.emplace(del.first);
	}

	snpTabOut << "chrom\tposition\tref\tvariant\tcount\tfrequency\talleleDepth\tsamples" << std::endl;

	std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());

	njh::sort(positions);
	for(const auto & pos : positions){
		if (njh::in(pos, insertionsFinalForVCF) || njh::in(pos, snpsFinal)) {
			std::vector<std::string> alts;
			std::vector<uint32_t> altsCounts;
			std::vector<double> altsFreqs;


			if(njh::in(pos, snpsFinal)){
				uint32_t snpCount = 0;
				for(const auto & b : snpsFinal.at(pos)){
					snpTabOut << region_.chrom_
							<< "\t" << pos
							<< "\t" << getBaseForGenomicRegion(pos)
							<< "\t" << std::string(1, b.first)
							<< "\t" << b.second
							<< "\t" << b.second/static_cast<double>(depthPerPosition.at(pos))
							<< "\t" << depthPerPosition.at(pos)
							<< "\t" << samplesPerPosition.at(pos).size() << std::endl;
					snpCount+= b.second;
					alts.emplace_back(std::string(1, b.first));
					altsCounts.emplace_back(b.second);
					altsFreqs.emplace_back(b.second/static_cast<double>(depthPerPosition.at(pos)));
				}
				snpTabOut << region_.chrom_
						<< "\t" << pos
						<< "\t" << getBaseForGenomicRegion(pos)
						<< "\t" << getBaseForGenomicRegion(pos)
						<< "\t" << depthPerPosition.at(pos) - snpCount
						<< "\t" << (depthPerPosition.at(pos) - snpCount)/static_cast<double>(depthPerPosition.at(pos))
						<< "\t" << depthPerPosition.at(pos)
						<< "\t" << samplesPerPosition.at(pos).size() << std::endl;
			}
		}
	}
}


void TranslatorByAlignment::VariantsInfo::setFinals(const RunPars & rPars){
//	std::cout << "rPars.lowVariantCutOff: " << rPars.lowVariantCutOff << std::endl;
//	std::cout << "rPars.occurrenceCutOff: " << rPars.occurrenceCutOff << std::endl;
//	std::cout << "totalPopCount: " << totalPopCount << std::endl;

	depthPerPosition.clear();
	snpsFinal.clear();
	deletionsFinal.clear();
	insertionsFinal.clear();
	variablePositons_.clear();

	//get counts per position
	for(const auto & pos : allBases){
		for(const auto & base : pos.second){
			//std::cout << "pos:" << pos.first << ",base:" << base.first << ",count:" << base.second << std::endl;
			depthPerPosition[pos.first] += base.second; //this add both deletion and match/mismatch because the addVariant calls adds - in query counts
		}
	}
//	std::cout << njh::bashCT::red;
//	for(const auto & pos : depthPerPosition){
//		std::cout << "pos: " << pos.first << ":depth: " << pos.second << std::endl;
//	}
//	std::cout << njh::bashCT::reset;

	//add depth counts for deleted sections, this is not necessary since the above addition of allBases will add deltions
//	for(const auto & pos : deletions){
//		for(const auto & del : pos.second){
//			for(const auto seqPos : iter::range(del.first.size())){
//				depthPerPosition[pos.first + seqPos] += del.second;
//			}
//		}
//	}
//	std::cout << njh::bashCT::blue;
//	for(const auto & pos : depthPerPosition){
//		std::cout << "pos: " << pos.first << ":depth: " << pos.second << std::endl;
//	}
//	std::cout << njh::bashCT::reset;

	//filter saps and indels by occurrence cut off
	for(auto & snp : snps){
		for(const auto & b : snp.second){
			//if(b.second < rPars.occurrenceCutOff || b.second/static_cast<double>(totalPopCount) < rPars.lowVariantCutOff){
			if(b.second < rPars.occurrenceCutOff || b.second/static_cast<double>(depthPerPosition[snp.first]) < rPars.lowVariantCutOff){
			continue;
			}
			snpsFinal[snp.first][b.first] = b.second;
			variablePositons_.emplace(snp.first);
		}
	}
	for(const auto & del : deletions){
		for(const auto & d : del.second){
			//if(d.second < rPars.occurrenceCutOff || d.second/static_cast<double>(totalPopCount) < rPars.lowVariantCutOff){
			if(d.second < rPars.occurrenceCutOff || d.second/static_cast<double>(depthPerPosition[del.first]) < rPars.lowVariantCutOff){
				continue;
			}

			deletionsFinal[del.first][d.first] = d.second;
			variablePositons_.emplace(del.first);
		}
	}
	for(const auto & ins : insertions){
		for(const auto & i : ins.second){
//		if(i.second < rPars.occurrenceCutOff || i.second/static_cast<double>(totalPopCount) < rPars.lowVariantCutOff){
			if(i.second < rPars.occurrenceCutOff || i.second/static_cast<double>(depthPerPosition[ins.first]) < rPars.lowVariantCutOff){
				continue;
			}
			insertionsFinal[ins.first][i.first] = i.second;
			variablePositons_.emplace(ins.first);
		}
	}
}


Bed3RecordCore TranslatorByAlignment::VariantsInfo::getVariableRegion() {
	if (!variablePositons_.empty()) {
		return Bed3RecordCore(seqBase_.name_,
				*std::min_element(variablePositons_.begin(), variablePositons_.end()),
				*std::max_element(variablePositons_.begin(), variablePositons_.end()) + 1 );
	}
	return Bed3RecordCore(seqBase_.name_, std::numeric_limits<uint32_t>::max(),
			std::numeric_limits<uint32_t>::max());
}


void TranslatorByAlignment::VariantsInfo::addVariantInfo(
		const std::string & alignedRefSeq,
		const std::string & alignedQuerySeq,
		uint32_t querySeqCount,
		const std::unordered_set<std::string> & samples,
		const comparison & comp,
		uint32_t offSetStart
		){
	uint32_t queryAlnStart = alignedQuerySeq.find_first_not_of("-");
	uint32_t queryAlnEnd = alignedQuerySeq.find_last_not_of("-");
	for(const auto seqPos : iter::range(queryAlnStart, queryAlnEnd + 1)){
		if('-' != alignedRefSeq[seqPos]){ //skip over insertions
			uint32_t seqChromPosition = getRealPosForAlnPos(alignedRefSeq, seqPos) + offSetStart;
			allBases[seqChromPosition][alignedQuerySeq[seqPos]] += querySeqCount;
			samplesPerPosition[seqChromPosition].insert(samples.begin(), samples.end());
		}
	}
	for(const auto & m : comp.distances_.mismatches_){
		snps[m.second.refBasePos + offSetStart][m.second.seqBase]+= querySeqCount;

	}
	for(const auto & gap : comp.distances_.alignmentGaps_){
		if(gap.second.ref_){
			//insertion
			insertions[gap.second.refPos_ + offSetStart][gap.second.gapedSequence_]+=querySeqCount;
		}else{
			//deletion
			deletions[gap.second.refPos_ + offSetStart][gap.second.gapedSequence_]+=querySeqCount;
			for(const auto pos : iter::range(gap.second.gapedSequence_.size())){
				samplesPerPosition[gap.second.refPos_ + offSetStart + pos].insert(samples.begin(), samples.end());
			}
		}
	}
}


std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> TranslatorByAlignment::translateBasedOnAlignment(
			const ReAlignedSeq & realigned,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			aligner & alignerObj){
	std::unordered_map<std::string, TranslateSeqRes> ret;
	for(const auto & transcript : currentGene.mRNAs_){
		auto currentTranscriptInfo = transcriptInfosForGene.at(transcript->getIDAttr());
		auto genePosInfoByGDna = currentTranscriptInfo->getInfosByGDNAPos();
		bool endsAtStopCodon = false;
		uint32_t transStart = 0;
		seqInfo balnSeq(realigned.querySeq_.name_);
		std::vector<uint32_t> codons;
		std::vector<GFFCore> cDNAIntersectedWith;
		for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
			if (realigned.gRegion_.overlaps(*cDna, 3)) {
				cDNAIntersectedWith.emplace_back(*cDna);
			}
		}
		if(cDNAIntersectedWith.size() == 0){

		} else {
			if (cDNAIntersectedWith.size() == 1
					&& realigned.gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
					&& realigned.gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
				balnSeq = realigned.querySeq_;
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(realigned.gRegion_.start_).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
					uint32_t gPos = realigned.gRegion_.end_ - 1;
					auto codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				} else {
					if (genePosInfoByGDna.at(realigned.gRegion_.end_ - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
					uint32_t gPos = realigned.gRegion_.start_;
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
					uint32_t gPos = std::min(cDnaStop, realigned.gRegion_.end_) - 1;
					auto codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				} else {
					auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
					uint32_t gPos = std::max(cDnaStart, realigned.gRegion_.start_);
					uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
						++transStart;
					}
				}
				std::vector<uint32_t> starts;
				std::vector<uint32_t> ends;
				for (const auto & cDna : cDNAIntersectedWith) {
					auto cDnaStart = cDna.start_ - 1;
					auto detStart = std::max(cDnaStart, realigned.gRegion_.start_);
					auto detStop = std::min(cDna.end_, realigned.gRegion_.end_);
					ends.emplace_back(detStop);
					starts.emplace_back(detStart);
					detStart -= realigned.gRegion_.start_;
					detStop -= realigned.gRegion_.start_;
					auto alnStart = getAlnPosForRealPos(realigned.alnRefSeq_.seq_,detStart);
					auto alnStop = getAlnPosForRealPos(realigned.alnRefSeq_.seq_, detStop - 1);
					balnSeq.append(
							realigned.alnQuerySeq_.getSubRead(alnStart,
									alnStop - alnStart + 1));
				}
				uint32_t cDnaStart = *std::min_element(starts.begin(), starts.end());
				uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(cDnaStart).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				} else {
					if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				}
				balnSeq.removeGaps();
			}
			if (currentGene.gene_->isReverseStrand()) {
				balnSeq.reverseComplementRead(false, true);
			}
			auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);
			MetaDataInName transMeta;
			transMeta.addMeta("transcript", transcript->getIDAttr());
			balnSeqTrans.name_ += transMeta.createMetaName();

			alignerObj.alignCacheGlobal(currentTranscriptInfo->protein_, balnSeqTrans);
			alignerObj.profilePrimerAlignment(currentTranscriptInfo->protein_, balnSeqTrans);
			TranslateSeqRes tRes;

			uint32_t cDnaLenRaw = len(balnSeq) - transStart;
			uint32_t cDnaLen = cDnaLenRaw - (cDnaLenRaw %3);
			uint32_t firstAmino = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-"));
			uint32_t lastAmino = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-"));
			lastAmino = std::min<uint32_t>(lastAmino,len(currentTranscriptInfo->protein_) - 1);
			auto aminoInfos = currentTranscriptInfo->getInfosByAAPos();
			tRes.firstAminoInfo_ = njh::mapAt(aminoInfos, firstAmino);
			tRes.lastAminoInfo_ = njh::mapAt(aminoInfos, lastAmino);
			tRes.cDna_ = balnSeq.getSubRead(transStart, cDnaLen);
			tRes.transcriptName_ = transcript->getIDAttr();
			tRes.translation_ = balnSeqTrans;
			tRes.refAlnTranslation_ = alignerObj.alignObjectA_.seqBase_;
			tRes.queryAlnTranslation_ = alignerObj.alignObjectB_.seqBase_;
			tRes.comp_ = alignerObj.comp_;

			ret[transcript->getIDAttr()] = tRes;
		}
	}
	return ret;
}


std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> TranslatorByAlignment::translateBasedOnAlignment(
		const BamTools::BamAlignment & bAln,
		const GeneFromGffs & currentGene,
		const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
		TwoBit::TwoBitFile & tReader,
		aligner & alignerObj,
		const BamTools::RefVector & refData){

	std::unordered_map<std::string, TranslateSeqRes> ret;

	auto results = std::make_shared<AlignmentResults>(bAln, refData, true);
	results->setRefSeq(tReader);
	results->setComparison(true);


	for(const auto & transcript : currentGene.mRNAs_){
		auto currentTranscriptInfo = transcriptInfosForGene.at(transcript->getIDAttr());
		auto genePosInfoByGDna = currentTranscriptInfo->getInfosByGDNAPos();
		bool endsAtStopCodon = false;
		uint32_t transStart = 0;
		seqInfo balnSeq(bAln.Name);
		std::vector<uint32_t> codons;
		std::vector<GFFCore> cDNAIntersectedWith;
		for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
			if (results->gRegion_.overlaps(*cDna)) {
				cDNAIntersectedWith.emplace_back(*cDna);
			}
		}
		if(cDNAIntersectedWith.size() == 0){

		} else {
			if (cDNAIntersectedWith.size() == 1
					&& results->gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
					&& results->gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
				balnSeq = *(results->alnSeq_);
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(results->gRegion_.start_).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
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
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
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
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				} else {
					if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				}
				balnSeq.removeGaps();
			}
			if (currentGene.gene_->isReverseStrand()) {
				balnSeq.reverseComplementRead(false, true);
			}

			auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);


			MetaDataInName transMeta;
			transMeta.addMeta("transcript", transcript->getIDAttr());
			balnSeqTrans.name_ += transMeta.createMetaName();

			alignerObj.alignCacheGlobal(currentTranscriptInfo->protein_, balnSeqTrans);
			alignerObj.profilePrimerAlignment(currentTranscriptInfo->protein_, balnSeqTrans);
			TranslateSeqRes tRes;

			uint32_t cDnaLenRaw = len(balnSeq) - transStart;
			uint32_t cDnaLen = cDnaLenRaw - (cDnaLenRaw %3);
			uint32_t firstAmino = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of("-"));
			uint32_t lastAmino = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of("-"));
			lastAmino = std::min<uint32_t>(lastAmino,len(currentTranscriptInfo->protein_) - 1);

			auto aminoInfos = currentTranscriptInfo->getInfosByAAPos();
			tRes.firstAminoInfo_ = njh::mapAt(aminoInfos, firstAmino);
			tRes.lastAminoInfo_ = njh::mapAt(aminoInfos, lastAmino);
			tRes.cDna_ = balnSeq.getSubRead(transStart, cDnaLen);
			tRes.transcriptName_ = transcript->getIDAttr();
			tRes.translation_ = balnSeqTrans;
			tRes.refAlnTranslation_ = alignerObj.alignObjectA_.seqBase_;
			tRes.queryAlnTranslation_ = alignerObj.alignObjectB_.seqBase_;
			tRes.comp_ = alignerObj.comp_;

			ret[transcript->getIDAttr()] = tRes;
		}
	}
	return ret;
}


TranslatorByAlignment::TranslatorByAlignment(const TranslatorByAlignmentPars & pars): pars_(pars){
	njh::sys::requireExternalProgramThrow("samtools");
	if(!pars_.useLastz_){
		njh::sys::requireExternalProgramThrow("bowtie2");
	}else{
		njh::sys::requireExternalProgramThrow("lastz");
	}
	if("" != pars_.knownAminoAcidMutationsFnp_){
		table knownAminoAcidChanges(pars_.knownAminoAcidMutationsFnp_, "\t", true);
		VecStr originalColNames = knownAminoAcidChanges.columnNames_;
		njh::for_each(knownAminoAcidChanges.columnNames_, [](std::string & col){
			njh::strToLower(col);
		});
		knownAminoAcidChanges.setColNamePositions();
		knownAminoAcidChanges.checkForColumnsThrow(VecStr{"transcriptid", "aaposition"}, __PRETTY_FUNCTION__);
		if(knownAminoAcidChanges.nRow() > 0){
			VecStr warnings;
			for(const auto & row : knownAminoAcidChanges){
				if(std::all_of(row.begin(), row.end(), [](const std::string & element){
					return "" == element;
				})){
					continue;
				}
				auto transciprtId = row[knownAminoAcidChanges.getColPos("transcriptid")];
				auto aaPosition = njh::StrToNumConverter::stoToNum<uint32_t>(row[knownAminoAcidChanges.getColPos("aaposition")]);
				if(njh::in(aaPosition, knownAminoAcidPositions_[transciprtId])){
					warnings.emplace_back(njh::pasteAsStr("already have aaposition ", aaPosition, " for transcriptID: ", transciprtId));
				}
				knownAminoAcidPositions_[transciprtId].emplace_back(aaPosition);
				if(originalColNames.size() > 2){
					for(const auto && colPos : iter::range(knownAminoAcidChanges.columnNames_.size())){
						if(!njh::in(knownAminoAcidChanges.columnNames_[colPos], VecStr{"transcriptid", "aaposition"})){
							metaDataAssociatedWithAminoacidPosition_[transciprtId][aaPosition].addMeta(originalColNames[colPos], row[colPos]);
						}
					}
				}
			}
			if(!warnings.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "found repeat amino acid positions "<< "\n";
				for(const auto & warn : warnings){
					ss << warn << std::endl;
				}
				throw std::runtime_error{ss.str()};
			}
		}
	}
}


TranslatorByAlignment::TranslatorByAlignmentResult TranslatorByAlignment::run(
		const SeqIOOptions & seqOpts,
					const std::unordered_map<std::string, std::unordered_set<std::string>> & sampCountsForHaps,
					const RunPars & rPars){

	uint32_t totalPopCount = 0;
	for(const auto & sampCount : sampCountsForHaps){
		totalPopCount += sampCount.second.size();
	}
	TranslatorByAlignmentResult ret;
	std::vector<bfs::path> fnpsToRemove;
	auto seqInputFnp = njh::files::make_path(pars_.workingDirtory_, "inputSeqs.fasta");
	fnpsToRemove.emplace_back(seqInputFnp);
	uint64_t seqMaxLen = 0;

	VecStr names;
	{
		//write out fasta file of input
		seqInfo seq;
		SeqInput reader(seqOpts);
		reader.openIn();
		SeqOutput writer(SeqIOOptions::genFastaOut(seqInputFnp));
		writer.openOut();
		uint32_t pos = 0;
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, seqMaxLen);
			names.emplace_back(seq.name_);
			seq.name_ = estd::to_string(pos);
			++pos;
			writer.write(seq);
		}
	}

//	for(const auto & pos : iter::range(names.size())){
//		std::cout << pos << ":" << names[pos] << std::endl;
//	}

	BioCmdsUtils bRunner(false);
	bRunner.RunFaToTwoBit(pars_.lzPars_.genomeFnp);
	if(!pars_.useLastz_){
		bRunner.RunBowtie2Index(pars_.lzPars_.genomeFnp);
	}

	auto gprefix = bfs::path(pars_.lzPars_.genomeFnp).replace_extension("");
	auto twoBitFnp = gprefix.string() + ".2bit";

	TwoBit::TwoBitFile tReader(twoBitFnp);

	auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(seqInputFnp);
	uniqueSeqInOpts.out_.outFilename_ = njh::files::make_path(pars_.workingDirtory_, "aligned_inputSeqs.sorted.bam");
	uniqueSeqInOpts.out_.outExtention_ = ".sorted.bam";
	fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_);
	fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_.string() + ".bai");

	uniqueSeqInOpts.out_.transferOverwriteOpts(seqOpts.out_);
	if(!pars_.useLastz_){
		auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, pars_.additionalBowtieArguments_);

		//auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 --end-to-end");
		BioCmdsUtils::checkRunOutThrow(bowtieRunOut, __PRETTY_FUNCTION__);
	}else{
		auto lastzRunOut = bRunner.lastzAlign(uniqueSeqInOpts, pars_.lzPars_);
		BioCmdsUtils::checkRunOutThrow(lastzRunOut, __PRETTY_FUNCTION__);
	}

	auto regionsCounter = GenomicRegionCounter::countRegionsInBam(uniqueSeqInOpts.out_.outName());
	auto ids = regionsCounter.getIntersectingGffIds(pars_.gffFnp_);
	ret.geneIds_ = ids;

	// get gene information
	auto geneInfoDir = njh::files::make_path(pars_.workingDirtory_, "geneInfos");
	if(pars_.keepTemporaryFiles_){
		njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
	}
	OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));

	std::unordered_map<std::string, VecStr> idToTranscriptName;
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes = GeneFromGffs::getGenesFromGffForIds(pars_.gffFnp_, ids);
	//std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>> geneTranscriptInfos;
	;
	uint64_t proteinMaxLen = 0;
	std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;

	for(const auto & gene : genes){
		genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
		for(const auto & transcript : gene.second->mRNAs_){
			idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
		}
		ret.transcriptInfosForGene_[gene.first] = gene.second->generateGeneSeqInfo(tReader, false);
		for(const auto & transcriptInfo : ret.transcriptInfosForGene_[gene.first]){
			ret.translationInfoForTranscirpt_[transcriptInfo.first] = transcriptInfo.second;
			readVec::getMaxLength(transcriptInfo.second->protein_, proteinMaxLen);
		}
		if(pars_.keepTemporaryFiles_){
			gene.second->writeOutGeneInfo(tReader, outOpts);
		}
	}

	aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2));
	aligner alignObjSeq(seqMaxLen + rPars.realnPars.extendAmount * 2, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2));

	std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> regionsToGeneIds;
	//targetName, GeneID, AA Position
	std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
	for(const auto & gCount : regionsCounter.counts_){
		for (const auto & g : genesByChrom[gCount.second.region_.chrom_]) {
			for(const auto & t : g.second->mRNAs_){
				regionsToGeneIds[gCount.second.region_.createUidFromCoords()][g.first].emplace(t->getIDAttr());
				alnRegionToGeneIds[gCount.second.region_.createUidFromCoords()].emplace_back(g.first);
			}
		}
	}
	BamTools::BamReader bReader;
	bReader.Open(uniqueSeqInOpts.out_.outName().string());
	checkBamOpenThrow(bReader, uniqueSeqInOpts.out_.outName());
	auto refData = bReader.GetReferenceData();
	BamTools::BamAlignment bAln;
	auto chromLengths = tReader.getSeqLens();

	struct MinMaxPos{
		MinMaxPos(){

		}
		size_t minPos_{std::numeric_limits<uint32_t>::max()};
		size_t maxPos_{0};
	};
	std::unordered_map<std::string, MinMaxPos> minMaxPositionsPerChrom;
	while (bReader.GetNextAlignment(bAln)) {
		if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
			bAln.Name = names[njh::StrToNumConverter::stoToNum<uint32_t>(bAln.Name)];
			auto balnGenomicRegion = GenomicRegion(bAln, refData);
			minMaxPositionsPerChrom[balnGenomicRegion.chrom_].minPos_ = std::min(minMaxPositionsPerChrom[balnGenomicRegion.chrom_].minPos_,balnGenomicRegion.start_);
			minMaxPositionsPerChrom[balnGenomicRegion.chrom_].maxPos_ = std::max(minMaxPositionsPerChrom[balnGenomicRegion.chrom_].minPos_,balnGenomicRegion.end_);
			auto results = ReAlignedSeq::genRealignment(bAln, refData, alignObjSeq, chromLengths, tReader, rPars.realnPars);
			ret.seqAlns_[bAln.Name].emplace_back(results);
			if (!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)) {
				continue;
			}
			for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
				const auto & currentGene = genes.at(g);
				const auto & currentGeneInfo = ret.transcriptInfosForGene_.at(g);
				auto translations = translateBasedOnAlignment(results, *currentGene, currentGeneInfo, alignObj);
				for(const auto & trans : translations){
					ret.translations_[bAln.Name].emplace(trans);
				}
			}
		}
	}
	if(!pars_.keepTemporaryFiles_){
		for(const auto & fnp : fnpsToRemove){
			if(bfs::exists(fnp)){
				bfs::remove(fnp);
			}
		}
	}
	//index snps
	for(const auto & seqName : ret.seqAlns_){
		for(const auto & aln : seqName.second){
			if(!njh::in(aln.gRegion_.chrom_, ret.seqVariants_)){
				Bed3RecordCore chromRegion(
										aln.gRegion_.chrom_,
										minMaxPositionsPerChrom[aln.gRegion_.chrom_].minPos_,
										minMaxPositionsPerChrom[aln.gRegion_.chrom_].maxPos_);
				ret.seqVariants_.emplace(aln.gRegion_.chrom_, VariantsInfo(chromRegion, GenomicRegion(chromRegion).extractSeq(tReader)));
			}
			for(uint32_t seqPos = 0; seqPos < aln.refSeq_.seq_.size(); ++ seqPos){
				ret.baseForPosition_[aln.gRegion_.chrom_][aln.gRegion_.start_ + seqPos] = aln.refSeq_.seq_[seqPos];
			}
			uint32_t popCount = njh::mapAt(sampCountsForHaps, aln.querySeq_.name_).size();
			ret.seqVariants_.at(aln.gRegion_.chrom_).addVariantInfo(
					aln.alnRefSeq_.seq_,
					aln.alnQuerySeq_.seq_,
					popCount,
					sampCountsForHaps.at(aln.querySeq_.name_),
					aln.comp_,
					aln.gRegion_.start_);
		}
	}

	for(auto & varPerChrom : ret.seqVariants_){
		//varPerChrom.second.setFinals(rPars, totalPopCount);
		varPerChrom.second.setFinals(rPars);
	}

	//index amino acid changes per transcript
	for(const auto & seqName : ret.translations_){
		for(const auto & transcript : seqName.second){
			if("" == ret.proteinForTranscript_[transcript.first]){
				ret.proteinForTranscript_[transcript.first] = njh::replaceString(transcript.second.refAlnTranslation_.seq_, "-", "");
			}
			auto popCount = sampCountsForHaps.at(seqName.first).size();
			if(!njh::in(transcript.first, ret.proteinVariants_)){
				ret.proteinVariants_.emplace(transcript.first, VariantsInfo{Bed3RecordCore(transcript.first, 0, transcript.second.refAlnTranslation_.seq_.size()),
					transcript.second.refAlnTranslation_});
			}
			ret.proteinVariants_.at(transcript.first).addVariantInfo(transcript.second.refAlnTranslation_.seq_,
					transcript.second.queryAlnTranslation_.seq_, popCount, sampCountsForHaps.at(seqName.first), transcript.second.comp_, 0);
		}
	}

	for(auto & varPerTrans : ret.proteinVariants_){
		varPerTrans.second.setFinals(rPars);
		//varPerTrans.second.setFinals(rPars, totalPopCount);
	}

	return ret;
}



}  // namespace njhseq
