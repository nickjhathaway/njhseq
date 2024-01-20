/*
 * GeneSeqInfo.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: nick
 */


#include "GeneSeqInfo.hpp"



namespace njhseq {

GeneSeqInfo::GeneSeqInfo(const seqInfo & cDna, seqInfo  gDna,
		GeneSeqInfoPars pars) :
		cDna_(cDna), gDna_(std::move(gDna)), pars_(std::move(pars)), infoTab_(VecStr { "cDnaPos", "chrom",
				"gDnaPos", "aaPos", "codonPos", "base", "aminoAcid" }) {

	protein_ = cDna.translateRet(false, false, 0, true);

}

void GeneSeqInfo::setCDnaAln(const seqInfo & cDnaAln) {
	cDnaAln_ = std::make_shared<seqInfo>(cDnaAln);
}

void GeneSeqInfo::setCDnaAlnByAligning(aligner & alignerObj){
	alignerObj.alignCacheGlobal(gDna_, cDna_);
	//check to see if alignment went ok
	for (const auto pos : iter::range(len(alignerObj.alignObjectA_))) {
		if (alignerObj.alignObjectA_.seqBase_.seq_[pos] != '-'
				&& alignerObj.alignObjectB_.seqBase_.seq_[pos] != '-') {
			if (alignerObj.alignObjectA_.seqBase_.seq_[pos]
					!= alignerObj.alignObjectB_.seqBase_.seq_[pos]) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error at aln position " << pos
						<< " cDNA doesn't match gDNA, "
								"check alignment parameters or make introns lowercase if they aren't already to aid in proper alignment"
						<< "\n";
				throw std::runtime_error { ss.str() };
			}
		}
	}

	cDnaAln_ = std::make_shared<seqInfo>(alignerObj.alignObjectB_.seqBase_);
}



void GeneSeqInfo::setTable(){
	if(nullptr == cDnaAln_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error  cDnaAln_ not set yet"  << "\n";
		throw std::runtime_error { ss.str() };
	}
	uint32_t posOffset = 0;
	if (pars_.oneBasedPos_) {
		posOffset = 1;
	}
	for (const auto pos : iter::range(len(cDna_))) {
		if(pars_.region_.reverseSrand_){
			infoTab_.addRow(pos + posOffset,
					pars_.region_.chrom_,
					pars_.region_.end_ - 1 - getAlnPosForRealPos(cDnaAln_->seq_, pos) + posOffset,
					(pos / 3) + posOffset,
					(pos % 3) + posOffset,
					cDna_.seq_[pos],
					len(protein_.seq_) > pos / 3 ? std::string(1, protein_.seq_[pos / 3]) : "NA");
		}else{
			infoTab_.addRow(pos + posOffset,
					pars_.region_.chrom_,
					getAlnPosForRealPos(cDnaAln_->seq_, pos) + pars_.region_.start_ + posOffset,
					(pos / 3) + posOffset,
					(pos % 3) + posOffset,
					cDna_.seq_[pos],
					len(protein_.seq_) > pos / 3 ? std::string(1, protein_.seq_[pos / 3]) : "NA");
		}
	}
	infosByAAPos_ = getInfosByAAPos();
}

Bed6RecordCore GeneSeqInfo::genBedFromAAPositions(const uint32_t aaStart,
		const uint32_t aaStop) const {
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(infoTab_.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error table hasn't been set yet"  << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (aaStart >= aaStop) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error aaStart, " << aaStart
				<< " is larger than or equal to aaStop, " << aaStop << "\n";
		throw std::runtime_error { ss.str() };
	}

	if (aaStart >= len(protein_) || aaStop > len(protein_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, out of range, protein length: "
				<< len(protein_) << ", aaStart: " << aaStart << ", aaStop: " << aaStop
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	// uint32_t posOffset = 0;
	// if (pars_.oneBasedPos_) {
	// 	posOffset = 1;
	// }

	auto aaSplit = infoTab_.splitTableOnColumn("aaPos");
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	Bed6RecordCore ret;
	ret.chrom_ = pars_.region_.chrom_;
	ret.strand_ = pars_.region_.reverseSrand_ ? '-' : '+';
	if(pars_.oneBasedPos_){
		ret.name_ = njh::pasteAsStr("[AA", aaStart, "-", "AA", aaStop, "]");
	}else{
		ret.name_ = njh::pasteAsStr("[AA", aaStart, "-", "AA", aaStop, ")");
	}

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(pars_.oneBasedPos_){
		/** @todo the one base does not currently work*/
		if(pars_.region_.reverseSrand_){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			ret.chromEnd_ =   vectorMaximum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStart)).getColumn("gDnaPos"))) + 1 - 1;
			ret.chromStart_ = vectorMinimum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStop)).getColumn("gDnaPos"))) - 1;
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}else{
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			ret.chromStart_ = vectorMinimum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStart)).getColumn("gDnaPos"))) - 1;
			ret.chromEnd_ =   vectorMaximum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStop)).getColumn("gDnaPos"))) + 1 - 1;
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}else{
		if(pars_.region_.reverseSrand_){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			ret.chromEnd_ =   vectorMaximum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStart)).getColumn("gDnaPos"))) + 1;
			//if in the reverse strand get the minimum location of the location prior to the amino acid stop
			ret.chromStart_ = vectorMinimum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStop - 1)).getColumn("gDnaPos")));
			//ret.chromStart_ = vectorMaximum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStop)).getColumn("gDnaPos"))) + 1;
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}else{
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			ret.chromStart_ = vectorMinimum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStart)).getColumn("gDnaPos")));
			//get the max gDNA pos of the amino acid prior to stop and then add one to get the non-inclusive genomic stop position
			ret.chromEnd_ =   vectorMaximum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStop - 1)).getColumn("gDnaPos"))) + 1;
			//ret.chromEnd_ =   vectorMinimum(vecStrToVecNum<uint32_t>(aaSplit.at(estd::to_string(aaStop)).getColumn("gDnaPos")));
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	ret.score_ = ret.length();
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return ret;
}

std::unordered_map<std::string, std::set<uint32_t>> GeneSeqInfo::genGenomicPositionsFromAAPositions(
		const std::vector<uint32_t> &aaPositions) const {
	std::unordered_map<std::string, std::set<uint32_t>> chromPositions;
	for (const auto &pos : aaPositions) {
		auto genomicLocationForAAPos = genBedFromAAPositions(pos, pos + 1);
		for (const auto gPos : iter::range(genomicLocationForAAPos.chromStart_,
				genomicLocationForAAPos.chromEnd_)) {
			chromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
		}
	}
	return chromPositions;
}

Bed6RecordCore GeneSeqInfo::genBedFromCDNAPositions(const uint32_t start,
		const uint32_t stop) const{
	if(infoTab_.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error table hasn't been set yet"  << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (start >= stop) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error start, " << start
				<< " is larger than stop, " << stop << "\n";
		throw std::runtime_error { ss.str() };
	}

	if (start > len(cDna_) || stop > len(cDna_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, out of range, protein length: "
				<< len(protein_) << ", start: " << start << ", stop: " << stop
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	uint32_t posOffset = 0;
	if (pars_.oneBasedPos_) {
		posOffset = 1;
	}
	auto cDnaPosSplit = infoTab_.splitTableOnColumn("cDnaPos");

	Bed6RecordCore ret;
	ret.chrom_ = pars_.region_.chrom_;
	ret.strand_ = pars_.region_.reverseSrand_ ? '-' : '+';
	ret.name_ = njh::pasteAsStr("CDNA", start, "-", "CDNA", stop);
	if(pars_.region_.reverseSrand_){
		ret.chromEnd_ = vectorMaximum(vecStrToVecNum<uint32_t>(cDnaPosSplit.at(estd::to_string(start)).getColumn("gDnaPos"))) + 1 - posOffset;
		ret.chromStart_ = vectorMinimum(vecStrToVecNum<uint32_t>(cDnaPosSplit.at(estd::to_string(stop)).getColumn("gDnaPos"))) - posOffset;
	}else{
		ret.chromStart_ = vectorMinimum(vecStrToVecNum<uint32_t>(cDnaPosSplit.at(estd::to_string(start)).getColumn("gDnaPos"))) - posOffset;
		ret.chromEnd_ = vectorMaximum(vecStrToVecNum<uint32_t>(cDnaPosSplit.at(estd::to_string(stop)).getColumn("gDnaPos"))) + 1 - posOffset;
	}
	ret.score_ = ret.length();
	return ret;
}

GeneSeqInfo::GenePosInfo rowToGenePosInfo(const table & infoTab, const VecStr & row){
	GeneSeqInfo::GenePosInfo info;
	info.gDNAPos_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[infoTab.getColPos("gDnaPos")]);
	info.cDNAPos_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[infoTab.getColPos("cDnaPos")]);
	info.aaPos_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[infoTab.getColPos("aaPos")]);
	info.codonPos_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[infoTab.getColPos("codonPos")]);

	info.base_ = row[infoTab.getColPos("base")].front();
	info.aa_ = row[infoTab.getColPos("aminoAcid")].front();
	return info;
}

std::unordered_map<uint32_t, GeneSeqInfo::GenePosInfo> GeneSeqInfo::getInfosByGDNAPos() const{
	std::unordered_map<uint32_t, GeneSeqInfo::GenePosInfo> ret;
	if(0 == infoTab_.nRow()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error  infoTable not set yet"  << "\n";
		throw std::runtime_error { ss.str() };
	}
	for(const auto & row : infoTab_.content_){
		auto info = rowToGenePosInfo(infoTab_, row);
		ret[info.gDNAPos_] = info;
	}
	return ret;
}



std::unordered_map<uint32_t, GeneSeqInfo::GenePosInfo> GeneSeqInfo::getInfosByCDNAPos() const{
	std::unordered_map<uint32_t, GeneSeqInfo::GenePosInfo> ret;
	if(0 == infoTab_.nRow()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error  infoTable not set yet"  << "\n";
		throw std::runtime_error { ss.str() };
	}
	for(const auto & row : infoTab_.content_){
		auto info = rowToGenePosInfo(infoTab_, row);
		ret[info.gDNAPos_] = info;
	}
	return ret;
}

std::unordered_map<uint32_t, std::tuple<GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo>> GeneSeqInfo::getInfosByAAPos() const{
	std::unordered_map<uint32_t, std::tuple<GenePosInfo,GenePosInfo,GenePosInfo>> ret;
	if(0 == infoTab_.nRow()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error  infoTable not set yet"  << "\n";
		throw std::runtime_error { ss.str() };
	}
	auto infoTabSplit = infoTab_.splitTableOnColumn("aaPos");

	for (const auto & aa : infoTabSplit) {
		if (3 != aa.second.nRow()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error for " << pars_.region_.uid_ << " amino acid split table should contain 3 rows not "
					<< aa.second.nRow() << "\n";
			aa.second.outPutContentOrganized(ss);
//			throw std::runtime_error { ss.str() };
		} else {
			auto codonSplit = aa.second.splitTableOnColumn("codonPos");
			if (!njh::in(std::string("0"), codonSplit)
					|| !njh::in(std::string("1"), codonSplit)
					|| !njh::in(std::string("2"), codonSplit)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", codon column should have 0, 1, and 2 not "
					 << njh::conToStr(getVectorOfMapKeys(codonSplit), ", ") << "\n";
				throw std::runtime_error { ss.str() };
			}
			auto codon0 = rowToGenePosInfo(infoTab_, codonSplit.at("0").content_.front());
			auto codon1 = rowToGenePosInfo(infoTab_, codonSplit.at("1").content_.front());
			auto codon2 = rowToGenePosInfo(infoTab_, codonSplit.at("2").content_.front());
			ret[njh::StrToNumConverter::stoToNum<uint32_t>(aa.first)] = std::make_tuple(codon0, codon1, codon2);
		}
	}
	return ret;
}



}  // namespace njhseq

