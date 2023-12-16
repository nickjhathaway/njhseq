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


void TranslatorByAlignment::TranslatorByAlignmentPars::setOptions(seqSetUp & setUp, bool requireGenome){
	setUp.setOption(knownAminoAcidMutationsFnp_, "--knownAminoAcidChangesFnp",
			"Known Amino Acid Changes, must have at least 2 columns, positions are 1-postion-based (first position is 1), 1)TranscriptID, 2)AAPosition ", false, "Translation Output");
	setUp.setOption(gffFnp_, "--gff,--gffFnp",
			"Gff file to intersect the final haplotypes with genes to get translations", requireGenome || !knownAminoAcidMutationsFnp_.empty(), "Translation Output");
	setUp.setOption(lzPars_.genomeFnp, "--genome,--genomeFnp",
			"Genome file so final haplotypes can be mapped to a genome", requireGenome || !gffFnp_.empty() || !knownAminoAcidMutationsFnp_.empty(), "Translation Output");
	if(njh::endsWith(lzPars_.genomeFnp.string(), ".2bit")){
		lzPars_.genomeFnp.replace_extension("fasta");
	}
	setUp.setOption(useLastz_, "--useLastz", "Use lastz for alignment", false,
			"Translation Output");
	if(!lzPars_.genomeFnp.empty() && bfs::exists(lzPars_.genomeFnp) && !bfs::is_regular_file(lzPars_.genomeFnp)){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr(lzPars_.genomeFnp, " should be a file, not a directory"));
	}

	setUp.setOption(aaExpand_, "--aaExpand",
				"Amount to expand the protein for when aligning ", false, "Translation Output");
	setUp.setOption(useFullProtein_, "--useFullProtein",
				"Just use the full protein for aligning, could be more time consuming ", false, "Translation Output");
	setUp.setOption(allowableStopCodons_, "--allowableStopCodons",
					"The number of stop codons to allow in a translation, if contains more than than the translation will be filtered off", false, "Translation Output");

	bool doNotWriteOutGeneInfos = false;
	setUp.setOption(doNotWriteOutGeneInfos, "--doNotWriteOutGeneInfos",
					"Do Not Write Out Gene Infos that intersect with seqs", false, "Translation Output");

	writeOutGeneInfos_  = !doNotWriteOutGeneInfos;
}


Bed6RecordCore TranslatorByAlignment::TranslateSeqRes::genBedRec() const{
	return Bed6RecordCore(transcriptName_,
						std::get<0>(firstAminoInfo_).aaPos_,
						std::get<0>(lastAminoInfo_).aaPos_ + 1,
						cDna_.name_,
						std::get<0>(lastAminoInfo_).aaPos_ + 1 - std::get<0>(firstAminoInfo_).aaPos_,
						'+');
}

TranslatorByAlignment::Codon TranslatorByAlignment::TranslateSeqRes::getCodonForAARefPos(
		uint32_t aaRefPos) const {
	if (aaRefPos < std::get<0>(firstAminoInfo_).aaPos_
			|| aaRefPos > std::get<0>(lastAminoInfo_).aaPos_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "position: " << aaRefPos
				<< " out of range, current range is:"
				<< std::get<0>(firstAminoInfo_).aaPos_ << "-"
				<< std::get<0>(lastAminoInfo_).aaPos_ << "\n";
		throw std::runtime_error { ss.str() };
	}

	uint32_t alnPosForLoc = getAlnPosForRealPos(refAlnTranslation_.seq_, aaRefPos);
	uint32_t queryCDNAPos = getRealPosForAlnPos(queryAlnTranslation_.seq_, alnPosForLoc) * 3;
	char aa = queryAlnTranslation_.seq_[alnPosForLoc];
//	std::cout << "aa: " << aa << std::endl;
//	std::cout << "alnPosForLoc: " << alnPosForLoc << std::endl;
//	std::cout << "queryAlnTranslation_.seq_.size(): " << queryAlnTranslation_.seq_.size() << std::endl;
//	std::cout << "queryCDNAPos: " << queryCDNAPos << std::endl;
//	std::cout << "cDna_.seq_.size(): " << cDna_.seq_.size() << std::endl;
	if('-' == aa){
		return Codon(aa, std::make_tuple('-', '-', '-'));
	}else{
		if (queryCDNAPos + 2 >= cDna_.seq_.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "" << queryCDNAPos + 2
					<< " is greater than cDna_.seq_.size(): "
					<< cDna_.seq_.size() << "\n";
			throw std::runtime_error { ss.str() };
		}
		return Codon(aa,
				std::make_tuple(
				cDna_.seq_[queryCDNAPos + 0],
				cDna_.seq_[queryCDNAPos + 1],
				cDna_.seq_[queryCDNAPos + 2])
				);
	}
}


TranslatorByAlignment::RunPars::RunPars(){
	realnPars.extendAmount = 100;
}

TranslatorByAlignment::VariantsInfo::VariantsInfo(const Bed3RecordCore & region, seqInfo  refSeq) : region_(region),
		seqBase_(std::move(refSeq)){

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

void VCFOutput::sortRecords() {
	njh::sort(records, [](const VCFRecord & r1, const VCFRecord & r2) {
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

void VCFOutput::writeOutFixedOnly(std::ostream&vcfOut) const {
	vcfOut << "##fileformat=" << vcfFormatVersion_ << std::endl;
	//write out infos
	for (const auto&info: infoEntries_) {
		vcfOut <<"##INFO=<"
		<< "ID=" << info.id_ << ","
		<< "Number=" << info.number_ << ","
		<< "Type=" << info.type_ << ","
		<< "Description=\"" << info.description_ << "\""
		<< ">"
		<< std::endl;
	}
	vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
	for(const auto & rec : records) {
		vcfOut << rec.chrom_
		<< "\t" << rec.pos_
		<< "\t" << rec.id_
		<< "\t" << rec.ref_
		<< "\t" << njh::conToStr(rec.alts_, ",")
		<< "\t" << rec.qual_
		<< "\t" << rec.filter_;
		std::string infoOut;
		for (const auto&info: infoEntries_) {
			if(!infoOut.empty()) {
				infoOut +=";";
			}
			infoOut += info.id_ + "=" + rec.info_.getMeta(info.id_);
		}
		vcfOut << "\t" << infoOut;
		vcfOut << std::endl;
	}
}

void VCFOutput::writeOutFixedAndSampleMeta(std::ostream&vcfOut) const {
	vcfOut << "##fileformat=" << vcfFormatVersion_ << std::endl;
	//write out infos
	for (const auto&info: infoEntries_) {
		vcfOut <<"##INFO=<"
		<< "ID=" << info.id_ << ","
		<< "Number=" << info.number_ << ","
		<< "Type=" << info.type_ << ","
		<< "Description=\"" << info.description_ << "\""
		<< ">"
		<< std::endl;
	}
	//write out formats
	for (const auto&info: formatEntries_) {
		vcfOut <<"##FORMAT=<"
		<< "ID=" << info.id_ << ","
		<< "Number=" << info.number_ << ","
		<< "Type=" << info.type_ << ","
		<< "Description=\"" << info.description_ << "\""
		<< ">"
		<< std::endl;
	}
	//check samples

	vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	std::set<std::string> allSamples;
	if(!records.empty()) {
		auto firstSetOfSamples = njh::vecToSet(getVectorOfMapKeys(records.front().sampleFormatInfos_));
		for(const auto & rec : records) {
			auto currentSetOfSamples = njh::vecToSet(getVectorOfMapKeys(rec.sampleFormatInfos_));
			if(firstSetOfSamples != currentSetOfSamples) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "samples different for record: " << rec.chrom_ << " " << rec.pos_ << " " << rec.ref_ << "\n";
				ss << "first set of samples: " << njh::conToStr(firstSetOfSamples, ",") << "\n";
				ss << "current set of samples: " << njh::conToStr(currentSetOfSamples, ",") << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		vcfOut << "\t" << njh::conToStr(firstSetOfSamples, "\t");
		allSamples = firstSetOfSamples;
	}
	vcfOut << std::endl;
	std::string formatOut;
	for(const auto & format : formatEntries_) {
		if(!formatOut.empty()) {
			formatOut +=":";
		}
		formatOut += format.id_;
	}
	if(!records.empty()) {
		for(const auto & rec : records) {
			vcfOut << rec.chrom_
			<< "\t" << rec.pos_
			<< "\t" << rec.id_
			<< "\t" << rec.ref_
			<< "\t" << njh::conToStr(rec.alts_, ",")
			<< "\t" << rec.qual_
			<< "\t" << rec.filter_;
			std::string infoOut;
			for (const auto&info: infoEntries_) {
				if(!infoOut.empty()) {
					infoOut +=";";
				}
				infoOut += info.id_ + "=" + rec.info_.getMeta(info.id_);
			}
			vcfOut << "\t" << infoOut;
			vcfOut << "\t" << formatOut;
			for(const auto & sampleName : allSamples) {
				const auto & sample = rec.sampleFormatInfos_.at(sampleName);
				std::string formatOutForSample;
				for(const auto & format : formatEntries_) {
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




VCFOutput TranslatorByAlignment::VariantsInfo::createVCFOutputFixed() const {
	VCFOutput ret;
	ret.infoEntries_.emplace_back(
		"DP", "1", "Integer", "Total Allele Depth, sum of AC with rest of depth being ref"
	);
	ret.infoEntries_.emplace_back(
		"NS", "1", "Integer", "Number of Samples With Data"
	);
	ret.infoEntries_.emplace_back(
		"AC", "A", "Integer", "Allele Count, number of microhaplotypes with variant"
	);
	ret.infoEntries_.emplace_back(
		"AF", "A", "Float", "Allele Frequency, calulated AC/DP"
	);
	//vcfOut << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read Depth for for the ref and alt alleles in the order listed, a count of 0 means not detected\">" << std::endl;
	//"FORMAT"
	//std::string format = "AD";
	std::unordered_set<uint32_t> positionsSet;
	for(const auto & snps : snpsFinal){
		positionsSet.emplace(snps.first);
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle deltions at position 0"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		deletionsFinalForVCF[del.first - 1] = del.second;
	}
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto & ins : insertionsFinalForVCF){
		positionsSet.emplace(ins.first);
	}
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto & del : deletionsFinalForVCF){
		positionsSet.emplace(del.first);
	}
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;

	std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
	njh::sort(positions);
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto & pos : positions){
		//add in SNPs and insertions since they have the same chrom,pos,ref
		if (njh::in(pos, insertionsFinalForVCF) || njh::in(pos, snpsFinal)) {
			std::vector<std::string> alts;
			std::vector<uint32_t> altsCounts;
			std::vector<double> altsFreqs;
			if(njh::in(pos, snpsFinal)){
				for(const auto & b : snpsFinal.at(pos)){
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
			VCFOutput::VCFRecord currentRecord;
			currentRecord.chrom_ = region_.chrom_;
			currentRecord.pos_ = pos + 1;
			currentRecord.id_ = ".";
			currentRecord.ref_ = getBaseForGenomicRegion(pos);
			currentRecord.alts_ = alts;
			currentRecord.qual_ = 40;
			currentRecord.filter_ = "PASS";
			currentRecord.info_.addMeta("DP", depthPerPosition.at(pos) );
			currentRecord.info_.addMeta("NS", samplesPerPosition.at(pos).size() );
			currentRecord.info_.addMeta("AC", njh::conToStr(altsCounts, ",") );
			currentRecord.info_.addMeta("AF", njh::conToStr(altsFreqs, ",") );
			ret.records.emplace_back(std::move(currentRecord));
		}
		//add in deletions
		if (njh::in(pos, deletionsFinalForVCF)) {
			for (const auto & d : deletionsFinalForVCF[pos]) {
				VCFOutput::VCFRecord currentRecord;
				currentRecord.chrom_ = region_.chrom_;
				currentRecord.pos_ = pos + 1;
				currentRecord.id_ = ".";
				currentRecord.ref_ = std::string(1, getBaseForGenomicRegion(pos)) + d.first;
				currentRecord.alts_ = VecStr{std::string(1, getBaseForGenomicRegion(pos))};
				currentRecord.qual_ = 40;
				currentRecord.filter_ = "PASS";
				currentRecord.info_.addMeta("DP", depthPerPosition.at(pos) );
				currentRecord.info_.addMeta("NS", samplesPerPosition.at(pos).size() );
				currentRecord.info_.addMeta("AC", d.second );
				currentRecord.info_.addMeta("AF", d.second/static_cast<double>(depthPerPosition.at(pos)) );
				ret.records.emplace_back(std::move(currentRecord));
			}
		}
	}
	return ret;
}

VCFOutput TranslatorByAlignment::VariantsInfo::writeVCF(std::ostream & vcfOut) const{
	auto vcfOutputInfo = createVCFOutputFixed();
	vcfOutputInfo.writeOutFixedOnly(vcfOut);
	return vcfOutputInfo;
// 	std::unordered_set<uint32_t> positionsSet;
// 	for(const auto & snps : snpsFinal){
// 		positionsSet.emplace(snps.first);
// 	}
// 	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 	std::map<uint32_t, std::map<std::string,uint32_t>> insertionsFinalForVCF;
// 	std::map<uint32_t, std::map<std::string,uint32_t>> deletionsFinalForVCF;
// 	for(const auto & ins : insertionsFinal){
// 		if(0 == ins.first ){
// 			std::stringstream ss;
// 			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle insertion at position 0"<< "\n";
// 			throw std::runtime_error{ss.str()};
// 		}
// 		insertionsFinalForVCF[ins.first - 1] = ins.second;
// 	}
// 	for(const auto & del : deletionsFinal){
// 		if(0 == del.first ){
// 			std::stringstream ss;
// 			ss << __PRETTY_FUNCTION__ << ", error " << "can't handle deltions at position 0"<< "\n";
// 			throw std::runtime_error{ss.str()};
// 		}
// 		deletionsFinalForVCF[del.first - 1] = del.second;
// 	}
// 	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 	for(const auto & ins : insertionsFinalForVCF){
// 		positionsSet.emplace(ins.first);
// 	}
// 	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 	for(const auto & del : deletionsFinalForVCF){
// 		positionsSet.emplace(del.first);
// 	}
// 	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
//
// 	vcfOut << "##fileformat=VCFv4.0" << std::endl;
// 	vcfOut << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth, sum of AC with rest of depth being ref\">" << std::endl;
// 	vcfOut << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
// 	vcfOut << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << std::endl;
// 	vcfOut << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">" << std::endl;
//
// 	//vcfOut << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read Depth for for the ref and alt alleles in the order listed, a count of 0 means not detected\">" << std::endl;
//
// 	vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
// 	//"FORMAT"
// 	//std::string format = "AD";
//
//
// 	std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
// 	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 	njh::sort(positions);
// 	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 	for(const auto & pos : positions){
// 		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 		if (njh::in(pos, insertionsFinalForVCF) || njh::in(pos, snpsFinal)) {
// 			vcfOut <<  region_.chrom_
// 					<< "\t" << pos + 1
// 					<< "\t" << "."
// 					<< "\t";
// 			std::vector<std::string> alts;
// 			std::vector<uint32_t> altsCounts;
// 			std::vector<double> altsFreqs;
// 			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 			vcfOut << getBaseForGenomicRegion(pos) << "\t";
// 			if(njh::in(pos, snpsFinal)){
// 				// uint32_t snpCount = 0;
// 				for(const auto & b : snpsFinal.at(pos)){
// 					// snpCount+= b.second;
// 					alts.emplace_back(std::string(1, b.first));
// 					altsCounts.emplace_back(b.second);
// 					altsFreqs.emplace_back(b.second/static_cast<double>(depthPerPosition.at(pos)));
// 				}
// 			}
// 			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 			if (njh::in(pos, insertionsFinalForVCF)) {
// 				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 				for (const auto & ins : insertionsFinalForVCF[pos]) {
// 					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// //					std::cout << "pos: " << pos << std::endl;
// //					std::cout << vectorMinimum(getVectorOfMapKeys(depthPerPosition)) << std::endl;
// //					std::cout << vectorMaximum(getVectorOfMapKeys(depthPerPosition)) << std::endl;
// //					std::cout << "ins.second: " << ins.second << std::endl;
// //					std::cout << "ins.first: " << ins.first << std::endl;
// 					alts.emplace_back(njh::pasteAsStr(getBaseForGenomicRegion(pos), ins.first));
// 					altsCounts.emplace_back(ins.second);
// 					altsFreqs.emplace_back(ins.second/static_cast<double>(depthPerPosition.at(pos)));
// 				}
// 			}
// 			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 			vcfOut << njh::conToStr(alts, ",")
// 			<< "\t40\tPASS\t";
// 			vcfOut
// 					<< "DP=" << depthPerPosition.at(pos) << ";"
// 					<< "NS=" << samplesPerPosition.at(pos).size() << ";"
// 					<< "AC=" << njh::conToStr(altsCounts, ",") << ";"
// 					<< "AF=" << njh::conToStr(altsFreqs, ",")
// 			<< std::endl;
// 		}
// 		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 		if (njh::in(pos, deletionsFinalForVCF)) {
// 			for (const auto & d : deletionsFinalForVCF[pos]) {
// 				vcfOut <<  region_.chrom_
// 						<< "\t" << pos + 1
// 						<< "\t" << "."
// 						<< "\t";
// 				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 				vcfOut << getBaseForGenomicRegion(pos) << d.first
// 				<< "\t" << getBaseForGenomicRegion(pos) << "\t";
// 				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 				vcfOut << "40\tPASS\t";
// 				vcfOut
// 						<< "DP=" << depthPerPosition.at(pos) << ";"
// 						<< "NS=" << samplesPerPosition.at(pos).size() << ";"
// 						<< "AC=" << d.second << ";"
// 						<< "AF=" << d.second/static_cast<double>(depthPerPosition.at(pos))
// 				<< std::endl;
// 				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 			}
// 		}
// 		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
// 	}
// 	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
}



VCFOutput TranslatorByAlignment::VariantsInfo::writeVCF(const OutOptions & vcfOutOpts) const {
	OutputStream vcfOut(vcfOutOpts);
	return writeVCF(vcfOut);
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

uint32_t TranslatorByAlignment::VariantsInfo::getFinalNumberOfSegratingSites() const{
	std::unordered_set<uint32_t> idPositions;
	for(const auto & snp : snpsFinal){
		idPositions.emplace(snp.first);
	}
	for(const auto & ins : insertionsFinal){
		idPositions.emplace(ins.first);
	}
	for(const auto & del : deletionsFinal){
		for(const auto & seq : del.second){
			idPositions.emplace(del.first + seq.first.size());
		}
	}
	return idPositions.size();
}

void TranslatorByAlignment::VariantsInfo::writeOutSNPsInfo(
		const OutOptions & outOpts,
		const std::string & name,
		const std::set<uint32_t> & snpPositions,
		bool oneBased){
	OutputStream out(outOpts);
	if(oneBased){
		out << njh::conToStr(TranslatorByAlignment::VariantsInfo::SNPHeaderAminoAcid(), "\t") << std::endl;
	}else{
		out << njh::conToStr(TranslatorByAlignment::VariantsInfo::SNPHeaderGenomic(), "\t") << std::endl;
	}
	writeOutSNPsInfo(out, name, snpPositions, oneBased);
}

void TranslatorByAlignment::VariantsInfo::writeOutSNPsInfo(std::ostream & out,
		const std::string & name,
		const std::set<uint32_t> & snpPositions,
		bool oneBased){

	for(const auto & snpPos : snpPositions){
		if(!njh::in(snpPos, allBases)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no info for snp position: " << snpPos << "\n";
			throw std::runtime_error{ss.str()};
		}
		for(const auto & aa : allBases[snpPos]){
			out << name
					<< "\t" << (oneBased ? snpPos +1 : snpPos)
					<< "\t" << getBaseForGenomicRegion(snpPos)
					<< "\t" << aa.first
					<< "\t" << aa.second
					<< "\t" << aa.second/static_cast<double>(depthPerPosition[snpPos])
					<< "\t" << depthPerPosition[snpPos]
					<< "\t" << samplesPerPosition[snpPos].size() << std::endl;
		}
	}
}

void TranslatorByAlignment::VariantsInfo::writeOutSNPsFinalInfo(std::ostream & out,
		const std::string & name,
		bool oneBased){
	writeOutSNPsInfo(out, name, njh::getSetOfMapKeys(snpsFinal), oneBased);
}

void TranslatorByAlignment::VariantsInfo::writeOutSNPsFinalInfo(
		const OutOptions & outOpts,
		const std::string & name,
		bool oneBased){
	OutputStream out(outOpts);
	if(oneBased){
		out << njh::conToStr(TranslatorByAlignment::VariantsInfo::SNPHeaderAminoAcid(), "\t") << std::endl;
	}else{
		out << njh::conToStr(TranslatorByAlignment::VariantsInfo::SNPHeaderGenomic(), "\t") << std::endl;
	}
	writeOutSNPsFinalInfo(out, name, oneBased);
}

void TranslatorByAlignment::VariantsInfo::writeOutSNPsAllInfo(std::ostream & out,
		const std::string & name,
		bool oneBased){
	writeOutSNPsInfo(out, name, njh::getSetOfMapKeys(allBases), oneBased);
}
void TranslatorByAlignment::VariantsInfo::writeOutSNPsAllInfo(
		const OutOptions & outOpts,
		const std::string & name,
		bool oneBased){
	OutputStream out(outOpts);
	if(oneBased){
		out << njh::conToStr(TranslatorByAlignment::VariantsInfo::SNPHeaderAminoAcid(), "\t") << std::endl;
	}else{
		out << njh::conToStr(TranslatorByAlignment::VariantsInfo::SNPHeaderGenomic(), "\t") << std::endl;
	}
	writeOutSNPsInfo(out, name, njh::getSetOfMapKeys(allBases), oneBased);

}


Bed3RecordCore TranslatorByAlignment::VariantsInfo::getVariableRegion() {
	if (!variablePositons_.empty()) {
		return Bed3RecordCore(region_.chrom_,
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
	uint32_t queryAlnStart = alignedQuerySeq.find_first_not_of('-');
	uint32_t queryAlnEnd = alignedQuerySeq.find_last_not_of('-');
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
			aligner & alignerObj) const{
	std::unordered_map<std::string, TranslateSeqRes> ret;
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto & transcript : currentGene.mRNAs_){
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto currentTranscriptInfo = njh::mapAt(transcriptInfosForGene, transcript->getIDAttr());
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto genePosInfoByGDna = currentTranscriptInfo->getInfosByGDNAPos();
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		// bool endsAtStopCodon = false;
		uint32_t transStart = 0;
		seqInfo balnSeq(realigned.querySeq_.name_);
		std::vector<uint32_t> codons;
		std::vector<GFFCore> cDNAIntersectedWith;
		for (const auto & cDna : njh::mapAt(currentGene.CDS_, transcript->getIDAttr())) {
			if (realigned.gRegion_.overlaps(*cDna, 3)) {
				cDNAIntersectedWith.emplace_back(*cDna);
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "cDNAIntersectedWith.size(): " << cDNAIntersectedWith.size() << std::endl;

		uint32_t cdnaGenomicStart = std::numeric_limits<uint32_t>::max();
		uint32_t cdnaGenomicEndInconclusive = std::numeric_limits<uint32_t>::max();

		if(cDNAIntersectedWith.empty()){
			continue;
		} else {
			if (cDNAIntersectedWith.size() == 1
					&& realigned.gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
					&& realigned.gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				balnSeq = realigned.querySeq_;
				if (currentGene.gene_->isReverseStrand()) {
					// if (njh::mapAt(genePosInfoByGDna,realigned.gRegion_.start_).cDNAPos_
					// 		== currentTranscriptInfo->cDna_.seq_.size() - 1) {
					// 	endsAtStopCodon = true;
					// }
					uint32_t gPos = realigned.gRegion_.end_ - 1;
					auto codon = njh::mapAt(genePosInfoByGDna,gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = njh::mapAt(genePosInfoByGDna,gPos).codonPos_;
						++transStart;
					}
					cdnaGenomicEndInconclusive = gPos;
					cdnaGenomicStart = realigned.gRegion_.start_;
				} else {
					// if (njh::mapAt(genePosInfoByGDna,realigned.gRegion_.end_ - 1).cDNAPos_
					// 		== currentTranscriptInfo->cDna_.seq_.size() - 1) {
					// 	endsAtStopCodon = true;
					// }
					uint32_t gPos = realigned.gRegion_.start_;
					uint32_t codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
						++transStart;
					}
					cdnaGenomicEndInconclusive = realigned.gRegion_.end_ - 1;
					cdnaGenomicStart = gPos;
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
					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
					auto cDnaStop = cDNAIntersectedWith.back().end_;
					uint32_t gPos = std::min(cDnaStop, realigned.gRegion_.end_) - 1;
					auto codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
						++transStart;
					}
				} else {
					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
					auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
					uint32_t gPos = std::max(cDnaStart, realigned.gRegion_.start_);
					uint32_t codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
						++transStart;
					}
					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::vector<uint32_t> starts;
				std::vector<uint32_t> ends;
				for (const auto & cDna : cDNAIntersectedWith) {
					auto cDnaStart = cDna.start_ - 1;
					auto detStart = std::max(cDnaStart, realigned.gRegion_.start_);
					auto detStop = std::min(cDna.end_, realigned.gRegion_.end_);
					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
					ends.emplace_back(detStop);
					starts.emplace_back(detStart);
					detStart -= realigned.gRegion_.start_;
					detStop -= realigned.gRegion_.start_;
//					realigned.alnRefSeq_.outPutSeqAnsi(std::cout);
//					realigned.alnQuerySeq_.outPutSeqAnsi(std::cout);
//					std::cout << "realigned.alnRefSeq_.seq_.size() : " << realigned.alnRefSeq_.seq_.size() << std::endl;
//					std::cout << "detStart    : " << detStart << std::endl;
//					std::cout << "detStop - 1 : " << detStop - 1 << std::endl;

					auto alnStart = getAlnPosForRealPos(realigned.alnRefSeq_.seq_,detStart);
					auto alnStop = getAlnPosForRealPos(realigned.alnRefSeq_.seq_, detStop - 1);
					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
					balnSeq.append(
							realigned.alnQuerySeq_.getSubRead(alnStart,
									alnStop - alnStart + 1));
				}
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				uint32_t cDnaStart = *std::min_element(starts.begin(), starts.end());
				uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
				cdnaGenomicStart = cDnaStart;
				cdnaGenomicEndInconclusive = cDnaStop -1 ;

				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << "cDnaStart:" << cDnaStart << std::endl;
//				std::cout << "cDnaStop:" << cDnaStop << std::endl;
//				std::cout << "realigned.gRegion_.start_:" << realigned.gRegion_.start_ << std::endl;
//				std::cout << "realigned.gRegion_.end_:" << realigned.gRegion_.end_ << std::endl;
//				std::cout << "cDNAIntersectedWith.front().start_:" << cDNAIntersectedWith.front().start_ << std::endl;
//				std::cout << "cDNAIntersectedWith.front().end_:" << cDNAIntersectedWith.front().end_ << std::endl;
//
//				std::cout << "cdnaGenomicStart          :" << cdnaGenomicStart << std::endl;
//				std::cout << "cdnaGenomicEndInconclusive:" << cdnaGenomicEndInconclusive << std::endl;

				// if (currentGene.gene_->isReverseStrand()) {
				// 	if (njh::mapAt(genePosInfoByGDna, cDnaStart).cDNAPos_
				// 			== currentTranscriptInfo->cDna_.seq_.size() - 1) {
				// 		endsAtStopCodon = true;
				// 	}
				// } else {
				// 	if (njh::mapAt(genePosInfoByGDna, cDnaStop - 1).cDNAPos_
				// 			== currentTranscriptInfo->cDna_.seq_.size() - 1) {
				// 		endsAtStopCodon = true;
				// 	}
				// }
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				balnSeq.removeGaps();
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			}
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (currentGene.gene_->isReverseStrand()) {
				balnSeq.reverseComplementRead(false, true);
			}
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			bool forceMStart = false;
			uint32_t aaStart = std::min(njh::mapAt(genePosInfoByGDna, cdnaGenomicStart).aaPos_, njh::mapAt(genePosInfoByGDna, cdnaGenomicEndInconclusive).aaPos_);
			uint32_t aaEnd =   std::max(njh::mapAt(genePosInfoByGDna, cdnaGenomicStart).aaPos_, njh::mapAt(genePosInfoByGDna, cdnaGenomicEndInconclusive).aaPos_);
			if(0 == aaStart){
				forceMStart = true;
			}
			auto balnSeqTrans = balnSeq.translateRet(false, false, transStart, forceMStart);
			MetaDataInName transMeta;
			transMeta.addMeta("transcript", transcript->getIDAttr());
			balnSeqTrans.name_ += transMeta.createMetaName();

			if(pars_.useFullProtein_){
				alignerObj.alignCacheGlobal(currentTranscriptInfo->protein_, balnSeqTrans);
			}else{

				if(aaStart > pars_.aaExpand_){
					aaStart -= pars_.aaExpand_;
				}else{
					aaStart  = 0;
				}
				if(len(currentTranscriptInfo->protein_) - aaEnd > pars_.aaExpand_){
					aaEnd += pars_.aaExpand_;
				}else{
					aaEnd = len(currentTranscriptInfo->protein_);
				}
				auto subProtein = currentTranscriptInfo->protein_.getSubRead(aaStart, aaEnd - aaStart);
				alignerObj.alignCacheGlobal(subProtein, balnSeqTrans);
				if(aaStart > 0){
					alignerObj.alignObjectA_.seqBase_.prepend(currentTranscriptInfo->protein_.getSubRead(0, aaStart));
					alignerObj.alignObjectB_.seqBase_.prepend(std::string(aaStart, '-'));
				}
				if(aaEnd != len(currentTranscriptInfo->protein_)){
					alignerObj.alignObjectA_.seqBase_.append(currentTranscriptInfo->protein_.getSubRead(aaEnd));
					alignerObj.alignObjectB_.seqBase_.append(std::string(len(currentTranscriptInfo->protein_) - aaEnd, '-'));
				}
			}
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			alignerObj.profilePrimerAlignment(currentTranscriptInfo->protein_,balnSeqTrans);
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			TranslateSeqRes tRes;
			uint32_t cDnaLenRaw = len(balnSeq) - transStart;
			uint32_t cDnaLen = cDnaLenRaw - (cDnaLenRaw %3);
			uint32_t firstAmino = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-'));
			uint32_t lastAmino = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-'));
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "firstAmino: " << firstAmino << std::endl;
//			std::cout << "lastAmino : " << lastAmino << std::endl;
			lastAmino = std::min<uint32_t>(lastAmino,len(currentTranscriptInfo->protein_) - 1);
			//auto aminoInfos = currentTranscriptInfo->getInfosByAAPos();
			tRes.firstAminoInfo_ = njh::mapAt(currentTranscriptInfo->infosByAAPos_, firstAmino);
			tRes.lastAminoInfo_ =  njh::mapAt(currentTranscriptInfo->infosByAAPos_, lastAmino);
			tRes.cDna_ = balnSeq.getSubRead(transStart, cDnaLen);
			tRes.transcriptName_ = transcript->getIDAttr();
			tRes.translation_ = balnSeqTrans;
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
			tRes.refAlnTranslation_ = alignerObj.alignObjectA_.seqBase_;
			tRes.queryAlnTranslation_ = alignerObj.alignObjectB_.seqBase_;
			tRes.comp_ = alignerObj.comp_;
			ret[transcript->getIDAttr()] = tRes;
			 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return ret;
}

TranslatorByAlignment::TranslatorByAlignment(TranslatorByAlignmentPars  pars): pars_(std::move(pars)){
	njh::sys::requireExternalProgramThrow("samtools");
	if(!pars_.useLastz_){
		njh::sys::requireExternalProgramThrow("bowtie2");
	}else{
		njh::sys::requireExternalProgramThrow("lastz");
	}
	VecStr warnings;
	if(!pars_.knownAminoAcidMutationsFnp_.empty()){
		if(!pars_.knownAminoAcidMutationsFnp_.empty()){
			if(pars_.gffFnp_.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "if supplying known amino acid positions than must also supply --gffFnp file"<< "\n";
				warnings.emplace_back(ss.str());
			}
		}

		table knownAminoAcidChanges(pars_.knownAminoAcidMutationsFnp_, "\t", true);
		VecStr originalColNames = knownAminoAcidChanges.columnNames_;
		njh::for_each(knownAminoAcidChanges.columnNames_, [](std::string & col){
			njh::strToLower(col);
		});
		knownAminoAcidChanges.setColNamePositions();
		knownAminoAcidChanges.checkForColumnsThrow(VecStr{"transcriptid", "aaposition"}, __PRETTY_FUNCTION__);
		if(knownAminoAcidChanges.nRow() > 0){
			for(const auto & row : knownAminoAcidChanges){
				if(std::all_of(row.begin(), row.end(), [](const std::string & element){
					return element.empty();
				})){
					continue;
				}
				auto transciprtId = row[knownAminoAcidChanges.getColPos("transcriptid")];
				auto aaPosition = njh::StrToNumConverter::stoToNum<uint32_t>(row[knownAminoAcidChanges.getColPos("aaposition")]);
				if(njh::in(aaPosition, knownAminoAcidPositions_[transciprtId])){
					warnings.emplace_back(njh::pasteAsStr("already have aaposition ", aaPosition, " for transcriptID: ", transciprtId));
				}
				knownAminoAcidPositions_[transciprtId].emplace(aaPosition);
				if(originalColNames.size() > 2){
					for(const auto && colPos : iter::range(knownAminoAcidChanges.columnNames_.size())){
						if(!njh::in(knownAminoAcidChanges.columnNames_[colPos], VecStr{"transcriptid", "aaposition"})){
							metaDataAssociatedWithAminoacidPosition_[transciprtId][aaPosition].addMeta(originalColNames[colPos], row[colPos]);
						}
					}
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



void TranslatorByAlignment::TranslatorByAlignmentResult::writeOutSeqAlnIndvVars(const OutOptions & outopts) const{
	OutputStream individualAlnedVariantInfoOut(outopts);
	individualAlnedVariantInfoOut << njh::conToStr(concatVecs(DistanceComp::BasicInfoHeader(), VecStr{"variantInPopulation"}), "\t") << std::endl;
	auto seqAlnNames = getVectorOfMapKeys(seqAlns_);
	njh::sort(seqAlnNames);
	for(const auto & seqAlnName : seqAlnNames){
		for(const auto & alns : seqAlns_.at(seqAlnName)){
			// alns.comp_.distances_.writeBasicInfo(individualAlnedVariantInfoOut, alns.gRegion_.genBed3RecordCore(), seqAlnName);
			const auto gPos = alns.gRegion_.genBed3RecordCore();
			uint32_t totalVariants = alns.comp_.distances_.mismatches_.size() + alns.comp_.distances_.alignmentGaps_.size();
			//mismatches
			for(const auto & m : alns.comp_.distances_.mismatches_){
				bool variantInPop = false;
				if(njh::in(gPos.chrom_, seqVariants_)) {
					variantInPop = njh::in(gPos.chromStart_ + m.second.refBasePos, seqVariants_.at(gPos.chrom_).snpsFinal) &&
						njh::in(m.second.seqBase, seqVariants_.at(gPos.chrom_).snpsFinal.at(gPos.chromStart_ + m.second.refBasePos)) &&
						seqVariants_.at(gPos.chrom_).snpsFinal.at(gPos.chromStart_ + m.second.refBasePos).at(m.second.seqBase) > 0;
				}
				individualAlnedVariantInfoOut << njh::conToStr(toVecStr(
						gPos.chrom_,
						gPos.chromStart_ + m.second.refBasePos,
						gPos.chromStart_ + m.second.refBasePos + 1,
						seqAlnName,
						"SNP",
						m.second.refBase,
						m.second.seqBase,
						totalVariants,
						variantInPop), "\t") << "\n";
			}
			//INDELs
			for(const auto & indel: alns.comp_.distances_.alignmentGaps_){
				if(indel.second.ref_){
					//gap is in reference so a insertion
					bool variantInPop = false;
					if(njh::in(gPos.chrom_, seqVariants_)) {
						variantInPop = njh::in(gPos.chromStart_ + indel.second.refPos_, seqVariants_.at(gPos.chrom_).insertionsFinal) &&
							njh::in(indel.second.gapedSequence_, seqVariants_.at(gPos.chrom_).insertionsFinal.at(gPos.chromStart_ + indel.second.refPos_)) &&
							seqVariants_.at(gPos.chrom_).insertionsFinal.at(gPos.chromStart_ + indel.second.refPos_).at(indel.second.gapedSequence_) > 0;
					}
					individualAlnedVariantInfoOut << njh::conToStr(toVecStr(
							gPos.chrom_,
							gPos.chromStart_ + indel.second.refPos_,
							gPos.chromStart_ + indel.second.refPos_ + 1,
							seqAlnName,
							"insertion",
							std::string(indel.second.gapedSequence_.size(), '-'),
							indel.second.gapedSequence_,
							totalVariants,
							variantInPop), "\t") << "\n";
				} else {
					//gap is in query so a deletion
					bool variantInPop = false;
					if(njh::in(gPos.chrom_, seqVariants_)) {
						variantInPop = njh::in(gPos.chromStart_ + indel.second.refPos_, seqVariants_.at(gPos.chrom_).deletionsFinal) &&
							njh::in(indel.second.gapedSequence_, seqVariants_.at(gPos.chrom_).deletionsFinal.at(gPos.chromStart_ + indel.second.refPos_)) &&
							seqVariants_.at(gPos.chrom_).deletionsFinal.at(gPos.chromStart_ + indel.second.refPos_).at(indel.second.gapedSequence_) > 0;
					}
					individualAlnedVariantInfoOut << njh::conToStr(toVecStr(
							gPos.chrom_,
							gPos.chromStart_ + indel.second.refPos_,
							gPos.chromStart_ + indel.second.refPos_ + indel.second.gapedSequence_.size(),
							seqAlnName,
							"deletion",
							indel.second.gapedSequence_,
							std::string(indel.second.gapedSequence_.size(), '-'),
							totalVariants,
							variantInPop), "\t") << "\n";
				}
			}
		}
	}
}





std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::genSeqSNPTypedStr() const{
	std::map<std::string, std::string> ret;



	for(const auto & seqAln : seqAlns_){
		std::string typed;
		for(const auto & aln : seqAln.second){
			if(!typed.empty()){
				typed += ";";
			}
			auto snpPositions = njh::getVecOfMapKeys(seqVariants_.at(aln.gRegion_.chrom_).snpsFinal);
			VecStr snpTyped;
			std::string currentType = aln.gRegion_.chrom_;
			for (const auto snpPos: snpPositions) {
				auto snpPosRel = snpPos - aln.gRegion_.start_;
				if (snpPos >= aln.gRegion_.start_ && snpPos < aln.gRegion_.end_) {
//					aln.alnRefSeq_.outPutSeqAnsi(std::cout);
//					aln.alnQuerySeq_.outPutSeqAnsi(std::cout);
//					std::cout << "aln.alnRefSeq_.seq_.size(): " << aln.alnRefSeq_.seq_.size() << std::endl;
//					std::cout << "aln.alnQuerySeq_.seq_.size(): " << aln.alnQuerySeq_.seq_.size() << std::endl;
//					std::cout << "aln.alnRefSeq_.seq_.size(): " << aln.refSeq_.seq_.size() << std::endl;
//					std::cout << "aln.alnQuerySeq_.seq_.size(): " << aln.querySeq_.seq_.size() << std::endl;
//					std::cout << "snpPosRel: " << snpPosRel << std::endl;
//					std::cout << "snpPos: " << snpPos << std::endl;
//					std::cout << aln.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
					snpTyped.emplace_back(
									njh::pasteAsStr(snpPos, aln.alnQuerySeq_.seq_[getAlnPosForRealPos(aln.alnRefSeq_.seq_, snpPosRel)]));
//					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				} else {
					snpTyped.emplace_back(njh::pasteAsStr(snpPos, "N"));
				}
			}
			typed += njh::pasteAsStr(currentType, "--", njh::conToStr(snpTyped, ":"));
		}
		ret[seqAln.first] = typed;
	}
	for(const auto & seqName : seqsUnableToBeMapped_){
		ret[seqName] = "Unmappable";
	}
	return ret;
}


std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::translated_genAATypedStr() const{
	std::map<std::string, std::string> ret;
	for(const auto & seqName : translated_fullAATypedWithCodonInfo_){
		std::unordered_map<std::string, VecStr> perTranscript;
		for(const auto & type : seqName.second){
			perTranscript[type.transcriptName_].emplace_back(njh::pasteAsStr(type.zeroBasedPos_ + 1, type.cod_.aa_));
		}
		std::string typed;
		for(const auto & trans : perTranscript){
			if(!typed.empty()){
				typed += ";";
			}
			typed = njh::pasteAsStr(trans.first, "--", njh::conToStr(trans.second, ":"));
		}
		ret[seqName.first] = typed;
	}
	return ret;
}



std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::translated_genAATypedStrOnlyKnowns() const{
	std::map<std::string, std::string> ret;
	for(const auto & seqName : translated_fullAATypedWithCodonInfo_){
		std::unordered_map<std::string, VecStr> perTranscript;
		for(const auto & type : seqName.second){
			if(type.knownMut_){
				perTranscript[type.transcriptName_].emplace_back(njh::pasteAsStr(type.zeroBasedPos_ + 1, type.cod_.aa_));
			}
		}
		std::string typed;
		for(const auto & trans : perTranscript){
			if(!typed.empty()){
				typed += ";";
			}
			typed = njh::pasteAsStr(trans.first, "--", njh::conToStr(trans.second, ":"));
		}
		ret[seqName.first] = typed;
	}
	return ret;
}

std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::translated_genAATypedStrOnlyPopVariant() const{
	std::map<std::string, std::string> ret;
	for(const auto & seqName : translated_variantAATypedWithCodonInfo_){
		std::unordered_map<std::string, VecStr> perTranscript;
		for(const auto & type : seqName.second){
			perTranscript[type.transcriptName_].emplace_back(njh::pasteAsStr(type.zeroBasedPos_ + 1, type.cod_.aa_));
		}
		std::string typed;
		for(const auto & trans : perTranscript){
			if(!typed.empty()){
				typed += ";";
			}
			typed = njh::pasteAsStr(trans.first, "--", njh::conToStr(trans.second, ":"));
		}
		ret[seqName.first] = typed;
	}
	return ret;
}



std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::genAATypedStr() const{
	std::map<std::string, std::string> ret;
	for(const auto & seqName : fullAATypedWithCodonInfo_){
		std::unordered_map<std::string, VecStr> perTranscript;
		for(const auto & type : seqName.second){
			perTranscript[type.transcriptName_].emplace_back(njh::pasteAsStr(type.zeroBasedPos_ + 1, type.cod_.aa_));
		}
		std::string typed;
		for(const auto & trans : perTranscript){
			if(!typed.empty()){
				typed += ";";
			}
			typed = njh::pasteAsStr(trans.first, "--", njh::conToStr(trans.second, ":"));
		}
		ret[seqName.first] = typed;
	}

	for(const auto & seqName : seqsUnableToBeMapped_){
		ret[seqName] = "Unmappable";
	}
	for(const auto & seqName : seqsTranslationFiltered_){
		ret[seqName] = "Untranslatable";
	}
	return ret;
}

std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::genAATypedStrOnlyKnowns() const{
	std::map<std::string, std::string> ret;
	for(const auto & seqName : fullAATypedWithCodonInfo_){
		std::unordered_map<std::string, VecStr> perTranscript;
		for(const auto & type : seqName.second){
			if(type.knownMut_){
				perTranscript[type.transcriptName_].emplace_back(njh::pasteAsStr(type.zeroBasedPos_ + 1, type.cod_.aa_));
			}
		}
		std::string typed;
		for(const auto & trans : perTranscript){
			if(!typed.empty()){
				typed += ";";
			}
			typed = njh::pasteAsStr(trans.first, "--", njh::conToStr(trans.second, ":"));
		}
		ret[seqName.first] = typed;
	}
	for(const auto & seqName : seqsUnableToBeMapped_){
		ret[seqName] = "Unmappable";
	}
	for(const auto & seqName : seqsTranslationFiltered_){
		ret[seqName] = "Untranslatable";
	}
	return ret;
}

std::map<std::string, std::string> TranslatorByAlignment::TranslatorByAlignmentResult::genAATypedStrOnlyPopVariant() const{
	std::map<std::string, std::string> ret;
	for(const auto & seqName : variantAATypedWithCodonInfo_){
		std::unordered_map<std::string, VecStr> perTranscript;
		for(const auto & type : seqName.second){
			perTranscript[type.transcriptName_].emplace_back(njh::pasteAsStr(type.zeroBasedPos_ + 1, type.cod_.aa_));
		}
		std::string typed;
		for(const auto & trans : perTranscript){
			if(!typed.empty()){
				typed += ";";
			}
			typed = njh::pasteAsStr(trans.first, "--", njh::conToStr(trans.second, ":"));
		}
		ret[seqName.first] = typed;
	}
	for(const auto & seqName : seqsUnableToBeMapped_){
		ret[seqName] = "Unmappable";
	}
	for(const auto & seqName : seqsTranslationFiltered_){
		ret[seqName] = "Untranslatable";
	}
	return ret;
}


void TranslatorByAlignment::TranslatorByAlignmentResult::writeOutAATypedInfo(const OutOptions &outOpts) const {
	OutputStream typedOut(outOpts);
	typedOut << "name\tfullTyped\tknownTyped\tvariantTyped" << std::endl;
	auto fullTyped = genAATypedStr();
	auto knownTyped = genAATypedStrOnlyKnowns();
	auto variantTyped = genAATypedStrOnlyPopVariant();

	for (const auto &seqName : getAllSeqNames()) {
		typedOut << seqName << "\t" << fullTyped[seqName] << "\t"
				<< knownTyped[seqName] << "\t" << variantTyped[seqName];
		typedOut << std::endl;
	}
}


void TranslatorByAlignment::TranslatorByAlignmentResult::writeOutTranslatedIndvVars(const OutOptions & outOpts, const std::unordered_map<std::string, std::set<uint32_t>> & knownAminoAcidPositions1Based) const {
	OutputStream individualVariantInfoOut(outOpts);
	individualVariantInfoOut << njh::conToStr(VecStr{"chrom", "start", "end", "queryName", "type", "refSeq", "querySeq", "refCodonSeq", "queryCodonSeq", "knownAAChange", "variantTotal", "genomicID", "variantInPopulation"}, "\t") << std::endl;
	auto seqNames = njh::getVecOfMapKeys(translations_);
	njh::sort(seqNames);
	std::unordered_map<std::string, std::set<uint32_t>> knownMutationsLocationsZeroBased;
	for(const auto & transcript : knownAminoAcidPositions1Based){
		for(const auto pos : transcript.second){
			knownMutationsLocationsZeroBased[transcript.first].emplace(pos - 1);
		}
	}

//	//this includes both known mutations as well as locations with significant variation
//	std::set<uint32_t> allLocations = getAllInterestingAAPosZeroBased(varPerTrans.first, ret);
//	for (auto & seqName : ret.translations_) {
//		if (njh::in(varPerTrans.first, seqName.second)) {
//			for (const auto & loc : allLocations) {
//				if(loc < std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ || loc > std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_){
//					//location is not within the aligned translation
//					continue;
//				}
//				auto codon = seqName.second[varPerTrans.first].getCodonForAARefPos(loc);
//				ret.fullAATypedWithCodonInfo_[seqName.first].emplace_back(
//						TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
//								njh::in(loc, knownMutationsLocationsZeroBased)));
//			}
//		}
//	}

	for(const auto & seqName : seqNames){
		for(const auto & transcript : translations_.at(seqName)){
			auto codonInfoByAAPos = translationInfoForTranscirpt_.at(transcript.first)->getInfosByAAPos();
			auto gPos = njh::mapAt(proteinVariants_, transcript.first).region_;
			uint32_t totalVariants = transcript.second.comp_.distances_.mismatches_.size() + transcript.second.comp_.distances_.alignmentGaps_.size();
			//mismatches
			for(const auto & m : transcript.second.comp_.distances_.mismatches_){

				//getting reference codon info
				std::string referenceCodonSeq = njh::pasteAsStr(
					std::get<0>(codonInfoByAAPos[m.second.refBasePos]).base_,
					std::get<1>(codonInfoByAAPos[m.second.refBasePos]).base_,
					std::get<2>(codonInfoByAAPos[m.second.refBasePos]).base_);

				//getting reference genomic info
				auto genomicLocation = translationInfoForTranscirpt_.at(transcript.first)->genBedFromAAPositions(m.second.refBasePos, m.second.refBasePos + 1);

				//getting query codon seq info
				auto codon = translations_.at(seqName).at(transcript.first).getCodonForAARefPos(m.second.refBasePos);
				std::string codonSeq = njh::pasteAsStr(std::get<0>(codon.bases_), std::get<1>(codon.bases_), std::get<2>(codon.bases_));
				bool variantInPop = false;
				if(njh::in(gPos.chrom_, proteinVariants_)) {
					variantInPop = njh::in(gPos.chromStart_ + m.second.refBasePos, proteinVariants_.at(gPos.chrom_).snpsFinal) &&
						njh::in(m.second.seqBase, proteinVariants_.at(gPos.chrom_).snpsFinal.at(gPos.chromStart_ + m.second.refBasePos)) &&
						proteinVariants_.at(gPos.chrom_).snpsFinal.at(gPos.chromStart_ + m.second.refBasePos).at(m.second.seqBase) > 0;
				}
				individualVariantInfoOut << njh::conToStr(toVecStr(
						gPos.chrom_,
						gPos.chromStart_ + m.second.refBasePos + 1,
						gPos.chromStart_ + m.second.refBasePos + 1,
						seqName,
						"SNP",
						m.second.refBase,
						m.second.seqBase,
						referenceCodonSeq,
						codonSeq,
						njh::boolToStr(njh::in(m.second.refBasePos, knownMutationsLocationsZeroBased[transcript.first])),
						totalVariants,
						genomicLocation.genUIDFromCoords(),
						variantInPop), "\t") << "\n";
			}
			//INDELs
			for(const auto & indel: transcript.second.comp_.distances_.alignmentGaps_){
				if(indel.second.ref_){
					//gap is in reference so a insertion
					//getting reference codon info
					std::string referenceCodonSeq = std::string(indel.second.gapedSequence_.size()*3, '-');

					//getting reference genomic info
					auto genomicLocation = translationInfoForTranscirpt_.at(transcript.first)->genBedFromAAPositions(indel.second.refPos_, indel.second.refPos_ + 1);
					bool variantInPop = false;
					if(njh::in(gPos.chrom_, proteinVariants_)) {
						variantInPop = njh::in(gPos.chromStart_ + indel.second.refPos_, proteinVariants_.at(gPos.chrom_).insertionsFinal) &&
							njh::in(indel.second.gapedSequence_, proteinVariants_.at(gPos.chrom_).insertionsFinal.at(gPos.chromStart_ + indel.second.refPos_)) &&
							proteinVariants_.at(gPos.chrom_).insertionsFinal.at(gPos.chromStart_ + indel.second.refPos_).at(indel.second.gapedSequence_) > 0;
					}
					individualVariantInfoOut << njh::conToStr(toVecStr(
							gPos.chrom_,
							gPos.chromStart_ + indel.second.refPos_ + 1,
							gPos.chromStart_ + indel.second.refPos_ + 1,
							seqName,
							"insertion",
							std::string(indel.second.gapedSequence_.size(), '-'),
							indel.second.gapedSequence_,
							referenceCodonSeq,
							translations_.at(seqName).at(transcript.first).cDna_.seq_.substr(indel.second.seqPos_ *3,indel.second.gapedSequence_.size() * 3),
							njh::boolToStr(false),
							totalVariants,
						  genomicLocation.genUIDFromCoords(),
						  variantInPop), "\t") << "\n";
				} else {
					//gap is in query so a deletion
					//getting reference codon info
					std::string referenceCodonSeq;
					for (const auto aaPos: iter::range<uint32_t>(indel.second.refPos_,
					                                             indel.second.refPos_ + indel.second.gapedSequence_.size())) {
						referenceCodonSeq += njh::pasteAsStr(
							std::get<0>(codonInfoByAAPos[aaPos]).base_,
							std::get<1>(codonInfoByAAPos[aaPos]).base_,
							std::get<2>(codonInfoByAAPos[aaPos]).base_);
					}


					//getting reference genomic info
					auto genomicLocation = translationInfoForTranscirpt_.at(transcript.first)->genBedFromAAPositions(indel.second.refPos_, indel.second.refPos_ + indel.second.gapedSequence_.size());
					bool variantInPop = false;
					if(njh::in(gPos.chrom_, proteinVariants_)) {
						variantInPop = njh::in(gPos.chromStart_ + indel.second.refPos_, proteinVariants_.at(gPos.chrom_).deletionsFinal) &&
							njh::in(indel.second.gapedSequence_, proteinVariants_.at(gPos.chrom_).deletionsFinal.at(gPos.chromStart_ + indel.second.refPos_)) &&
							proteinVariants_.at(gPos.chrom_).deletionsFinal.at(gPos.chromStart_ + indel.second.refPos_).at(indel.second.gapedSequence_) > 0;
					}
					individualVariantInfoOut << njh::conToStr(toVecStr(
							gPos.chrom_,
							gPos.chromStart_ + indel.second.refPos_ + 1,
							gPos.chromStart_ + indel.second.refPos_ + indel.second.gapedSequence_.size(),
							seqName,
							"deletion",
							indel.second.gapedSequence_,
							std::string(indel.second.gapedSequence_.size(), '-'),
							referenceCodonSeq,
							std::string(indel.second.gapedSequence_.size()*3, '-'),
							njh::boolToStr(false),
							totalVariants,
							genomicLocation.genUIDFromCoords(),
							variantInPop), "\t") << "\n";
				}
			}
		}
	}
}


void TranslatorByAlignment::TranslatorByAlignmentResult::writeSeqLocationsTranslation(std::ostream & out) const {
	auto seqNames = getVectorOfMapKeys(translations_);
	njh::sort(seqNames);
	for(const auto & seqName : seqNames){
		for(const auto & seqs : translations_.at(seqName)){
			out << seqs.second.genBedRec().toDelimStrWithExtra()<< std::endl;
		}
	}
	auto filteredSeqName = getVectorOfMapKeys(translations_);
	njh::sort(filteredSeqName);
	for(const auto & filtName : filteredSeqName){
		if(!njh::in(filtName, filteredSeqName)){
			out << "*"
					<< "\t" << "*"
					<< "\t" << "*"
					<< "\t" << filtName
					<< "\t" << "*"
					<< "\t" << "*" << std::endl;
		}
	}
	for(const auto & pop : iter::sorted(seqsUnableToBeMapped_)){
		out << "*"
				<< "\t" << "*"
				<< "\t" << "*"
				<< "\t" << pop
				<< "\t" << "*"
				<< "\t" << "*" << std::endl;
	}
}


std::set<std::string> TranslatorByAlignment::TranslatorByAlignmentResult::getAllSeqNames() const {
	std::set<std::string> ret;
	njh::addVecToSet(seqsUnableToBeMapped_, ret);
	njh::addVecToSet(seqsTranslationFiltered_, ret);
	njh::addVecToSet(njh::getVecOfMapKeys(translations_), ret);
	return ret;
}


void TranslatorByAlignment::TranslatorByAlignmentResult::writeSeqLocations(std::ostream & out ) const {
	auto seqAlnNames = getVectorOfMapKeys(seqAlns_);
	njh::sort(seqAlnNames);
	for(const auto & seqName : seqAlnNames){
		for(const auto & loc : seqAlns_.at(seqName)){
			out << loc.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
	}
	for(const auto & pop : iter::sorted(seqsUnableToBeMapped_)){
		out << "*"
				<< "\t" << "*"
				<< "\t" << "*"
				<< "\t" << pop
				<< "\t" << "*"
				<< "\t" << "*" << std::endl;
	}
}


std::set<uint32_t> TranslatorByAlignment::getAllInterestingAAPosZeroBased(const std::string & transcript, const TranslatorByAlignmentResult & results){
	std::set<uint32_t> ret;
	//add in known locations
	//known mutation locations are one-based positioned
	if(njh::in(transcript, knownAminoAcidPositions_)){
		for(const auto pos : knownAminoAcidPositions_[transcript]){
			ret.emplace(pos - 1);
		}
	}
	//add in locations with variations in final set
	if(njh::in(transcript, results.proteinVariants_)){
		njh::addVecToSet(njh::getVecOfMapKeys(results.proteinVariants_.at(transcript).snpsFinal), ret);
	}
	return ret;
}



TranslatorByAlignment::TranslatorByAlignmentResult TranslatorByAlignment::run(
		const SeqIOOptions & seqOpts,
		const std::unordered_map<std::string, std::unordered_set<std::string>> & sampCountsForHaps,
		const RunPars & rPars){
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
	TranslatorByAlignmentResult ret;
	std::vector<bfs::path> fnpsToRemove;
	auto seqInputFnp = njh::files::make_path(pars_.workingDirtory_, "inputSeqs.fasta");
	fnpsToRemove.emplace_back(seqInputFnp);
	uint64_t seqMaxLen = 500;
	uint32_t averageLen = 0;
	VecStr names;
	{

		//write out fasta file of input, replacing seq name so can be mapped
		seqInfo seq;
		SeqInput reader(seqOpts);
		reader.openIn();
		SeqOutput writer(SeqIOOptions::genFastaOut(seqInputFnp));
		writer.openOut();
		uint32_t pos = 0;
		std::vector<uint32_t> allReadLens;
		while(reader.readNextRead(seq)){
			allReadLens.emplace_back(len(seq));
			readVec::getMaxLength(seq, seqMaxLen);
			names.emplace_back(seq.name_);
			seq.name_ = estd::to_string(pos);
			++pos;
			writer.write(seq);
		}
		averageLen = static_cast<uint32_t>(std::round(vectorMean(allReadLens)));
	}
	seqMaxLen = seqMaxLen + rPars.realnPars.extendAmount * 2;
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;

	BioCmdsUtils bRunner(false);
	auto twoBitConversionOutput = bRunner.RunFaToTwoBit(pars_.lzPars_.genomeFnp);
	// BioCmdsUtils::checkRunOutThrow(twoBitConversionOutput, __PRETTY_FUNCTION__);
	if(!pars_.useLastz_){
		auto bowtie2RunOutput = bRunner.RunBowtie2Index(pars_.lzPars_.genomeFnp);
		// BioCmdsUtils::checkRunOutThrow(bowtie2RunOutput, __PRETTY_FUNCTION__);
	}

	auto gprefix = bfs::path(pars_.lzPars_.genomeFnp).replace_extension("");
	auto twoBitFnp = gprefix.string() + ".2bit";

	TwoBit::TwoBitFile tReader(twoBitFnp);

	auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(seqInputFnp);

	{

		uniqueSeqInOpts.out_.outFilename_ = njh::files::make_path(pars_.workingDirtory_, "aligned_inputSeqs.sorted.bam");
		uniqueSeqInOpts.out_.outExtention_ = ".sorted.bam";
		fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_);
		fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_.string() + ".bai");

		uniqueSeqInOpts.out_.transferOverwriteOpts(seqOpts.out_);
		//map to genome
		if(!pars_.useLastz_){
			auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, pars_.additionalBowtieArguments_);
			//auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 --end-to-end");
			BioCmdsUtils::checkRunOutThrow(bowtieRunOut, __PRETTY_FUNCTION__);
		}else{
			auto lastzRunOut = bRunner.lastzAlign(uniqueSeqInOpts, pars_.lzPars_);
			BioCmdsUtils::checkRunOutThrow(lastzRunOut, __PRETTY_FUNCTION__);
		}
		// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//count the regions mapped
		auto regionsCounter = GenomicRegionCounter::countRegionsInBam(uniqueSeqInOpts.out_.outName());

		auto ids = regionsCounter.getIntersectingGffIds(pars_.gffFnp_);
		ret.geneIds_ = ids;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "ids.size() : " << ids.size() << std::endl;

		// get gene information
		auto geneInfoDir = njh::files::make_path(pars_.workingDirtory_, "geneInfos");
		if(pars_.keepTemporaryFiles_ || pars_.writeOutGeneInfos_){
			njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
		}
		OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));

		std::unordered_map<std::string, VecStr> idToTranscriptName;
		std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> rawGenes = GeneFromGffs::getGenesFromGffForIds(pars_.gffFnp_, ids);
//		std::cout << "rawGenes.size(): " << rawGenes.size() << std::endl;

		std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;


		for(const auto & gene : rawGenes){
			bool failFilter = false;
			for(const auto & transcript : gene.second->mRNAs_){
				if(njh::in(njh::strToLowerRet(transcript->type_), VecStr{"rrna", "trna", "snorna","snrna","ncrna"}) ){
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					transcript->writeGffRecord(std::cout);
					failFilter = true;
					break;
				}
			}

			if(!failFilter){
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				genes[gene.first] = gene.second;
			}
		}
//		std::cout << "genes.size(): " << genes.size() << std::endl;

		//std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>> geneTranscriptInfos;
		uint64_t proteinMaxLen = seqMaxLen;
		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;
		// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(const auto & gene : genes){
			genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
			for(const auto & transcript : gene.second->mRNAs_){
				idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
			}
			ret.transcriptInfosForGene_[gene.first] = gene.second->generateGeneSeqInfo(tReader, false);
			for(const auto & transcriptInfo : ret.transcriptInfosForGene_[gene.first]){
				ret.translationInfoForTranscirpt_[transcriptInfo.first] = transcriptInfo.second;
				readVec::getMaxLength(transcriptInfo.second->protein_, proteinMaxLen);
	//			ret.proteinForTranscript_[transcriptInfo.first] = transcriptInfo.second->protein_.seq_;
				ret.proteinVariants_.emplace(transcriptInfo.first,
										VariantsInfo{Bed3RecordCore(transcriptInfo.first, 0, transcriptInfo.second->protein_.seq_.size()), transcriptInfo.second->protein_});
			}
			if(pars_.keepTemporaryFiles_ || pars_.writeOutGeneInfos_){
				gene.second->writeOutGeneInfo(tReader, outOpts);
			}
		}
		// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		aligner alignObj(proteinMaxLen, gapScoringParameters(6,1,0,0,0,0), substituteMatrix(2,-2));
		//aligner alignObjSeq(seqMaxLen + rPars.realnPars.extendAmount * 2, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2));

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
		// //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		BamTools::BamReader bReader;
		bReader.Open(uniqueSeqInOpts.out_.outName().string());
		checkBamOpenThrow(bReader, uniqueSeqInOpts.out_.outName());
		auto refData = bReader.GetReferenceData();
		BamTools::BamAlignment bAln;
		auto chromLengths = tReader.getSeqLens();

		struct MinMaxPos{
			MinMaxPos()= default;
			size_t minPos_{std::numeric_limits<uint32_t>::max()};
			size_t maxPos_{0};
		};
//		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::unordered_map<std::string, MinMaxPos> minMaxPositionsPerChrom;
		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;cDNAIntersectedWith.size
				bAln.Name = names[njh::StrToNumConverter::stoToNum<uint32_t>(bAln.Name)];
				auto balnGenomicRegion = GenomicRegion(bAln, refData);
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				minMaxPositionsPerChrom[balnGenomicRegion.chrom_].minPos_ = std::min(minMaxPositionsPerChrom[balnGenomicRegion.chrom_].minPos_,balnGenomicRegion.start_);
				minMaxPositionsPerChrom[balnGenomicRegion.chrom_].maxPos_ = std::max(minMaxPositionsPerChrom[balnGenomicRegion.chrom_].maxPos_,balnGenomicRegion.end_);
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << "rPars.realnPars.extendAmount: " << rPars.realnPars.extendAmount << std::endl;
//				std::cout << "averageLen: " << averageLen << std::endl;
//				std::cout << "bAln.QueryBases.size(): " << bAln.QueryBases.size() << std::endl;
//				std::cout << "bAln.GetEndPosition() - bAln.Position: " << bAln.GetEndPosition() - bAln.Position << std::endl;
//				std::cout << "diff: " << uAbsdiff(averageLen, bAln.QueryBases.size()) << std::endl;
//				std::cout << "diff with aln: " << uAbsdiff(averageLen, bAln.GetEndPosition() - bAln.Position) << std::endl;
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				auto reAlignParsCopy = rPars.realnPars;
				if(uAbsdiff(averageLen, bAln.GetEndPosition() - bAln.Position) > reAlignParsCopy.extendAmount){
					reAlignParsCopy.extendAmount = reAlignParsCopy.extendAmount + uAbsdiff(averageLen, bAln.GetEndPosition() - bAln.Position);
				}
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				ReAlignedSeq results;
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				if(uAbsdiff(averageLen, bAln.GetEndPosition() - bAln.Position) > 100){
//					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
					//adjust aligner giving more weight to matches to help with large deletions or insertions
					aligner alignObjAdjusted(proteinMaxLen, gapScoringParameters(6,1,0,0,0,0), substituteMatrix(10,-2));
					results = ReAlignedSeq::genRealignment(bAln, refData, alignObjAdjusted, chromLengths, tReader, reAlignParsCopy);
//					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}else{
//					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << njh::bashCT::red;
//					std::cout << "alignObj.parts_.maxSize_: " << alignObj.parts_.maxSize_ << std::endl;
//					std::cout << njh::bashCT::reset;
					results = ReAlignedSeq::genRealignment(bAln, refData, alignObj, chromLengths, tReader, reAlignParsCopy);
//					std::cout << results.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//					if(bAln.Name == "Pf13-2840639-2840945.398[HapPopUIDCount=1]"){
//						exit(1);
//					}
//					 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}
//				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				minMaxPositionsPerChrom[balnGenomicRegion.chrom_].minPos_ = std::min(minMaxPositionsPerChrom[results.gRegion_.chrom_].minPos_,results.gRegion_.start_);
				minMaxPositionsPerChrom[balnGenomicRegion.chrom_].maxPos_ = std::max(minMaxPositionsPerChrom[results.gRegion_.chrom_].maxPos_,results.gRegion_.end_);
				ret.seqAlns_[bAln.Name].emplace_back(results);

				if (!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)) {
					continue;
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				for (const auto & g : njh::mapAt(alnRegionToGeneIds,balnGenomicRegion.createUidFromCoords())) {
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					const auto & currentGene = njh::mapAt(genes, g);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					const auto & currentGeneInfo = njh::mapAt(ret.transcriptInfosForGene_, g);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> translations;
					try {
//            std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
						translations = translateBasedOnAlignment(results, *currentGene, currentGeneInfo, alignObj);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//            std::cout << "translations.size(): " << translations.size() << std::endl;
						for(const auto & trans : translations){
							auto queryTransStart = trans.second.queryAlnTranslation_.seq_.find_first_not_of('-');
							if('-' == trans.second.refAlnTranslation_.seq_[queryTransStart]){
								//probably should do a more intensive check here fo
								ret.filteredOffTranslations_[bAln.Name].emplace(trans);
//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
							}else if(countOccurences(trans.second.queryAlnTranslation_.seq_, "*") > pars_.allowableStopCodons_){
//                trans.second.queryAlnTranslation_.outPutSeq(std::cout);
								ret.filteredOffTranslations_[bAln.Name].emplace(trans);
//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
							}else{
								ret.translations_[bAln.Name].emplace(trans);
//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
							}
						}
					} catch (std::exception & e) {
//						std::cerr << e.what() << std::endl;
						//Generally if an exception happen while attempting to translate alignment is mangled
						if(!njh::in(bAln.Name, ret.translations_)){
							ret.seqsTranslationFiltered_.emplace_back(bAln.Name);
						}
					}
				}
			}else if(!bAln.IsMapped()){
				ret.seqsUnableToBeMapped_.emplace_back(names[njh::StrToNumConverter::stoToNum<uint32_t>(bAln.Name)]);
			}
		}
		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//list the seqs no tran
		for(const auto & filteredOff : ret.filteredOffTranslations_){
			if(!njh::in(filteredOff.first, ret.translations_)){
				ret.seqsTranslationFiltered_.emplace_back(filteredOff.first);
			}
		}



		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//index snps
		for(const auto & positons : minMaxPositionsPerChrom){
			Bed3RecordCore chromRegion(
					positons.first,
					positons.second.minPos_,
					positons.second.maxPos_);
			auto refSeq = GenomicRegion(chromRegion).extractSeq(tReader);
			ret.seqVariants_.emplace(chromRegion.chrom_, VariantsInfo(chromRegion, refSeq));
	//		for(uint32_t seqPos = 0; seqPos < refSeq.seq_.size(); ++seqPos){
	//			ret.baseForPosition_[chromRegion.chrom_][chromRegion.chromStart_ + seqPos] = refSeq.seq_[seqPos];
	//		}
		}
		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(const auto & seqName : ret.seqAlns_){
			for(const auto & aln : seqName.second){
				uint32_t popCount = njh::mapAt(sampCountsForHaps, aln.querySeq_.name_).size();
				ret.seqVariants_.at(aln.gRegion_.chrom_).addVariantInfo(
						aln.alnRefSeq_.seq_,
						aln.alnQuerySeq_.seq_,
						popCount,
						njh::mapAt(sampCountsForHaps, aln.querySeq_.name_),
						aln.comp_,
						aln.gRegion_.start_);
			}
		}
		//set finals for the snps
		for(auto & varPerChrom : ret.seqVariants_){
			varPerChrom.second.setFinals(rPars);
		}


		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//index amino acid changes per transcript
		for(const auto & seqName : ret.translations_){
			for(const auto & transcript : seqName.second){
				if(countOccurences(transcript.second.queryAlnTranslation_.seq_, "*") > 1){
					//should log which ones have messed up translations
				}else{
					auto popCount = njh::mapAt(sampCountsForHaps, seqName.first).size();
					ret.proteinVariants_.at(transcript.first).addVariantInfo(
							transcript.second.refAlnTranslation_.seq_,
							transcript.second.queryAlnTranslation_.seq_,
							popCount,
							njh::mapAt(sampCountsForHaps, seqName.first),
							transcript.second.comp_,
							0);
				}
			}
		}
		 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(auto & varPerTrans : ret.proteinVariants_){
			varPerTrans.second.setFinals(rPars);
		}
		//type sequences for significant protein variation including codon info
		for(auto & varPerTrans : ret.proteinVariants_){

			std::set<uint32_t> knownMutationsLocationsZeroBased;
			for(const auto pos : knownAminoAcidPositions_[varPerTrans.first]){
				knownMutationsLocationsZeroBased.emplace(pos - 1);
			}
			//this includes both known mutations as well as locations with significant variation
			std::set<uint32_t> allLocations = getAllInterestingAAPosZeroBased(varPerTrans.first, ret);
			for (auto & seqName : ret.translations_) {
				if (njh::in(varPerTrans.first, seqName.second)) {
					for (const auto & loc : allLocations) {
						if(loc < std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ || loc > std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_){
							//location is not within the aligned translation
							continue;
						}
						auto codon = seqName.second[varPerTrans.first].getCodonForAARefPos(loc);
						ret.fullAATypedWithCodonInfo_[seqName.first].emplace_back(
								TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
										njh::in(loc, knownMutationsLocationsZeroBased)));
						ret.translated_fullAATypedWithCodonInfo_[njh::pasteAsStr(seqName.first, "[transcript=", varPerTrans.first, "]")].emplace_back(
								TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
										njh::in(loc, knownMutationsLocationsZeroBased)));
						if(njh::in(loc, varPerTrans.second.snpsFinal)){
							ret.variantAATypedWithCodonInfo_[seqName.first].emplace_back(
									TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
											njh::in(loc, knownMutationsLocationsZeroBased)));
							ret.translated_variantAATypedWithCodonInfo_[njh::pasteAsStr(seqName.first, "[transcript=", varPerTrans.first, "]")].emplace_back(
									TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
											njh::in(loc, knownMutationsLocationsZeroBased)));
						}
					}
				}
			}
		}

		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> translated_fullAATypedWithCodonInfo_;
		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> translated_variantAATypedWithCodonInfo_;

		//sort
		for(auto & seqName : ret.fullAATypedWithCodonInfo_){
			njh::sort(seqName.second, [](const TranslatorByAlignment::AAInfo & info1, const TranslatorByAlignment::AAInfo & info2){
				if(info1.transcriptName_ == info2.transcriptName_){
					return info1.zeroBasedPos_ < info2.zeroBasedPos_;
				}else{
					return info1.transcriptName_ < info2.transcriptName_;
				}
			});
		}
	}
	fnpsToRemove.emplace_back(seqInputFnp);
	//remove the temporary files
	if(!pars_.keepTemporaryFiles_){
		for(const auto & fnp : fnpsToRemove){
			if(bfs::exists(fnp)){
				bfs::remove(fnp);
			}
		}
	}
//	//remove input seqs
//	if(!pars_.keepTemporaryFiles_){
//		bfs::remove(seqInputFnp);
//		bfs::remove(uniqueSeqInOpts.out_.outName());
//		bfs::remove(uniqueSeqInOpts.out_.outName().string() + ".bai");
//	}
	 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return ret;
}


std::unordered_map<std::string, std::set<uint32_t>> TranslatorByAlignment::readInAAPositions(const bfs::path & knownAminoAcidChangesFnp){
	njh::files::checkExistenceThrow(knownAminoAcidChangesFnp, __PRETTY_FUNCTION__);
	table knownAminoAcidChanges(knownAminoAcidChangesFnp, "\t", true);
	njh::for_each(knownAminoAcidChanges.columnNames_, [](std::string & col){
		njh::strToLower(col);
	});
	knownAminoAcidChanges.setColNamePositions();
	knownAminoAcidChanges.checkForColumnsThrow(VecStr{"transcriptid", "aaposition"}, __PRETTY_FUNCTION__);
	std::unordered_map<std::string, std::set<uint32_t>> knownMutationsLocationsMap;
	if(knownAminoAcidChanges.nRow() > 0){
		for(const auto & row : knownAminoAcidChanges){
			if(std::all_of(row.begin(), row.end(), [](const std::string & element){
				return "" ==element;
			})){
				continue;
			}
			knownMutationsLocationsMap[row[knownAminoAcidChanges.getColPos("transcriptid")]].emplace(njh::StrToNumConverter::stoToNum<uint32_t>(row[knownAminoAcidChanges.getColPos("aaposition")]));
		}
	}
	return knownMutationsLocationsMap;
}



}  // namespace njhseq
