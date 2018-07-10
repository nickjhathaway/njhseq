/*
 * GeneFromGffs.cpp
 *
 *  Created on: Aug 7, 2017
 *      Author: nick
 */

#include "GeneFromGffs.hpp"
#include "bibseq/objects/BioDataObject/reading.hpp"
#include "bibseq/objects/BioDataObject/BioDataFileIO.hpp"

namespace bibseq {

GeneFromGffs::GeneFromGffs(const std::vector<std::shared_ptr<GFFCore>> & geneRecords){
	VecStr alreadyAddedIDs;
	//just from what i've seen from chado annotations
	VecStr acceptableMrnaLikeRecords = {"mRNA", "ncRNA", "rRNA", "snoRNA", "tRNA", "snRNA"};
	for(const auto & record : geneRecords){
		if(record->hasAttr("ID") && bib::in(record->getAttr("ID"), alreadyAddedIDs)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error already have record with ID "
					<< record->getAttr("ID")  << "\n";
			throw std::runtime_error{ss.str()};
		}
		if("gene" == record->type_ ){
			if (nullptr != gene_) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__
						<< ", error already have gene record with name "
						<< gene_->getAttr("ID")
						<< ", found another record with feature equaled to gene, named: "
						<< record->getAttr("ID") << "\n";
				throw std::runtime_error{ss.str()};
			}
			gene_ = record;
		}else if(bib::in(record->type_, acceptableMrnaLikeRecords)){
			mRNAs_.emplace_back(record);
		}else if("CDS" == record->type_){
			CDS_[record->getAttr("Parent")].emplace_back(record);
		}else if("exon" == record->type_){
			exons_[record->getAttr("Parent")].emplace_back(record);
		}else if("polypeptide" == record->type_){
			polypeptides_[record->getAttr("Derives_from")].emplace_back(record);
		}else {
			others_[record->getAttr("Parent")].emplace_back(record);
		}
	}
	std::stringstream missingMessage;
	bool failed = false;
	//check that CDS and exons have parents that are in mRNAs
	for(const auto & CD : CDS_){
		bool found = false;
		for(const auto & mRNA : mRNAs_){
			if(mRNA->getAttr("ID") == CD.first){
				found = true;
				break;
			}
		}
		if(!found){
			missingMessage << "CDS_ has records for " << CD.first << " but not fond in mRNAs_" << "\n";
			missingMessage << "options: " << bib::conToStr(bib::convert<std::shared_ptr<GFFCore>,std::string>(mRNAs_, [](const std::shared_ptr<GFFCore> & gffRecord){
				return gffRecord->getAttr("ID");
			}), ",") << "\n";
			failed = true;
		}
	}
	for(const auto & exon : exons_){
		bool found = false;
		for(const auto & mRNA : mRNAs_){
			if(mRNA->getAttr("ID") == exon.first){
				found = true;
				break;
			}
		}
		if(!found){
			missingMessage << "exons_ has records for " << exon.first << " but not fond in mRNAs_" << "\n";
			missingMessage << "options: " << bib::conToStr(bib::convert<std::shared_ptr<GFFCore>,std::string>(mRNAs_, [](const std::shared_ptr<GFFCore> & gffRecord){
				return gffRecord->getAttr("ID");
			}), ",") << "\n";
			failed = true;
		}
	}
	//check to see that mRNAs have CDS and exon records
	//if there are no exons, then create them from the CDS
	for(const auto & mRNA : mRNAs_){
		if (!bib::in(mRNA->getAttr("ID"), CDS_)
				&& !bib::in(mRNA->getAttr("ID"), exons_)) {
			missingMessage << mRNA->getAttr("ID") << " missing from CDS_ and exons_, gene records should have either CDS or exon features" << "\n";
			missingMessage << "options: " << bib::conToStr(bib::getVecOfMapKeys(CDS_), ",") << "\n";
			failed = true;
		} else if (!bib::in(mRNA->getAttr("ID"), CDS_)
				&& bib::in(mRNA->getAttr("ID"), exons_)) {
			//this is a lazy way of doing this but most of the other functions assume CDS_ always has records
			CDS_[mRNA->getAttr("ID")] = exons_[mRNA->getAttr("ID")];
		}
		//apparently some gff files only have CDS and not exons so the bellow is now commented out
//		if(!bib::in(mRNA->getAttr("ID"), exons_)){
//			missingMessage << mRNA->getAttr("ID") << " missing from exons_" << "\n";
//			missingMessage << "options: " << bib::conToStr(bib::getVecOfMapKeys(exons_), ",") << "\n";
//			failed = true;
//		}
	}

	if(failed){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error missing records: " << "\n";
		ss << missingMessage.str();
		throw std::runtime_error{ss.str()};
	}
}



void GeneFromGffs::writeGffRecords(std::ostream & out) const{
	gene_->writeGffRecord(out);
	auto writeOutRecs = [&out](const std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> & records, const std::string & id){
		if(bib::in(id, records)){
			for(const auto & rec : records.at(id) ){
				rec->writeGffRecord(out);
			}
		}
	};
	for(const auto & rna : mRNAs_){
		rna->writeGffRecord(out);
		auto rnaID = rna->getIDAttr();
		writeOutRecs(CDS_, rnaID);
		writeOutRecs(exons_, rnaID);
		writeOutRecs(polypeptides_, rnaID);
		writeOutRecs(others_, rnaID);
	}
}


std::string GeneFromGffs::getOneGeneDetailedName() const{
	auto allNames = getGeneDetailedName();
	std::set<std::string> uniqueNames;
	for(const auto & name : allNames){
		if("na" != bib::strToLowerRet(name.second)){
			uniqueNames.insert(name.second);
		}
	}
	return bib::conToStr(uniqueNames, ", ");
}

std::unordered_map<std::string, std::string> GeneFromGffs::getGeneDetailedName() const {
	std::unordered_map<std::string, std::string> ret;
	for (const auto & m : mRNAs_) {
		std::string name = "";
		if ("chado" == gene_->source_) {
			if (bib::in(m->getIDAttr(), polypeptides_)) {
				if (polypeptides_.at(m->getIDAttr()).size() == 1
						&& polypeptides_.at(m->getIDAttr()).front()->hasAttr("product")) {
					auto productToks = tokenizeString(polypeptides_.at(m->getIDAttr()).front()->getAttr("product"), ";");
					for(const auto & productTok : productToks){
						auto subToks = tokenizeString(productTok, "=");
						if(2 == subToks.size() && "term" == subToks[0]){
							if("" != name){
								name += ", ";
							}
							name += subToks[1];
						}
					}
				}
			}
		} else if ("PlasmoDB" == gene_->source_) {
			if (gene_->hasAttr("description")) {
				name = gene_->getAttr("description");
			}
		}
		if("" == name){
			name = "NA";
		}
		ret[m->getIDAttr()] = name;
	}

	return ret;
}

uint32_t GeneFromGffs::numOfTranscripts() const {
	return mRNAs_.size();
}


std::unordered_map<std::string, std::vector<Bed6RecordCore>> GeneFromGffs::getIntronExonBedLocs() const{
	auto infoTab = getIntronExonTables();
	std::unordered_map<std::string, std::vector<Bed6RecordCore>> ret;
	for(const auto & transcript : infoTab){
		for(const auto & row : transcript.second){
			std::stringstream line;
			line << bib::conToStr(VecStr{
				row[transcript.second.getColPos("chrom")],
					row[transcript.second.getColPos("start")],
					row[transcript.second.getColPos("end")],
					transcript.first + "-" + row[transcript.second.getColPos("name")],
					row[transcript.second.getColPos("length")],
					row[transcript.second.getColPos("strand")]}, "\t");
			ret[transcript.first].emplace_back(line.str());
		}
	}
	return ret;
}

void GeneFromGffs::writeOutGeneInfo(TwoBit::TwoBitFile & tReader, const OutOptions & outPrefix) const{
	auto gsInfos = generateGeneSeqInfo(tReader, false);

	for(const auto & transcript : mRNAs_){
		//idToTranscriptName[id].emplace_back(transcript->getIDAttr());
		GenomicRegion mRnaRegion(*transcript) ;
		OutOptions transcriptOut(bfs::path(outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + "_basePositions"), ".tab.txt");
		transcriptOut.transferOverwriteOpts(outPrefix);
		auto tOutFile = transcriptOut.openFile();
		auto gsInfo = gsInfos[transcript->getIDAttr()];
		gsInfo->infoTab_.addColumn(VecStr{transcript->getIDAttr()}, "transcript");
		gsInfo->infoTab_.addColumn(VecStr{gene_->getIDAttr()}, "GeneID");
		gsInfo->infoTab_.addColumn(VecStr{std::string(1, transcript->strand_)}, "strand");
		gsInfo->infoTab_.outPutContents(*tOutFile, "\t");
		auto gDNAOpts = SeqIOOptions::genFastaOut(    outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + "_gDNA");
		auto cDNAOpts = SeqIOOptions::genFastaOut(    outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + "_cDNA");
		auto proteinOpts = SeqIOOptions::genFastaOut( outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + "_protein");
		auto tableOpts = TableIOOpts::genTabFileOut(  outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + "_exonIntronPositions", true);
		auto bedOpts = OutOptions(bfs::path(          outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + "_exonIntronPositions.bed"));
		auto transcriptBedOpts = OutOptions(bfs::path(outPrefix.outFilename_.string() + "_" + transcript->getAttr("ID") + ".bed"));
		bedOpts.transferOverwriteOpts(outPrefix);
		transcriptBedOpts.transferOverwriteOpts(outPrefix);
		tableOpts.out_.transferOverwriteOpts(outPrefix);
		OutputStream transcriptBedOut(transcriptBedOpts);
		transcriptBedOut << GenomicRegion(*transcript).genBedRecordCore().toDelimStrWithExtra() <<std::endl;
		auto exonIntronPositions = getIntronExonTables();
		auto exonIntronBeds = getIntronExonBedLocs();
		exonIntronPositions[transcript->getIDAttr()].outPutContents(tableOpts);
		BioDataFileIO<Bed6RecordCore> reader{IoOptions(bedOpts)};
		reader.openWrite(exonIntronBeds[transcript->getIDAttr()], [](const Bed6RecordCore & record, std::ostream & out){
			out << record.toDelimStrWithExtra() << std::endl;
		});
		gDNAOpts.out_.transferOverwriteOpts(outPrefix);
		cDNAOpts.out_.transferOverwriteOpts(outPrefix);
		proteinOpts.out_.transferOverwriteOpts(outPrefix);
		gsInfo->cDna_.name_ = transcript->getAttr("ID") + "_CodingDNA";
		gsInfo->gDna_.name_ = transcript->getAttr("ID") + "_GenomicDNA";
		SeqOutput::write(std::vector<seqInfo>{gsInfo->gDna_}, gDNAOpts);
		SeqOutput::write(std::vector<seqInfo>{gsInfo->cDna_}, cDNAOpts);
		gsInfo->protein_.name_ = transcript->getAttr("ID") + "_protein";
		if('*' == gsInfo->protein_.seq_.back()){
			gsInfo->protein_.trimBack(gsInfo->protein_.seq_.size() - 1);
		}
		SeqOutput::write(std::vector<seqInfo>{gsInfo->protein_}, proteinOpts);
	}
}

std::unordered_map<std::string,std::shared_ptr<GeneSeqInfo>> GeneFromGffs::generateGeneSeqInfo(TwoBit::TwoBitFile & tReader, bool oneBased) const{
	std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> ret;
	for(const auto & transcript : mRNAs_){
		GenomicRegion mRnaRegion(*transcript) ;
		auto gDna = mRnaRegion.extractSeq(tReader);
		auto cdsRegions = gffPtrsToGenomicRegs(CDS_.at(transcript->getAttr("ID")));
		seqInfo cDna;
		seqInfo cDnaAln;
		if(cdsRegions.size() > 1){
			sortGRegionsByStart(cdsRegions);
			if(mRnaRegion.reverseSrand_){
				cDna = cdsRegions.back().extractSeq(tReader);
				cDnaAln = cDna;
				for(const auto & pos : iter::range<uint32_t>(1, cdsRegions.size())){
					//add in gaps
					uint32_t gapSize = cdsRegions[cdsRegions.size() - pos].start_ - cdsRegions[cdsRegions.size() - 1 - pos].start_ - cdsRegions[cdsRegions.size() - 1 - pos].getLen();
					cDnaAln.append(std::string(gapSize, '-'));
					cDna.append(cdsRegions[cdsRegions.size() - 1 - pos].extractSeq(tReader));
					cDnaAln.append(cdsRegions[cdsRegions.size() - 1 - pos].extractSeq(tReader));
				}
			}else{
				cDna = cdsRegions.front().extractSeq(tReader);
				cDnaAln = cDna;
				for(const auto & pos : iter::range<uint32_t>(1, cdsRegions.size())){
					//add in gaps
					uint32_t gapSize = cdsRegions[pos].start_ - cdsRegions[pos - 1].start_ - cdsRegions[pos - 1].getLen();
					cDnaAln.append(std::string(gapSize, '-'));
					cDna.append(cdsRegions[pos].extractSeq(tReader));
					cDnaAln.append(cdsRegions[pos].extractSeq(tReader));
				}
			}
		}else if(cdsRegions.size() == 1){
			cDna = cdsRegions.front().extractSeq(tReader);
			cDnaAln = cDna;
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error no CDS regions found for " << transcript->getAttr("ID") << "\n";
			throw std::runtime_error{ss.str()};
		}
		GeneSeqInfo::GeneSeqInfoPars giInfoPar;
		giInfoPar.oneBasedPos_ = oneBased;
		giInfoPar.region_ = mRnaRegion;
		ret[transcript->getIDAttr()]  = std::make_shared<GeneSeqInfo>(cDna, gDna, giInfoPar);
		ret[transcript->getIDAttr()]->setCDnaAln(cDnaAln);
		ret[transcript->getIDAttr()]->setTable();
	}
	return ret;
}

std::unordered_map<std::string, table> GeneFromGffs::getIntronExonTables() const{
	std::unordered_map<std::string, table> ret;

	for (const auto & transcript : mRNAs_) {
		table transTab(VecStr { "GeneID", "transcript", "type", "chrom", "start",
				"end","name", "length", "strand", "cDNAStart", "cDNAStop" });
		if (CDS_.at(transcript->getAttr("ID")).size() == 1) {
			auto exon1Ptr = CDS_.at(transcript->getAttr("ID")).front();
			auto exon1 = GenomicRegion(*exon1Ptr);
			transTab.addRow(gene_->getIDAttr(),
					transcript->getIDAttr(),
					"exon",
					transcript->seqid_,
					exon1.start_,
					exon1.end_,
					"exon1",
					exon1.getLen(),
					transcript->strand_,
					0,
					exon1.getLen());
		} else {
			auto exonGenRegs = gffPtrsToGenomicRegs(CDS_.at(transcript->getIDAttr()));
			bib::sort(exonGenRegs, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
				return reg1.start_ < reg2.start_;
			});
			struct CDnaPos {
				uint32_t start_;
				uint32_t end_;
			};
			std::unordered_map<std::string, CDnaPos> exonCDNAPositions;
			if(transcript->isReverseStrand()){
				uint32_t exonCount = 1;
				uint32_t cDnaStart = 0;
				for(auto  & exonGenReg : iter::reversed(exonGenRegs)){
					exonGenReg.uid_ = "exon" + estd::to_string(exonCount);
					exonCDNAPositions[exonGenReg.uid_].start_ = cDnaStart;
					exonCDNAPositions[exonGenReg.uid_].end_ = cDnaStart + exonGenReg.getLen();
					cDnaStart = cDnaStart + exonGenReg.getLen();
					++exonCount;
				}
				uint32_t intronCount = 1;
				uint32_t exonTotal = exonGenRegs.size();
				for(const auto pos : iter::range(exonTotal - 1)){
					exonGenRegs.emplace_back(
							"intron" + estd::to_string(intronCount),
							transcript->seqid_,
							exonGenRegs[exonTotal - 1 - pos - 1].end_,
							exonGenRegs[exonTotal - 1 - pos].start_,
							transcript->isReverseStrand());
					++intronCount;
				}
			} else {
				uint32_t exonCount = 1;
				uint32_t cDnaStart = 0;
				for(auto  & exonGenReg : exonGenRegs){
					exonGenReg.uid_ = "exon" + estd::to_string(exonCount);
					exonCDNAPositions[exonGenReg.uid_].start_ = cDnaStart;
					exonCDNAPositions[exonGenReg.uid_].end_ = cDnaStart + exonGenReg.getLen();
					cDnaStart = cDnaStart + exonGenReg.getLen();
					++exonCount;
				}
				uint32_t intronCount = 1;
				for(const auto pos : iter::range(exonGenRegs.size() - 1)){
					exonGenRegs.emplace_back(
							"intron" + estd::to_string(intronCount),
							transcript->seqid_,
							exonGenRegs[pos].end_,
							exonGenRegs[pos + 1].start_,
							transcript->isReverseStrand());
					++intronCount;
				}
			}
			bib::sort(exonGenRegs, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
				return reg1.start_ < reg2.start_;
			});

			for(const auto & reg : exonGenRegs){

				transTab.addRow(gene_->getIDAttr(),
						transcript->getIDAttr(),
						(bib::containsSubString(reg.uid_, "exon") ? "exon" : "intron"),
						transcript->seqid_,
						reg.start_,
						reg.end_,
						reg.uid_,
						reg.getLen(),
						transcript->strand_,
						(bib::containsSubString(reg.uid_, "exon") ? estd::to_string(exonCDNAPositions[reg.uid_].start_) : "NA"),
						(bib::containsSubString(reg.uid_, "exon") ? estd::to_string(exonCDNAPositions[reg.uid_].end_) : "NA"));
			}
		}
		ret[transcript->getIDAttr()] = transTab;
	}
	return ret;
}


std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> GeneFromGffs::getGenesFromGffForIds(
		const bfs::path & gffFnp,
		const std::set<std::string> & ids){

	std::unordered_map<std::string, std::set<std::string>> parents;
	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> gffRecs;
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;
	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(gffFnp))};
	std::vector<std::shared_ptr<GFFCore>> cache;
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	VecStr deirvedFromRecordFeatures {"polypeptide"};
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	while (nullptr != gRecord) {
		if(gRecord->hasAttr("ID") && bib::in(gRecord->getAttr("ID"), ids) ){
			if("gene" != gRecord->type_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error feature type needs to be gene, not " << gRecord->type_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			auto currentId = gRecord->getAttr("ID");
			parents[currentId].insert(gRecord->getAttr("ID"));
			gffRecs[currentId].emplace_back(std::make_shared<GFFCore>(*gRecord));
		}else if(gRecord->hasAttr("Parent")){
			auto currentParent = gRecord->getAttr("Parent");
			for(auto & p : parents){
				if(bib::in(currentParent, p.second)){
					p.second.insert(gRecord->getAttr("ID"));
					gffRecs[p.first].emplace_back(std::make_shared<GFFCore>(*gRecord));
				}
			}
		} else if(bib::in(gRecord->type_, deirvedFromRecordFeatures) ){
			cache.emplace_back(std::make_shared<GFFCore>(*gRecord));
		}
		if("gene" == gRecord->type_){
			for(const auto & fromCache : cache){
				//grab any possible polypeptides or misc records,
				//I have found that sometimes these come before the gene record so you have to keep a cache around since the last cache
				//there might be a better way of doing this
				if(fromCache->hasAttr("Derives_from")){
					auto currentParent = fromCache->getAttr("Derives_from");
					for(auto & p : parents){
						if(bib::in(currentParent, p.second)){
							gffRecs[p.first].emplace_back(std::make_shared<GFFCore>(*fromCache));
						}
					}
				}
			}
			cache.clear();
		}
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (bib::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			bib::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		++count;
	}

	for(const auto & geneGffs : gffRecs){
		genes[geneGffs.first] = std::make_shared<GeneFromGffs>(geneGffs.second);
	}
	return genes;
}

} /* namespace bibseq */
