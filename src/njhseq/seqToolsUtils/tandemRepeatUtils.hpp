#pragma once

/*
 * tandemRepeatUtils.hpp
 *
 *  Created on: May 31, 2020
 *      Author: nick
 */



#include "njhseq/programUtils/seqSetUp.hpp"
#include "njhseq/objects/helperObjects/motif.hpp"
#include "njhseq/IO/SeqIO/SeqInput.hpp"
#include "njhseq/IO/OutputStream.hpp"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"




namespace njhseq {

class SimpleTandemRepeatFinder {
public:
	struct RepeatStartLenRepeats {
		RepeatStartLenRepeats(const uint32_t start, const uint32_t len, const uint32_t repeats): start_(start), len_(len),
			repeats_(repeats) {
		}

		uint32_t start_;
		uint32_t len_;
		uint32_t repeats_;
	};
	static RepeatStartLenRepeats getRepeatInfo(const std::string & seq, const motif & m);

	struct SimpleTRFinderLocsPars {
		std::vector<char> alphabet{'A', 'G', 'C', 'T'};
		OutOptions outOpts{bfs::path(""), ".bed"};
		uint32_t minRepeatUnitSize = 1;
		uint32_t maxRepeatUnitSize = 3;
		uint32_t lengthCutOff = 12;
		uint32_t minNumRepeats = 2;
		bool searchAllUnits = false;
		bool doNotAddFlankingSeq = false;
		bool addFullSeqToOuput = false;
		uint32_t numThreads = 1;

		bool verbose = false;
		bool debug = false;

		void setDefaultOpts(seqSetUp& setUp);
	};
	explicit SimpleTandemRepeatFinder(SimpleTRFinderLocsPars  pars);
	SimpleTRFinderLocsPars pars_;

	void runSimpleTRFinderLocs(const SeqIOOptions& seqInput) const;

	[[nodiscard]] VecStr genAllUnitsPossible() const;
	struct MinimalUnitsAndAltMotifs {
		std::shared_ptr<VecStr> allUnits;
		std::unordered_map<std::string, std::vector<motif>> altMots;
	};

	[[nodiscard]] MinimalUnitsAndAltMotifs genMinimalUnitsNeededForSearch() const;
};





}  // namespace njhseq
