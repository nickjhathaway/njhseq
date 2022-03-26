#pragma once

/*
 * nhmmscanOutput.hpp
 *
 *  Created on: Feb 19, 2022
 *      Author: nick
 */



#include "njhseq/objects/BioDataObject/HmmerObjs/HmmerSeqHitTab.hpp"

#include "njhseq/IO/IOOptions/OutOptions.hpp"

namespace njhseq {

class nhmmscanOutput{
public:

	VecStr header_;
	std::unordered_map<std::string, std::string> parameterInfo_;

	class Hit : public HmmerSeqHitTab {
	public:
		// Hit example
		/*
		 * >> stevors3d7
    score  bias    Evalue   hmmfrom    hmm to     alifrom    ali to      envfrom    env to      mod len      acc
   ------ ----- ---------   -------   -------    --------- ---------    --------- ---------    ---------    ----
 !   51.6   0.0   1.9e-16        37       257 ..      3238      3007 ..      3258      2999 ..      1076    0.53

  Alignment:
  score: 51.6 bits
  stevors3d7   37 tttttttgaAttaAaAatgaatatatAt.taccttAAaatgtTattatttacctttttaataaatacatTattattaccacAttaTgtAtgtaaaAaaa.............aTat........AttattattacaataacataaataTAtaAtaaaaacaattattacAtttt....caatatatatt
cattaattc...ttttttttttttttTtTAggaaaAttatctaaataaccattaTAAtatA 257
                  t ttt   +At aAaAatg   atatAt ta+ t A   t tT  tat ta   tt t +ta ata+at a tatta+   AttaT        aAaa              aTa+        At at at acaata  at aat TA aA  aaa   at at a+At tt    ca+t  a  tt
catt+ t +   t tt tt t tt tTtTA     Attat  a ataa  a t TAA atA
           0 3238 TATTTACCTATAAAAAATGGTGATATATaTATGTGA---TATTTATATGTAGA-TTATTTTATATATATAAATATTATATTATTATA-------AAAATatttattagatccATACctgcgatgATGATAATAACAATATAATTAATGTAAAAG-AAA---ATAATAATATGTTatagCAGTGGACTTT
CATTTTTATtgtTATTCTTATATTATTTTA---TGATTATACATATAATTAATTTAAAATA 3007
                  457777888999********9999999976777654...55555565555543.444556666777777777777766666666653.......222220111111121122222111111111222233333333333333333333333332.333...4445555555554222256666777777
777777666433444555556666677776...5677778878777777766666543333 PP
		 */
//		std::string modelName_;
//		double score_{std::numeric_limits<double>::max()};
//		double bias_{std::numeric_limits<double>::max()};
//		double Evalue_{std::numeric_limits<double>::max()};
//		uint32_t hmmfrom_{std::numeric_limits<uint32_t>::max()};
//		uint32_t hmmto_{std::numeric_limits<uint32_t>::max()};
		std::string hmmEdgeInfo_; //.., [., .], or []
//		uint32_t alifrom_{std::numeric_limits<uint32_t>::max()};
//		uint32_t alito_{std::numeric_limits<uint32_t>::max()};
		std::string aliEdgeInfo_; //.., [., .], or []
//		uint32_t envfrom_{std::numeric_limits<uint32_t>::max()};
//		uint32_t envto_{std::numeric_limits<uint32_t>::max()};
		std::string envEdgeInfo_; //.., [., .], or []
//		uint32_t modlen_{std::numeric_limits<uint32_t>::max()};
		double acc_{std::numeric_limits<double>::max()};//ranges from 0 to 100


		std::string modelAln_;
		std::string alnAgreement_;
		std::string queryAln_;
		std::string aln_posterior_probability_;
		/*
		 * posterior probability (essentially the expected accuracy) of each aligned residue. A 0 means 0-5%, 1 means 5-15%, and so on; 9 means 85-95%, and a * means 95-100% posterior probability
		 */

		[[nodiscard]] uint32_t sumOfPosteriorProbability() const;
		[[nodiscard]] double averagePP() const;

		[[nodiscard]] double percentPerfectHit() const;

		[[nodiscard]] double percentGappedHit() const;

		bool overlaps_env(const Hit & otherHit, uint32_t minOverlap = 1) const;

		static VecStr getOutputDetHeader();
		VecStr getOutputDet() const;
		virtual Json::Value toJson() const;
		virtual ~Hit(){

		}
	};
	class QueryResults{
	public:
		std::string queryName_;
		uint64_t queryLen_{std::numeric_limits<uint64_t>::max()};

		std::vector<Hit> hits_;

		uint32_t targetModels_{std::numeric_limits<uint32_t>::max()};
		uint64_t targetModNodes_{std::numeric_limits<uint64_t>::max()};
		uint64_t residuesSearched_{std::numeric_limits<uint64_t>::max()};
		uint64_t residuesPass_SSV_filter_{std::numeric_limits<uint64_t>::max()};
		uint64_t residuesPass_bias_filter_{std::numeric_limits<uint64_t>::max()};
		uint64_t residuesPass_Vit_filter_{std::numeric_limits<uint64_t>::max()};
		uint64_t residuesPass_Fwd_filter_{std::numeric_limits<uint64_t>::max()};

		std::string cpuRunInfo_;
		double Mc_per_sec_{std::numeric_limits<double>::max()};
		double hitsFrac_{std::numeric_limits<double>::max()};

		[[nodiscard]] Json::Value toJson() const;

		static void sortHitsByEvaluesScores(std::vector<Hit> &hits);

		static std::vector<Hit> getNonOverlapHits(const std::vector<Hit> &hits );
		static std::vector<Hit> getNonOverlapHits(std::vector<Hit> &hits, const std::function<bool(const Hit&, const Hit&)> & sortFunc);

	};

	std::vector<QueryResults> qResults_;


	[[nodiscard]] Json::Value toJson() const;


	static nhmmscanOutput parseRawOutput(const bfs::path & input);
	static nhmmscanOutput parseRawOutput(const bfs::path & input, const std::unordered_map<uint32_t, std::string> & seqKey);


	static bool run_hmmpress_ifNeed(const bfs::path & hmmModelFnp);

	void outputCustomHitsTable(const OutOptions & outOpts) const;

	/**
	 * @brief rename query, only works if the input was named with query index
	 * @param seqKey index to name
	 */
	void renameQuery(const std::unordered_map<uint32_t, std::string> & seqKey);

	struct PostProcessHitsPars{
		uint32_t minLength = 0;
		double accCutOff = 0;
		double scoreCutOff = 0;
		uint32_t hmmStartFilter = std::numeric_limits<uint32_t>::max();
	};
	struct PostProcessHitsRes{
		std::unordered_map<std::string, std::vector<nhmmscanOutput::Hit>> filteredHitsByQuery_;
		std::unordered_map<std::string, std::vector<nhmmscanOutput::Hit>> filteredNonOverlapHitsByQuery_;
	};

	PostProcessHitsRes postProcessHits(const PostProcessHitsPars & pars);

	
	void writeInfoFiles(const PostProcessHitsRes & postProcessResults, const bfs::path & outputDir) const;
};



}  // namespace njhseq



