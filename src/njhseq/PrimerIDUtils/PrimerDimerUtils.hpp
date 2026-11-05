#pragma once

// Created by Nicholas Hathaway on 5/23/23.
/* 
    
*/

#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

#include <array>
#include <functional>


namespace njhseq {

class PrimerDimerUtils{
public:
	static std::unordered_map<std::string, double> matchScores;
	static std::unordered_map<std::string, double> singleMismatchScores;


	double doubleMismatchScore_ = 0.2;
	double endBonus_ = -0.5;
	uint32_t endBonusLen_ = 4;
	std::unordered_map<std::string, double> scoringMatrix_;
	uint32_t proofreadingMaxUnpaired3Prime_ = 3; //how far the stacked stretch can sit from the nearer 3' end
	bool debug_{false};

	PrimerDimerUtils();

private:
	std::array<double, 256> proofreadingDoubletScores_{};

	[[nodiscard]] double computeDimerScore(const std::string &shorter, const std::string &longer) const;
	[[nodiscard]] std::vector<std::vector<double>> computeFullMatrix(const std::vector<seqInfo> & primers, uint32_t numThreads,
			const std::function<double(const std::string &, const std::string &)> & scorer) const;
public:
	[[nodiscard]] double computeDimerScoreTop(const std::string &str1, const std::string &str2) const;

	/**@brief A score for polymerases that can remove unpaired 3' bases before extending, so a 3' end doesn't have to be paired
	 *
	 * The most negative contiguous run of nearest neighbour stacks over all registers, counting a register only when its
	 * best run is within proofreadingMaxUnpaired3Prime_ bases of a 3' end. 0 if no register qualifies.
	 */
	[[nodiscard]] double computeProofreadingDimerScore(const std::string &str1, const std::string &str2) const;

	[[nodiscard]] std::vector<std::vector<double>> computeFullScoreMatrix(const std::vector<seqInfo> & primers, uint32_t numThreads = 1) const;

	[[nodiscard]] std::vector<std::vector<double>> computeFullProofreadingScoreMatrix(const std::vector<seqInfo> & primers, uint32_t numThreads = 1) const;

	static void writeMatrix(const std::vector<std::vector<double>> & scores, std::ostream & out, const std::vector<seqInfo> & primers, bool outputRawMatrix = false);
};



}  // namespace njhseq

