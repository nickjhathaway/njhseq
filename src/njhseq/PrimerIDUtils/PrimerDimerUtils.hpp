#pragma once

// Created by Nicholas Hathaway on 5/23/23.
/* 
    
*/

#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"


namespace njhseq {

class PrimerDimerUtils{
public:
	static std::unordered_map<std::string, double> matchScores;
	static std::unordered_map<std::string, double> singleMismatchScores;


	double doubleMismatchScore_ = 0.2;
	double endBonus_ = -0.5;
	uint32_t endBonusLen_ = 4;
	std::unordered_map<std::string, double> scoringMatrix_;
	bool debug_{false};

	PrimerDimerUtils();

private:
	[[nodiscard]] double computeDimerScore(const std::string &shorter, const std::string &longer) const;
public:
	[[nodiscard]] double computeDimerScoreTop(const std::string &str1, const std::string &str2) const;

	[[nodiscard]] std::vector<std::vector<double>> computeFullScoreMatrix(const std::vector<seqInfo> & primers) const;

	static void writeMatrix(const std::vector<std::vector<double>> & scores, std::ostream & out, const std::vector<seqInfo> & primers, bool outputRawMatrix = false);
};



}  // namespace njhseq

