//
// Created by Nicholas Hathaway on 5/23/23.
//

#include "PrimerDimerUtils.hpp"
#include "njhseq/helpers/seqUtil.hpp"
#include "njhseq/concurrency/PairwisePairFactory.hpp"

namespace njhseq {
PrimerDimerUtils::PrimerDimerUtils() {

	//create scoring matrix
	std::vector<char> alphabet = {'A', 'C', 'G', 'T'};
	std::vector<std::string> doubles;
	doubles.reserve(alphabet.size() * alphabet.size());
	for (const auto &alph1: alphabet) {
		for (const auto &alph2: alphabet) {
			doubles.emplace_back(njh::pasteAsStr(alph1, alph2));
		}
	}
	for (const auto &doub1: doubles) {
		for (const auto &doub2: doubles) {
			scoringMatrix_[njh::pasteAsStr(doub1, "/", doub2)] = doubleMismatchScore_;
		}
	}
	//update match scores
	for (const auto &match: matchScores) {
		scoringMatrix_[match.first] = match.second;
	}
	//update single match scores
	for (const auto &singleMatch: singleMismatchScores) {
		scoringMatrix_[singleMatch.first] = singleMatch.second;
	}
}

double PrimerDimerUtils::computeDimerScore(const std::string &shorter, const std::string &longer) const{
	uint32_t best_start = 0;
	double best_score = 10.0;
	std::vector<bool> best_matches;

//		bool print = shorter == "GTGGATCAACGACGAGCACA" && longer == "ATACACACACACATGCATGC";
	bool print = debug_;
	for (const auto i: iter::range(len(longer) - 1)) {
		std::string matching;
		double currentScore = 0;
		std::vector<bool> matches;
		uint32_t lastJ = 0;
		for (const auto j: iter::range(len(shorter) - 1)) {
			std::string current = longer.substr(i + j, 2) + "/" + shorter.substr(j, 2);
			currentScore += scoringMatrix_.at(current);
			matches.emplace_back(longer[i + j] == shorter[j]);
			if (print) {
				std::cout << "\ti: " << i << std::endl;
				std::cout << "\tj: " << j << std::endl;
				//std::cout << "longer: " << longer << std::endl;
				std::cout << "\tl_bases: " << longer.substr(i + j, 2) << std::endl;
				std::cout << "\ts_bases: " << shorter.substr(j, 2) << std::endl;
				std::cout << "\tcurrent: " << current << std::endl;
				std::cout << "\tscore: " << scoringMatrix_.at(current) << std::endl;
				std::cout << "\ti + j: " << i + j << std::endl;
				std::cout << "\t    j: " << j << std::endl;
				std::cout << "\tself.rc_map[s[j]]:" << shorter[j] << std::endl;
				std::cout << "\tl[i + j]         :" << longer[i + j] << std::endl;
				std::cout << std::endl;
//					std::cout << "(i + j + 2): " << (i + j ) << std::endl;
//					std::cout << "(len(longer) ): " << (len(longer)) << std::endl;
//					std::cout << "(i + j + 2) ==  (len(longer) ): " << njh::colorBool((i + j + 2) ==  (len(longer) )) << std::endl;
			}
			lastJ = j;
			if ((i + j + 2) == (len(longer))) {
				break;
			}
		}
		//add last matching
		matches.emplace_back(longer[i + lastJ + 1] == shorter[lastJ + 1]);
		// calculate liner extension bonus, score gets an additional -0.5 for every matching bases at the ends if there is overlap

		//front
		uint32_t matchingLeft = 0;
		if (i > 0) {
			for (const auto &matchEnum: iter::enumerate(matches)) {
				if (matchEnum.index >= endBonusLen_) {
					break;
				}
				if (matchEnum.element) {
					++matchingLeft;
				} else {
					break;
				}
			}
		}
		//back
		uint32_t matchingRight = 0;
		if (len(matches) < len(shorter)) {
			uint32_t index = 0;
			for (auto iter = matches.rbegin(); iter < matches.rend(); ++iter) {

				if (index >= endBonusLen_) {
					break;
				}
				if (*iter) {
					++matchingRight;
				} else {
					break;
				}
				++index;
			}
		}
		double linear_extension_bonus = (matchingLeft + matchingRight) * endBonus_;
		if (print) {

			std::cout << "matches: " << njh::conToStr(matches) << std::endl;
			std::cout << "currentScore: " << currentScore << std::endl;
			std::cout << "len(matches): " << len(matches) << std::endl;
			std::cout << "len(shorter): " << len(shorter) << std::endl;
			std::cout << "i > 0: " << njh::colorBool(i > 0) << std::endl;
			std::cout << "len(matches) < len(shorter): " << njh::colorBool(len(matches) < len(shorter)) << std::endl;
			std::cout << "\tmatchingRight: " << matchingRight << std::endl;
			std::cout << "\tmatchingLeft: " << matchingLeft << std::endl;
			std::cout << "linear_extension_bonus: " << linear_extension_bonus << std::endl;
		}
		currentScore += linear_extension_bonus;
		if (print) {
			std::cout << "currentScore: " << currentScore << std::endl;
		}
		if (currentScore < best_score) {
			best_score = currentScore;
			best_matches = matches;
			best_start = i;
		}
	}
	if (print) {
		std::cout << "shorter: " << shorter << ", longer: " << longer << std::endl;
		std::cout << "best_score: " << best_score << std::endl;
		std::cout << "best_matches: " << njh::conToStr(best_matches, ",") << std::endl;
		std::cout << "best_start: " << best_start << std::endl;

		std::cout << std::endl;

		std::cout << seqInfo(longer, longer).getSeqAnsi() << std::endl;
		std::cout << std::string(best_start, ' ');
		for (const auto &m: best_matches) {
			if (m) {
				std::cout << '|';
			} else {
				std::cout << ' ';
			}
		}
		std::cout << std::endl;

		std::cout << std::string(best_start, ' ');
		std::cout << seqInfo(shorter, shorter).getSeqAnsi() << std::endl;
	}
	return best_score;
}

double PrimerDimerUtils::computeDimerScoreTop(const std::string &str1, const std::string &str2) const{
	if (str1.size() >= str2.size()) {
		return (computeDimerScore(seqUtil::reverseComplement(str2, "DNA"), str1));
	} else {
		return (computeDimerScore(seqUtil::reverseComplement(str1, "DNA"), str2));
	}
}


std::unordered_map<std::string, double>PrimerDimerUtils::matchScores = std::unordered_map<std::string, double>{
				{"AT/AT", -0.88},
				{"TA/TA", -0.6},
				{"AA/AA", -1.02},
				{"TT/TT", -1.02},
				{"AC/AC", -1.46},
				{"GT/GT", -1.46},
				{"CA/CA", -1.46},
				{"TG/TG", -1.46},
				{"TC/TC", -1.32},
				{"GA/GA", -1.32},
				{"AG/AG", -1.29},
				{"CT/CT", -1.29},
				{"CG/CG", -2.17},
				{"GC/GC", -2.24},
				{"GG/GG", -1.83},
				{"CC/CC", -1.83}
};

std::unordered_map<std::string, double>PrimerDimerUtils::singleMismatchScores = std::unordered_map<std::string, double>{
				{"AA/AT", 0.61},
				{"AA/AG", 0.88},
				{"AA/AC", 0.14},
				{"AC/AT", 0.77},
				{"AC/AG", 1.33},
				{"AC/AA", 0.64},
				{"AG/AT", 0.02},
				{"AG/AC", -0.13},
				{"AG/AA", 0.71},
				{"AT/AG", 0.73},
				{"AT/AC", 0.07},
				{"AT/AA", 0.69},
				{"CA/CT", 0.43},
				{"CA/CG", 0.75},
				{"CA/CC", 0.03},
				{"CC/CT", 0.79},
				{"CC/CG", 0.7},
				{"CC/CA", 0.62},
				{"CG/CT", 0.11},
				{"CG/CC", -0.11},
				{"CG/CA", -0.47},
				{"CT/CG", 0.4},
				{"CT/CC", -0.32},
				{"CT/CA", -0.12},
				{"GA/GT", 0.17},
				{"GA/GG", 0.81},
				{"GA/GC", -0.25},
				{"GC/GT", 0.47},
				{"GC/GG", 0.79},
				{"GC/GA", 0.62},
				{"GG/GT", -0.52},
				{"GG/GC", -1.11},
				{"GG/GA", 0.08},
				{"GT/GG", 0.98},
				{"GT/GC", -0.59},
				{"GT/GA", 0.45},
				{"TA/TT", 0.69},
				{"TA/TG", 0.92},
				{"TA/TC", 0.42},
				{"TC/TT", 1.33},
				{"TC/TG", 1.05},
				{"TC/TA", 0.97},
				{"TG/TT", 0.74},
				{"TG/TC", 0.44},
				{"TG/TA", 0.43},
				{"TT/TG", 0.75},
				{"TT/TC", 0.34},
				{"TT/TA", 0.68}};

std::vector<std::vector<double>> PrimerDimerUtils::computeFullScoreMatrix(const std::vector<seqInfo> & primers, uint32_t numThreads) const{
	std::vector<std::vector<double>> scores(primers.size(), std::vector<double>(primers.size(), 0.0));
	PairwisePairFactory pFac(primers.size());

	std::function<void()> computeDim = [&pFac, &scores, this, &primers]() {
		PairwisePairFactory::PairwisePairVec pairs;
		while (pFac.setNextPairs(pairs, 300)) {
			for(const auto & pair : pairs.pairs_){
				scores[pair.row_][pair.col_] = computeDimerScoreTop(primers[pair.row_].seq_, primers[pair.col_].seq_);
				scores[pair.col_][pair.row_] = scores[pair.row_][pair.col_];
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(computeDim, numThreads);
	//pair factory doesn't produce the self comparison so take a quick run through the diagonal
	for (const auto i: iter::range(primers.size())) {
		scores[i][i] = computeDimerScoreTop(primers[i].seq_, primers[i].seq_);
	}
//	for (const auto i: iter::range(primers.size())) {
//		std::cout << "\ri:" << i << ", " << static_cast<double>(i) / primers.size() * 100.0;
//		std::cout.flush();
//		for (const auto j: iter::range(i, primers.size())) {
//			scores[i][j] = computeDimerScoreTop(primers[i].seq_, primers[j].seq_);
//			scores[j][i] = scores[i][j];
//		}
//	}
//	std::cout << std::endl;
	return scores;
}

void PrimerDimerUtils::writeMatrix(const std::vector<std::vector<double>> & scores, std::ostream & out, const std::vector<seqInfo> & primers, bool outputRawMatrix ){
	if(outputRawMatrix){
		for(const auto & row : scores){
			out << njh::conToStr(row, ",") << std::endl;
		}
	}else{
		out << "primer_name";
		for(const auto & primer : primers){
			out << "\t" << primer.name_;
		}
		out << std::endl;
		for(const auto & rowEnum : iter::enumerate(scores)){
			out << primers[rowEnum.index].name_ << "\t" << njh::conToStr(rowEnum.element, "\t") << std::endl;
		}
	}
}


}  // namespace njhseq
