#pragma once
/*
 * substituteMatrix.hpp
 *
 *  Created on: Apr 11, 2014
 *      Author: nickhathaway
 */

#include "bibseq/utils.h"


namespace bibseq {


class substituteMatrix {
public:
	//constructor
	/**
	 * Simple 1 for match -1 for mismatch
	 */
	substituteMatrix(){
		setWithSimple(1,-1);
	}
	//constructor
	/**
	 * Given match and mismatch values
	 */
	substituteMatrix(int32_t match, int32_t mismatch){
		setWithSimple(match, mismatch);
	}
	//constructor
	/**
	 * Given a map of maps of score for char to char
	 */
	substituteMatrix(const std::map<char, std::map<char, int32_t>> & mapScores){
		setWithMap(mapScores);
	}
	//constructor
	/**
	 * Given a unordered_map of unordered_maps of score for char to char
	 */
	substituteMatrix(const std::unordered_map<char, std::unordered_map<char, int32_t>> & mapScores){
		setWithUnoMap(mapScores);
	}
	//constructor
	/**
	 * Given the file name of scoring matrix file
	 */
	substituteMatrix(const std::string & filename){
		setWithFilename(filename);
	}
	//constructor
	/**
	 * Given an array with scores for a basic DNA sub matrix
	 */
	substituteMatrix(const int matchMatrix[4][4]){
		setWtihDNAArray(matchMatrix);
	}
	//constructor
	/**
	 * Given an std::array of std::array just put it in
	 */
	substituteMatrix(const std::array<std::array<int32_t, 127>, 127> & mat):mat_(mat){

	}
	//members
	std::array<std::array<int32_t, 127>, 127> mat_;

	//functions

	void setWithZeros();
	void setWithSimple(int32_t match, int32_t mismatch);
	void setWithCaseInsensitive(int32_t match, int32_t mismatch);
	void setWithMap(const std::map<char, std::map<char, int32_t>> & mapScores);
	void setWithUnoMap(const std::unordered_map<char, std::unordered_map<char, int32_t>> & mapScores);
	void setWtihDNAArray(const int matchMatrix[4][4]);
	void setWtihBlosum62();
	void setWithPam250();
	void setWithDegenScoring(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringLessN();
	void setWithDegenScoringCaseInsen(int32_t matchScore, int32_t mismatchScore);
	void setWithFilename(const std::string & filename);

	std::vector<char> determineLetters()const;
	//print
	void printScores(const std::vector<char> & alphbet, std::ostream & out)const ;
	void printScores(std::ostream & out)const;


	//static creating matrixes
	static substituteMatrix createDegenScoreMatrix(int32_t matchScore, int32_t mismatchScore);
	static substituteMatrix createDegenScoreMatrixCaseInsensitive(int32_t matchScore, int32_t mismatchScore);
	static substituteMatrix createDegenScoreMatrixLessN();
	static substituteMatrix createScoreMatrix(int32_t matchScore, int32_t mismatchScore,
			bool degenerativeScoring, bool caseInsensitive);
};



} /* namespace bib */


#ifndef NOT_HEADER_ONLY
#include "substituteMatrix.cpp"
#endif
