#pragma once
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
/*
 * substituteMatrix.hpp
 *
 *  Created on: Apr 11, 2014
 *      Author: nickhathaway
 */

#include "njhseq/utils.h"


namespace njhseq {


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

	void setWithMap(const std::map<char, std::map<char, int32_t>> & mapScores);
	void setWithUnoMap(const std::unordered_map<char, std::unordered_map<char, int32_t>> & mapScores);
	void setWtihDNAArray(const int matchMatrix[4][4]);
	void setWtihBlosum62();
	void setWithPam250();
	void setWithSimple(int32_t match, int32_t mismatch);
	void setWithSimpleCaseInsense(int32_t match, int32_t mismatch);

	void setWithDegenScoringCaseSense(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringCaseInsen(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringLessNCaseInsen(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringLessNCaseSense(int32_t matchScore, int32_t mismatchScore);

	void setWithDegenScoringNoNInRefCaseSense(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringNoNInRefCaseInsen(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringNoNInRefCaseSenseLessN(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringNoNInRefCaseInsenLessN(int32_t matchScore, int32_t mismatchScore);

	void setWithDegenScoringNoNInSeqCaseSense(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringNoNInSeqCaseInsen(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringNoNInSeqCaseSenseLessN(int32_t matchScore, int32_t mismatchScore);
	void setWithDegenScoringNoNInSeqCaseInsenLessN(int32_t matchScore, int32_t mismatchScore);

	void setWithFilename(const std::string & filename);

	std::vector<char> determineLetters()const;
	//print
	void printScores(const std::vector<char> & alphbet, std::ostream & out)const ;
	void printScores(std::ostream & out)const;

	std::array<std::array<double, 127>, 127> getRates(const std::vector<char> & alphbet)const;

  /**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const ;

	//static creating matrixes
	static substituteMatrix createDegenScoreMatrix(int32_t matchScore, int32_t mismatchScore);
	static substituteMatrix createDegenScoreMatrixCaseInsensitive(int32_t matchScore, int32_t mismatchScore);
	static substituteMatrix createDegenScoreMatrixLessN(int32_t matchScore, int32_t mismatchScore);
	static substituteMatrix createDegenScoreMatrixNoNInRef(int32_t matchScore, int32_t mismatchScore);
	static substituteMatrix createScoreMatrix(int32_t matchScore, int32_t mismatchScore,
			bool degenerativeScoring, bool degenLessN, bool caseInsensitive);
};



} /* namespace njhseq */



