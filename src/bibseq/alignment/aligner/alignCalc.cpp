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
 * alnCalc.cpp
 *
 *  Created on: Nov 29, 2015
 *      Author: nick
 */


#include "alignCalc.hpp"


namespace njhseq {

void alignCalc::runSmithSave(const std::string& objA, const std::string& objB,
                         alnParts& parts) {
  // std::cout << "doing smith reg" << std::endl;
  // std::cout<<"mark smith non simple"<<std::endl;
  /*if (currentSetUp_ != "smith") {
    initializeSmith();
  }*/
  // Create the alignment score matrix to do the alignment, a column for each
  // letter in sequence b and a row for each letter in sequence a
  parts.ScoreMatrix_[0][0].leftInherit = 0;
  parts.ScoreMatrix_[0][0].upInherit = 0;
  parts.ScoreMatrix_[0][0].diagInherit = 0;
  parts.ScoreMatrix_[0][0].upInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].leftInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].diagInheritPtr = '\0';
  // initialize first column:
  for (uint32_t i = 1; i < parts.maxSize_; ++i) {
    parts.ScoreMatrix_[i][0].upInherit = 0;
    parts.ScoreMatrix_[i][0].leftInherit = 0;
    parts.ScoreMatrix_[i][0].diagInherit = 0;
    parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
    parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
    parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
  }
  // initialize first row:
  for (uint32_t j = 1; j < parts.maxSize_; ++j) {
    parts.ScoreMatrix_[0][j].upInherit = 0;
    parts.ScoreMatrix_[0][j].leftInherit = 0;
    parts.ScoreMatrix_[0][j].diagInherit = 0;
    parts.ScoreMatrix_[0][j].upInheritPtr = '\0';
    parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';
    parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
  }
  // to find the best score
  // objectA=objA;
  // objectB=objB;
  parts.lHolder_.addFromFile_ = false;
  parts.lHolder_.gapInfos_.clear();
  uint32_t bestJ = 0;
  uint32_t bestI = 0;
  int32_t bestValue = 0;
  // empty the alignment strings and qualities vectors to reset for the new
  // alignment
  // alignObjectA_.clear();
  // alignObjectB_.clear();
  // get the lenth of the strings to create the alignment score matrix
  uint32_t lena = objA.size() + 1;
  uint32_t lenb = objB.size() + 1;
  for (uint32_t i = 1; i < lena; ++i) {
    for (uint32_t j = 1; j < lenb; ++j) {
      char ptrFlag;
      // first set the upInherit score. If we are in the first row, this will
      // be
      // the value of the above cell's leftInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max of the three scores in the cell above,
      // with the appropriate penalty applied (either a
      // parts.gapScores_.gapOpen or
      // parts.gapScores_.gapExtend_).
      parts.ScoreMatrix_[i][j].upInherit =
          smithMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                           parts.gapScores_.gapExtend_,
                       parts.ScoreMatrix_[i - 1][j].leftInherit -
                           parts.gapScores_.gapOpen_,
                       parts.ScoreMatrix_[i - 1][j].diagInherit -
                           parts.gapScores_.gapOpen_,
                       ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

      // next set the leftInherit score. If we are in the first column, this
      // will be the value of the left cell's upInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max score of the three scores in the cell
      // to
      // the left, with the appropriate penalty applied.
      parts.ScoreMatrix_[i][j].leftInherit = smithMaximum(
          parts.ScoreMatrix_[i][j - 1].upInherit - parts.gapScores_.gapOpen_,
          parts.ScoreMatrix_[i][j - 1].leftInherit -
              parts.gapScores_.gapExtend_,
          parts.ScoreMatrix_[i][j - 1].diagInherit -
              parts.gapScores_.gapOpen_,
          ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      // int match = scoringArray[objA.seqBase_.seq_[i -
      // 1]-'A'][objB.seqBase_.seq_[j - 1]-'A'];
      int match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match + smithMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                               parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                               parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                               ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
      int32_t tempValue =
          smithMaximum(parts.ScoreMatrix_[i][j].upInherit,
                       parts.ScoreMatrix_[i][j].leftInherit,
                       parts.ScoreMatrix_[i][j].diagInherit, ptrFlag);
      if (tempValue > bestValue) {
        bestValue = tempValue;
        bestI = i;
        bestJ = j;
      }
    }
  }
  // set the i (row) cursor and j (column) cursor to the bottom right corner

  int icursor = bestI;
  int jcursor = bestJ;
  // std::cout<<"icursor is "<<icursor<<std::endl;
  // std::cout<<"jcursor is "<<jcursor<<std::endl;

  // tracerNext holds to where to go next in the matrix, will be (D) diagonal,
  // (U) up, or (L) left depending on the maximum score determined during the
  // matrix set up.
  char tracerNext = ' ';

  // get the alignment score from the  bottom right cell and set the tacer to
  // where to go next
  // keep tracing back until at the begining of either sequence
  // Traceback algorithm follows. Score is the max of all three scores stored
  // in
  // the bottom right cell.
  // Alignments are constructed by following the correct pointer backwards at
  // each stage.
  // Since alignment strings are constructed in reverse, we must call the
  // reverse() funcion after traceback.
  parts.score_ = smithMaximum(
      parts.ScoreMatrix_[icursor][jcursor].upInherit,
      parts.ScoreMatrix_[icursor][jcursor].leftInherit,
      parts.ScoreMatrix_[icursor][jcursor].diagInherit, tracerNext);
  double currentValue = parts.score_;
  uint32_t gapBSize = 0;
  uint32_t gapASize = 0;
  // while (icursor!=0 || jcursor!=0 || (int)currentValue!=0) {
  while ((icursor != 0 || jcursor != 0) && currentValue != 0) {
    if (tracerNext == 'U') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      currentValue = parts.ScoreMatrix_[icursor][jcursor].upInherit;
      if (currentValue != 0) {
        if (tracerNext != 'U' && tracerNext != 'B') {
          parts.lHolder_.gapInfos_.emplace_back(
              gapInfo(jcursor, gapBSize, false));
          gapBSize = 0;
        }
        icursor--;
      }
    } else if (tracerNext == 'L') {
      ++gapASize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
      currentValue = parts.ScoreMatrix_[icursor][jcursor].leftInherit;
      if (currentValue != 0) {
        if (tracerNext != 'L') {
          parts.lHolder_.gapInfos_.emplace_back(
              gapInfo(icursor, gapASize, true));
          gapASize = 0;
        }
        jcursor--;
      }
    } else if (tracerNext == 'D') {
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].diagInheritPtr;
      currentValue = parts.ScoreMatrix_[icursor][jcursor].diagInherit;
      icursor--;
      jcursor--;
      // overlaps++;
    }
    // if ambigous traceback (can go either left or up), we give precedence
    // to
    // an 'up' traceback. This will not affect the score, of course.
    else if (tracerNext == 'B') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      currentValue = parts.ScoreMatrix_[icursor][jcursor].upInherit;
      if (currentValue != 0) {
        if (tracerNext != 'U' && tracerNext != 'B') {
          parts.lHolder_.gapInfos_.emplace_back(
              gapInfo(jcursor, gapBSize, false));
          gapBSize = 0;
        }
        icursor--;
      }
    } else {
      std::cerr << "ERROR!!!!!" << std::endl;
    }
  }
  parts.lHolder_.localAStart_ = icursor;
  parts.lHolder_.localASize_ = bestI - icursor;
  parts.lHolder_.localBStart_ = jcursor;
  parts.lHolder_.localBSize_ = bestJ - jcursor;
  parts.lHolder_.score_ = parts.score_;
  // rearrangeLocal(objA, objB);
  // alnInfoLocal_.addFromFile_ = false;
  // alnHolders_.locaparts.lHolder__[parts.gapScores_.getIdentifer()]
  //  .addAlnInfo(objA.seq_, objB.seq_, alnInfoLocal_);
}




void alignCalc::runNeedleSave(const std::string& objA, const std::string& objB,
                          alnParts& parts) {
  parts.gHolder_.gapInfos_.clear();
  parts.gHolder_.addFromFile_ = false;
  // Create the alignment score matrix to do the alignment, a column for each
  // letter in sequence b and a row for each letter in sequence a
  parts.ScoreMatrix_[0][0].leftInherit = 0;
  parts.ScoreMatrix_[0][0].upInherit = 0;
  parts.ScoreMatrix_[0][0].diagInherit = 0;
  parts.ScoreMatrix_[0][0].upInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].leftInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].diagInheritPtr = '\0';
  // get the length of the strings to create the alignment score matrix
  uint32_t lena = objA.size() + 1;
  uint32_t lenb = objB.size() + 1;

  // initialize first column:
  {
  		const uint32_t i = 1;
		parts.ScoreMatrix_[i][0].upInherit = - parts.gapScores_.gapLeftQueryOpen_;
		parts.ScoreMatrix_[i][0].leftInherit = 0;
		parts.ScoreMatrix_[i][0].diagInherit = 0;
		parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
		parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
		parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
  }
  for (uint32_t i = 2; i < lena; ++i) {
    parts.ScoreMatrix_[i][0].upInherit =
        parts.ScoreMatrix_[i - 1][0].upInherit -
        parts.gapScores_.gapLeftQueryExtend_;
    parts.ScoreMatrix_[i][0].leftInherit = 0;
    parts.ScoreMatrix_[i][0].diagInherit = 0;
    parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
    parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
    parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
  }
  // initialize first row:
  {
  		const uint32_t j = 1;
		parts.ScoreMatrix_[0][j].upInherit = 0;
		parts.ScoreMatrix_[0][j].leftInherit = - parts.gapScores_.gapLeftRefOpen_;
		parts.ScoreMatrix_[0][j].diagInherit = 0;
		parts.ScoreMatrix_[0][j].upInheritPtr = '\0';
		parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';
		parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
  }
  for (uint32_t j = 2; j < lenb; ++j) {
    parts.ScoreMatrix_[0][j].upInherit = 0;
    parts.ScoreMatrix_[0][j].leftInherit =
        parts.ScoreMatrix_[0][j - 1].leftInherit -
        parts.gapScores_.gapLeftRefExtend_;
    parts.ScoreMatrix_[0][j].diagInherit = 0;
    parts.ScoreMatrix_[0][j].upInheritPtr = '\0';
    parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';
    parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
  }

  //i = 1, j = 1
  {
  	const uint32_t i = 1;
  	const uint32_t j = 1;
		//up inherit is always left since it will have to inherit from the left
		parts.ScoreMatrix_[i][j].upInherit =
				parts.ScoreMatrix_[i - 1][j].leftInherit - parts.gapScores_.gapOpen_;
		parts.ScoreMatrix_[i][j].upInheritPtr = 'L';
  	//left inherit is always coming from up
    parts.ScoreMatrix_[i][j].leftInherit =
        parts.ScoreMatrix_[i][j - 1].upInherit -
        parts.gapScores_.gapOpen_;
    parts.ScoreMatrix_[i][j].leftInheritPtr = 'U';
    //diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit = parts.ScoreMatrix_[i - 1][j - 1].leftInherit + match;
    parts.ScoreMatrix_[i][j].diagInheritPtr = 'L';
  }
  {
    //i = 1
    const uint32_t i = 1;
		for (uint32_t j = 2; j < lenb - 1; ++j) {
			char ptrFlag;
			//up inherit is always left since it will have to inherit from the left
			parts.ScoreMatrix_[i][j].upInherit =
					parts.ScoreMatrix_[i - 1][j].leftInherit - parts.gapScores_.gapOpen_;
			parts.ScoreMatrix_[i][j].upInheritPtr = 'L';

			//regular left inherit
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;

      //diag inherit will also have to be from the left
			int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          parts.ScoreMatrix_[i - 1][j - 1].leftInherit + match;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'L';
		}
  }

  //j = 1
  {
		const uint32_t j = 1;
    for (uint32_t i = 2; i < lena - 1 ; ++i) {
    	char ptrFlag;
    	//regular up inherit
      parts.ScoreMatrix_[i][1].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//left inherit is always coming from up
      parts.ScoreMatrix_[i][j].leftInherit =
          parts.ScoreMatrix_[i][j - 1].upInherit -
          parts.gapScores_.gapOpen_;
      parts.ScoreMatrix_[i][j].leftInheritPtr = 'U';
      //diag inherit is always coming from up
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit = parts.ScoreMatrix_[i - 1][j - 1].upInherit + match;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'U';
    }
  }
  //i = 1, j = lenb - 1
  {
  	char ptrFlag;
  	//
  	const uint32_t j = lenb - 1;
  	const uint32_t i = 1;
    parts.ScoreMatrix_[i][j ].upInherit =
                    parts.ScoreMatrix_[i - 1][j ].leftInherit -
                    parts.gapScores_.gapRightQueryOpen_;
    parts.ScoreMatrix_[i][j ].upInheritPtr = 'L';
		//regular left inherit
    parts.ScoreMatrix_[i][j].leftInherit =
        needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                          parts.gapScores_.gapOpen_,
                      parts.ScoreMatrix_[i][j - 1].leftInherit -
                          parts.gapScores_.gapExtend_,
                      parts.ScoreMatrix_[i][j - 1].diagInherit -
                          parts.gapScores_.gapOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
    //diag inherit will also have to be from the left
		int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        parts.ScoreMatrix_[i - 1][j - 1].leftInherit + match;
    parts.ScoreMatrix_[i][j].diagInheritPtr = 'L';
  }
  //j = 1, i = lena - 1
  {
  	char ptrFlag;
  	const uint32_t i = lena - 1;
  	const uint32_t j = 1;
  	//regular up inherit
    parts.ScoreMatrix_[i][1].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
  	//left
    parts.ScoreMatrix_[i][j].leftInherit =
        parts.ScoreMatrix_[i][j - 1].upInherit -
        parts.gapScores_.gapRightRefOpen_;
    parts.ScoreMatrix_[i][j].leftInheritPtr = 'U';
    //diag inherit is always coming from up
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit = parts.ScoreMatrix_[i - 1][j - 1].upInherit + match;
    parts.ScoreMatrix_[i][j].diagInheritPtr = 'U';
  }
  for (uint32_t i = 2; i < lena - 1 ; ++i) {
    for (uint32_t j = 2; j < lenb - 1 ; ++j) {
    	char ptrFlag;
      // std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
      //char ptrFlag;
      // first set the upInherit score. If we are in the first row, this will
      // be
      // the value of the above cell's leftInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max of the three scores in the cell above,
      // with the appropriate penalty applied (either a
      // parts.gapScores_.gapOpen or
      // parts.gapScores_.gapExtend_).
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

      // next set the leftInherit score. If we are in the first column, this
      // will be the value of the left cell's upInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max score of the three scores in the cell
      // to
      // the left, with the appropriate penalty applied.
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      //diag, match or mismatch
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //i = lena - 1
  {
  	const uint32_t i = lena - 1;
    for (uint32_t j = 2; j < lenb - 1; ++j) {
    	char ptrFlag;
    	//normal up
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//end left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapRightRefOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapRightRefExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapRightRefOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
    	//normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //j = lenb - 1
  {
  	const uint32_t j = lenb - 1;
    for (uint32_t i = 2; i < lena - 1; ++i) {
    	char ptrFlag;
    	//end up
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapRightQueryExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//regular left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
    	//normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //i = lena - 1, j = lenb -1
  {
  	char ptrFlag;
  	const uint32_t i = lena - 1;
  	const uint32_t j = lenb - 1;
  	//end up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapRightQueryExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
  	//end left
    parts.ScoreMatrix_[i][j].leftInherit =
        needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      parts.ScoreMatrix_[i][j - 1].leftInherit -
                          parts.gapScores_.gapRightRefExtend_,
                      parts.ScoreMatrix_[i][j - 1].diagInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  	//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  }

  int icursor = lena - 1;
  int jcursor = lenb - 1;

  // tracerNext holds to where to go next in the matrix, will be (D) diagonal,
  // (U) up, or (L) left depending on the maximum score determined during the
  // matrix set up.
  char tracerNext = ' ';

  // get the alignment score from the  bottom right cell and set the tacer to
  // where to go next
  // keep tracing back until at the begining of either sequence
  // Traceback algorithm follows. Score is the max of all three scores stored
  // in
  // the bottom right cell.
  // Alignments are constructed by following the correct pointer backwards at
  // each stage.
  parts.score_ = needleMaximum(
      parts.ScoreMatrix_[icursor][jcursor].upInherit,
      parts.ScoreMatrix_[icursor][jcursor].leftInherit,
      parts.ScoreMatrix_[icursor][jcursor].diagInherit, tracerNext);
  parts.gHolder_.score_ = parts.score_;
  uint32_t gapBSize = 0;
  uint32_t gapASize = 0;
  // std::cout <<"rnv2" << std::endl;
  while (icursor != 0 || jcursor != 0) {
    if (tracerNext == 'U') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      if (tracerNext != 'U' && tracerNext != 'B') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    } else if (tracerNext == 'L') {
      ++gapASize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
      if (tracerNext != 'L') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(icursor, gapASize, true));
        gapASize = 0;
      }
      --jcursor;
    } else if (tracerNext == 'D') {
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].diagInheritPtr;
      --icursor;
      --jcursor;
    }
    // if ambigous traceback (can go either left or up), we give precedence
    // to an 'up' traceback.
    else if (tracerNext == 'B') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      if (tracerNext != 'U' && tracerNext != 'B') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    }
  }
  if ((tracerNext == 'U' || tracerNext == 'B') && gapBSize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(njhseq::gapInfo(jcursor, gapBSize, false));
  } else if (tracerNext == 'L' && gapASize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(njhseq::gapInfo(icursor, gapASize, true));
  }
}



void alignCalc::runNeedleOnlyEndGapsSave(const std::string& objA, const std::string& objB,
                          alnParts& parts) {
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;

  parts.gHolder_.gapInfos_.clear();
  parts.gHolder_.addFromFile_ = false;
  // Create the alignment score matrix to do the alignment, a column for each
  // letter in sequence b and a row for each letter in sequence a
  parts.ScoreMatrix_[0][0].leftInherit = 0;
  parts.ScoreMatrix_[0][0].upInherit = 0;
  parts.ScoreMatrix_[0][0].diagInherit = 0;
  parts.ScoreMatrix_[0][0].upInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].leftInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].diagInheritPtr = '\0';
  // initialize first column:
  for (uint32_t i = 1; i < parts.maxSize_; ++i) {
    if (i == 1) {
      parts.ScoreMatrix_[i][0].upInherit = - parts.gapScores_.gapLeftQueryOpen_;
      parts.ScoreMatrix_[i][0].leftInherit = 0;
      parts.ScoreMatrix_[i][0].diagInherit = 0;
      parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
      parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
      parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
    } else {
      parts.ScoreMatrix_[i][0].upInherit =
          parts.ScoreMatrix_[i - 1][0].upInherit -
          parts.gapScores_.gapLeftQueryExtend_;
      parts.ScoreMatrix_[i][0].leftInherit = 0;
      parts.ScoreMatrix_[i][0].diagInherit = 0;
      parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
      parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
      parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
    }
  }
  // initialize first row:
  for (uint32_t j = 1; j < parts.maxSize_; ++j) {
    if (j == 1) {
      parts.ScoreMatrix_[0][j].upInherit = 0;
      parts.ScoreMatrix_[0][j].leftInherit = - parts.gapScores_.gapLeftRefOpen_;
      parts.ScoreMatrix_[0][j].diagInherit = 0;
      parts.ScoreMatrix_[0][j].upInheritPtr = '\0';
      parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';
      parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
    } else {
      parts.ScoreMatrix_[0][j].upInherit = 0;
      parts.ScoreMatrix_[0][j].leftInherit =
          parts.ScoreMatrix_[0][j - 1].leftInherit -
          parts.gapScores_.gapLeftRefExtend_;
      parts.ScoreMatrix_[0][j].diagInherit = 0;
      parts.ScoreMatrix_[0][j].upInheritPtr = '\0';
      parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';
      parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
    }
  }
  // get the length of the strings to create the alignment score matrix
  uint32_t lena = objA.size() + 1;
  uint32_t lenb = objB.size() + 1;
  //i = 1, j = 1
  {
  	const uint32_t i = 1;
  	const uint32_t j = 1;
    //diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit = parts.ScoreMatrix_[i - 1][j - 1].leftInherit + match;
    parts.ScoreMatrix_[i][j].diagInheritPtr = 'L';
  }
  {
    //i = 1
    const uint32_t i = 1;
		for (uint32_t j = 2; j < lenb - 1; ++j) {
      //diag inherit will also have to be from the left
			int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          parts.ScoreMatrix_[i - 1][j - 1].leftInherit + match;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'L';
		}
  }

  //j = 1
  {
		const uint32_t j = 1;
    for (uint32_t i = 2; i < lena - 1 ; ++i) {
      //diag inherit is always coming from up
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit = parts.ScoreMatrix_[i - 1][j - 1].upInherit + match;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'U';
    }
  }
  //i = 1, j = lenb - 1
  {
  //	char ptrFlag;
  	//
  	const uint32_t j = lenb - 1;
  	const uint32_t i = 1;
    //diag inherit will also have to be from the left
		int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        parts.ScoreMatrix_[i - 1][j - 1].leftInherit + match;
    parts.ScoreMatrix_[i][j].diagInheritPtr = 'L';
    //up inherit has to be left, only way to go
    parts.ScoreMatrix_[i][j ].upInherit =
                    parts.ScoreMatrix_[i - 1][j ].leftInherit -
                    parts.gapScores_.gapRightQueryOpen_;
    parts.ScoreMatrix_[i][j ].upInheritPtr = 'L';


  }
  //j = 1, i = lena - 1
  {
		//	char ptrFlag;
		const uint32_t i = lena - 1;
		const uint32_t j = 1;
  		//left is always up as that is the only to go from here
    parts.ScoreMatrix_[i][j].leftInherit =
        parts.ScoreMatrix_[i][j - 1].upInherit -
        parts.gapScores_.gapRightRefOpen_;
    parts.ScoreMatrix_[i][j].leftInheritPtr = 'U';
    //diag inherit is always coming from up
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit = parts.ScoreMatrix_[i - 1][j - 1].upInherit + match;
    parts.ScoreMatrix_[i][j].diagInheritPtr = 'U';
  }

  for (uint32_t i = 2; i < lena - 1 ; ++i) {
    for (uint32_t j = 2; j < lenb - 1 ; ++j) {
    	//char ptrFlag;
      //diag, match or mismatch
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =  match + parts.ScoreMatrix_[i - 1][j - 1].diagInherit;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'D';
    }
  }

  //i = lena - 1, last row
  {
  	const uint32_t i = lena - 1;
		for (uint32_t j = 2; j < lenb - 1; ++j) {
			char ptrFlag;
			//end left
			if (parts.ScoreMatrix_[i][j - 1].leftInherit
					- parts.gapScores_.gapRightRefExtend_
					> parts.ScoreMatrix_[i][j - 1].diagInherit
							- parts.gapScores_.gapRightRefOpen_) {
				parts.ScoreMatrix_[i][j].leftInherit =
						parts.ScoreMatrix_[i][j - 1].leftInherit
								- parts.gapScores_.gapRightRefExtend_;
				ptrFlag = 'L';
			} else {
				parts.ScoreMatrix_[i][j].leftInherit =
						parts.ScoreMatrix_[i][j - 1].diagInherit
								- parts.gapScores_.gapRightRefOpen_;
				ptrFlag = 'D';
			}
			parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
    	//normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit = match + parts.ScoreMatrix_[i - 1][j - 1].diagInherit;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'D';
    }
  }

  //j = lenb - 1, last column
  {
  	const uint32_t j = lenb - 1;
    for (uint32_t i = 2; i < lena - 1; ++i) {
    		char ptrFlag;
    		//end up
    		if(parts.ScoreMatrix_[i - 1][j].upInherit -
            parts.gapScores_.gapRightQueryExtend_ >
    			parts.ScoreMatrix_[i - 1][j].diagInherit -
    		      parts.gapScores_.gapRightQueryOpen_){
    			ptrFlag = 'U';
    			parts.ScoreMatrix_[i][j].upInherit = parts.ScoreMatrix_[i - 1][j].upInherit -
              parts.gapScores_.gapRightQueryExtend_;
    		}else{
    			parts.ScoreMatrix_[i][j].upInherit = parts.ScoreMatrix_[i - 1][j].diagInherit -
    		      parts.gapScores_.gapRightQueryOpen_;
    			ptrFlag = 'D';
    		}
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

    	//normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match + parts.ScoreMatrix_[i - 1][j - 1].diagInherit;
      parts.ScoreMatrix_[i][j].diagInheritPtr = 'D';
    }
  }

	//i = lena - 1, j = lenb -1, very last cell
	{
		char ptrFlag;
		const uint32_t i = lena - 1;
		const uint32_t j = lenb - 1;
		//end up
		if (parts.ScoreMatrix_[i - 1][j].upInherit
				- parts.gapScores_.gapRightQueryExtend_
				> parts.ScoreMatrix_[i - 1][j].diagInherit
						- parts.gapScores_.gapRightQueryOpen_) {
			ptrFlag = 'U';
			parts.ScoreMatrix_[i][j].upInherit =
					parts.ScoreMatrix_[i - 1][j].upInherit
							- parts.gapScores_.gapRightQueryExtend_;
		} else {
			parts.ScoreMatrix_[i][j].upInherit =
					parts.ScoreMatrix_[i - 1][j].diagInherit
							- parts.gapScores_.gapRightQueryOpen_;
			ptrFlag = 'D';
		}
		parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

		//end left
		if (parts.ScoreMatrix_[i][j - 1].leftInherit
				- parts.gapScores_.gapRightRefExtend_
				> parts.ScoreMatrix_[i][j - 1].diagInherit
						- parts.gapScores_.gapRightRefOpen_) {
			parts.ScoreMatrix_[i][j].leftInherit =
					parts.ScoreMatrix_[i][j - 1].leftInherit
							- parts.gapScores_.gapRightRefExtend_;
			ptrFlag = 'L';
		} else {
			parts.ScoreMatrix_[i][j].leftInherit =
					parts.ScoreMatrix_[i][j - 1].diagInherit
							- parts.gapScores_.gapRightRefOpen_;
			ptrFlag = 'D';
		}
		parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
		//normal diag
		int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
		parts.ScoreMatrix_[i][j].diagInherit = match
				+ parts.ScoreMatrix_[i - 1][j - 1].diagInherit;
		parts.ScoreMatrix_[i][j].diagInheritPtr = 'D';
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;

  int icursor = lena - 1;
  int jcursor = lenb - 1;

  // tracerNext holds to where to go next in the matrix, will be (D) diagonal,
  // (U) up, or (L) left depending on the maximum score determined during the
  // matrix set up.
  char tracerNext = ' ';

  // get the alignment score from the  bottom right cell and set the tacer to
  // where to go next
  // keep tracing back until at the begining of either sequence
  // Traceback algorithm follows. Score is the max of all three scores stored
  // in
  // the bottom right cell.
  // Alignments are constructed by following the correct pointer backwards at
  // each stage.
  parts.score_ = needleMaximum(
      parts.ScoreMatrix_[icursor][jcursor].upInherit,
      parts.ScoreMatrix_[icursor][jcursor].leftInherit,
      parts.ScoreMatrix_[icursor][jcursor].diagInherit, tracerNext);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;

  parts.gHolder_.score_ = parts.score_;
  uint32_t gapBSize = 0;
  uint32_t gapASize = 0;
  if(tracerNext == 'B'){
  		tracerNext = 'U';
  }
  // std::cout <<"rnv2" << std::endl;
  while (icursor != 0 || jcursor != 0) {
		if (tracerNext == 'U') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      if (tracerNext != 'U') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    } else if (tracerNext == 'L') {
      ++gapASize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
      if (tracerNext != 'L') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(icursor, gapASize, true));
        gapASize = 0;
      }
      --jcursor;
    } else if (tracerNext == 'D') {
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].diagInheritPtr;
      --icursor;
      --jcursor;
    }
  }
  if ((tracerNext == 'U') && gapBSize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(njhseq::gapInfo(jcursor, gapBSize, false));
  } else if (tracerNext == 'L' && gapASize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(njhseq::gapInfo(icursor, gapASize, true));
  }
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;

}


alignCalc::MatCursor alignCalc::runNeedleDiagonalSaveInit(
									const std::string& objA,
									const std::string& objB,
									uint32_t alignmentBlockSize,
                   alnParts& parts) noexcept {

  // Create the alignment score matrix to do the alignment, a column for each
  // letter in sequence b and a row for each letter in sequence a
  parts.ScoreMatrix_[0][0].leftInherit = 0;
  parts.ScoreMatrix_[0][0].upInherit = 0;
  parts.ScoreMatrix_[0][0].diagInherit = 0;
  parts.ScoreMatrix_[0][0].upInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].leftInheritPtr = '\0';
  parts.ScoreMatrix_[0][0].diagInheritPtr = '\0';
  // get the length of the strings to create the alignment score matrix
  const uint32_t lena = std::min<uint32_t>(alignmentBlockSize, objA.size() + 1);
  const uint32_t lenb = std::min<uint32_t>(alignmentBlockSize, objB.size() + 1);

  uint32_t stopa = lena == objA.size() + 1 ? lena - 1 : lena;
  uint32_t stopb = lenb == objB.size() + 1 ? lenb - 1 : lenb;

//  std::cout << __PRETTY_FUNCTION__ << std::endl;
//  std::cout << "alignmentBlockSize: " << alignmentBlockSize << std::endl;
//  std::cout << "objA.size() + 1:    " << objA.size() + 1 << std::endl;
//  std::cout << "lena:               " << lena << std::endl;
//  std::cout << "objB.size() + 1:    " << objB.size() + 1 << std::endl;
//  std::cout << "lenb:               " << lenb << std::endl;
//  std::cout << "std::numeric_limits<int32_t>::lowest(): " << std::numeric_limits<int32_t>::lowest() << std::endl;
//	std::cout << "stopa:              " << stopa << std::endl;
//	std::cout << "stopb:              " << stopb << std::endl;
//	std::cout << "obja: " << objA << std::endl;
//	std::cout << "objb: " << objB << std::endl;
//  std::cout << std::endl;
  // initialize first column:
  {
  		const uint32_t i = 1;
		parts.ScoreMatrix_[i][0].upInherit = - parts.gapScores_.gapLeftQueryOpen_;
		parts.ScoreMatrix_[i][0].upInheritPtr = 'U';

		parts.ScoreMatrix_[i][0].leftInherit = std::numeric_limits<int32_t>::lowest()/2;
		parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';

		parts.ScoreMatrix_[i][0].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
		parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
  }
  for (uint32_t i = 2; i < lena; ++i) {
    parts.ScoreMatrix_[i][0].upInherit =
        parts.ScoreMatrix_[i - 1][0].upInherit -
        parts.gapScores_.gapLeftQueryExtend_;
    parts.ScoreMatrix_[i][0].upInheritPtr = 'U';

    parts.ScoreMatrix_[i][0].leftInherit = std::numeric_limits<int32_t>::lowest()/2;
    parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
    parts.ScoreMatrix_[i][0].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
    parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';

  }
  // initialize first row:
  {
  		const uint32_t j = 1;

		parts.ScoreMatrix_[0][j].leftInherit = - parts.gapScores_.gapLeftRefOpen_;
		parts.ScoreMatrix_[0][j].upInheritPtr = '\0';

		parts.ScoreMatrix_[0][j].upInherit = std::numeric_limits<int32_t>::lowest()/2;
		parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';

		parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
		parts.ScoreMatrix_[0][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;

  }
  for (uint32_t j = 2; j < lenb; ++j) {
    parts.ScoreMatrix_[0][j].leftInherit =
        parts.ScoreMatrix_[0][j - 1].leftInherit -
        parts.gapScores_.gapLeftRefExtend_;
    parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';

    parts.ScoreMatrix_[0][j].upInherit = std::numeric_limits<int32_t>::lowest()/2;
    parts.ScoreMatrix_[0][j].upInheritPtr = '\0';

    parts.ScoreMatrix_[0][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
    parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
  }





  for (uint32_t i = 1; i < stopa ; ++i) {
    for (uint32_t j = 1; j < stopb ; ++j) {
    	char ptrFlag;
      // std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
      //char ptrFlag;
      // first set the upInherit score. If we are in the first row, this will
      // be
      // the value of the above cell's leftInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max of the three scores in the cell above,
      // with the appropriate penalty applied (either a
      // parts.gapScores_.gapOpen or
      // parts.gapScores_.gapExtend_).
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

      // next set the leftInherit score. If we are in the first column, this
      // will be the value of the left cell's upInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max score of the three scores in the cell
      // to
      // the left, with the appropriate penalty applied.
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      //diag, match or mismatch
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }


  // i = lena - 1
  // end gap
  //
  if(stopa == lena - 1){
  	  const uint32_t i = lena - 1;
    for (uint32_t j = 1; j < stopb; ++j) {
     	char ptrFlag;
    	  //normal up
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
      	//end left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapRightRefOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapRightRefExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapRightRefOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      	//normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //j = lenb - 1
  //end gap
  if(stopb == lenb - 1){
   	const uint32_t j = lenb - 1;
    for (uint32_t i = 1; i < stopa; ++i) {
			char ptrFlag;
			//end up
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapRightQueryExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//regular left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      //normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //i = lena - 1, j = lenb -1
  //end
  if(stopa == lena - 1 && stopb == lenb - 1) {
  		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " putting in end gap " << std::endl;
		char ptrFlag;
		const uint32_t i = lena - 1;
		const uint32_t j = lenb - 1;
		//end up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapRightQueryExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    //end left
    parts.ScoreMatrix_[i][j].leftInherit =
        needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      parts.ScoreMatrix_[i][j - 1].leftInherit -
                          parts.gapScores_.gapRightRefExtend_,
                      parts.ScoreMatrix_[i][j - 1].diagInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  		//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  } else if(stopa == lena - 1 && stopb != lenb - 1) {
		char ptrFlag;
		const uint32_t i = stopa;
		const uint32_t j = stopb;
   	//normal up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    //end left
    parts.ScoreMatrix_[i][j].leftInherit =
        needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      parts.ScoreMatrix_[i][j - 1].leftInherit -
                          parts.gapScores_.gapRightRefExtend_,
                      parts.ScoreMatrix_[i][j - 1].diagInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  		//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  } else if(stopa != lena - 1 && stopb == lenb - 1) {
		char ptrFlag;
		const uint32_t i = stopa;
		const uint32_t j = stopb;
		//end up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapRightQueryExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
	  //regular left
		parts.ScoreMatrix_[i][j].leftInherit =
				needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
													parts.gapScores_.gapOpen_,
											parts.ScoreMatrix_[i][j - 1].leftInherit -
													parts.gapScores_.gapExtend_,
											parts.ScoreMatrix_[i][j - 1].diagInherit -
													parts.gapScores_.gapOpen_,
											ptrFlag);
		parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  		//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  }

  int icursor = lena - 1;
  int jcursor = lenb - 1;

  if(lena == objA.size() + 1 && lenb == objB.size() + 1){
  		return MatCursor{icursor, jcursor, lena, lenb};
  } else {
    // tracerNext holds to where to go next in the matrix, will be (D) diagonal,
    // (U) up, or (L) left depending on the maximum score determined during the
    // matrix set up.
    char tracerNext = ' ';
    // get the alignment score from the  bottom right cell and set the tacer to
    // where to go next
    // keep tracing back until at the beginning of either sequence
    // Traceback algorithm follows. Score is the max of all three scores stored
    // in
    // the bottom right cell.
    // Alignments are constructed by following the correct pointer backwards at
    // each stage.
    needleMaximum(
        parts.ScoreMatrix_[icursor][jcursor].upInherit,
        parts.ScoreMatrix_[icursor][jcursor].leftInherit,
        parts.ScoreMatrix_[icursor][jcursor].diagInherit, tracerNext);
    // std::cout <<"rnv2" << std::endl;
    while (icursor != 0 || jcursor != 0) {
//    	std::cout << "within icursor: " << icursor << "/" << lena << " within jcursor: " << jcursor<< "/" << lenb << " tracerNext: " << tracerNext << std::endl;

      if (tracerNext == 'U') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
        if (tracerNext != 'U' && tracerNext != 'B') {
        }
        --icursor;
      } else if (tracerNext == 'L') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
        --jcursor;
      } else if (tracerNext == 'D') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].diagInheritPtr;
        break;
      }
      // if ambigous traceback (can go either left or up), we give precedence
      // to an 'up' traceback.
      else if (tracerNext == 'B') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
        --icursor;
      }
//			if (tracerNext == '\0') {
//				std::cout << "null: within icursor: " << icursor << "/" << lena
//						<< " within jcursor: " << jcursor << "/" << lenb << " tracerNext: "
//						<< tracerNext << std::endl;
//				exit(1);
//			}
    }
    if(0 == icursor && 0 == jcursor){
    		return MatCursor{static_cast<int32_t>(lena - 1), static_cast<int32_t>(lenb - 1), lena, lenb};
    }else{
    		return MatCursor{icursor, jcursor, lena, lenb};
    }
  }
}

alignCalc::MatCursor alignCalc::runNeedleDiagonalSaveStep(
									const std::string& objA,
									const std::string& objB,
									uint32_t alignmentBlockSize,
									uint32_t alignmentBlockWalkbackSize,
									const MatCursor & lastCursors,
                   alnParts& parts) noexcept {
//  bool print = false;
//  if("AGTGAAGAAGCAATGGTGCCTTTGGAAAGTATGTAGCTCTTAGTTAAGCTTTTATCTGGAAGGTGCTGGAGAAGGTCTGGGCACAGAGCTGGTTGGCTCCATGGTTTTGTGACCTGTGCCATGATTTGAGAAGTGCCCATAGGTAGAAATCCAGGATGACTATGGGGCTCGTTTTGTATGTTTCCTGCCTCTCAAGGATTACAGACCTGTCTCTCTGTTGTCCAATGACTCTAAAGAATTGCTTCATAGATTTTGTTCAGTTTCACAGCTTTCTATAAGAGAAGTTAAATATGATACCTGTTACTCCTTCATGGCCAGAACTGGAAGCTAAGCATGGGTTTGATGTCACGCAGACCTGGGTTTGGACCCCAGCTCTGCCCCTTAATGGCTCTGTGACCTAAAGCAGGTTATTTCATCTCTTTGAGCCTCAGTGTCTTTACCTAAGAAACAGCAGGTTATTGTGAGGACCAAATAAAACAAAGCACACTTAGC" == objA &&
//  	 "TCAACATCACTAGGCAATAGGAATTTTTCAGCTCCATTTTAATCTTAAGGGACCACTGGTAATTTTAATTGGGTCTTGGATATTATATATTAAGATTTTGAGATATAATTTTAGGCTATGGGTGATTTTATTTTTTTTCCTGAGAAAGTGAAGTTTGCTTATTTTAGGCAACTGGGGTGAGGACATTAGTAATCCAAGGCCACTTTAAGCCAATCAGGTATTTATGTAATTTGAAGCTGGCCTTGATTCCCTACTATGGCTAGTTTGTTTTTAGTTTACTCTCACTCTTGGGGTATAGCCCTTCAAGTTCCCAGCCAAAAGCCTGGGGTGTTGGCCAAGGTTTGTTCTTCATTGAGGGCCCTGAGTACTAGTTTCTATCTGCCCAGATCCATAAGACTTGGTAAAAGCTCTCTTCAGCTCCTCAGCCTCTCAGCCCGGATGATGCATCAGACACACCTGGGACCCTCGGAGCTCCTGCTGTC" == objB){
//  	print = true;
//  }
//	bool print = false;
//	if ("CCCAGAGCTAGCACAGGCAGAACCCTCATGTGAATCAAAAAGGACACCACTCCGTCCCTGCTAACCCTGTCCTCTGTCTCTTCATTTTCCACCTAGACTGATATTGCCCCCAAACATTAAGCTTTCTTGTATTACCTTCATGAT"
//			== objB
//			&& "TATGTATACCTTCCTAAGATTAGCGGTGTCCTCTTGATTGGGGGATGGTTATAAATAGCTTTTAAGGTGATTTCCAGGCCTAACATTCTATTATTGACTTATAACTGCACTTTAAAGGCTTAAACACACACAATCCCATTGAATCCTCAAGATACTTCCATTGGCAGGCAAGTATTTTCCTTG"
//					== objA) {
//		print = true;
//	}
  const uint32_t starta = std::max(lastCursors.icursor_ + 1 > static_cast<int64_t>(alignmentBlockWalkbackSize) ? lastCursors.icursor_ + 1 - alignmentBlockWalkbackSize : 0, lastCursors.lastStartA_);
  const uint32_t startb = std::max(lastCursors.jcursor_ + 1 > static_cast<int64_t>(alignmentBlockWalkbackSize) ? lastCursors.jcursor_ + 1 - alignmentBlockWalkbackSize : 0, lastCursors.lastStartB_);

  // Create the alignment score matrix to do the alignment, a column for each
  // letter in sequence b and a row for each letter in sequence a
  // get the length of the strings to create the alignment score matrix
  const uint32_t lena = std::min<uint32_t>(starta + alignmentBlockSize, objA.size() + 1);
  const uint32_t lenb = std::min<uint32_t>(startb + alignmentBlockSize, objB.size() + 1);

//  const uint32_t fillStartA = startb < lastCursors.lastStartB_ ? starta : lastCursors.lena_;
//  const uint32_t fillStartB = starta < lastCursors.lastStartA_ ? startb : lastCursors.lenb_;

//  bool aFull = fillStartA == lena;
//  bool bFull = fillStartB == lenb;
  bool aFull = lastCursors.lena_ == lena;
  bool bFull = lastCursors.lenb_ == lenb;
//  bool aFull = false;
//  bool bFull = false;

  uint32_t stopa = (lena == objA.size() + 1) ? lena - 1 : lena;
  uint32_t stopb = (lenb == objB.size() + 1) ? lenb - 1 : lenb;
//  std::cout  << "lena: " << lena << std::endl;
//  std::cout  << "lenb: " << lenb << std::endl;
//  std::cout  << "lastCursors.lena_: " << lastCursors.lena_ << std::endl;
//  std::cout  << "lastCursors.lenb_: " << lastCursors.lenb_ << std::endl;
//  std::cout  << "fillStartA: " << fillStartA << std::endl;
//  std::cout  << "fillStartB: " << fillStartB << std::endl;
//  std::cout  << "aFull: " << njh::colorBool(aFull) << std::endl;
//  std::cout  << "bFull: " << njh::colorBool(bFull) << std::endl;
//  std::cout << std::endl;
//  if(print){
//  		std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//  	  std::cout  << "lena: " << lena << std::endl;
//  	  std::cout  << "lenb: " << lenb << std::endl;
////  	  std::cout  << "fillStartA: " << fillStartA << std::endl;
////  	  std::cout  << "fillStartB: " << fillStartB << std::endl;
//  	  std::cout  << "lastCursors.lena_: " << lastCursors.lena_ << std::endl;
//  	  std::cout  << "lastCursors.lenb_: " << lastCursors.lenb_ << std::endl;
//  	  std::cout  << "aFull: " << njh::colorBool(aFull) << std::endl;
//  	  std::cout  << "bFull: " << njh::colorBool(bFull) << std::endl;
//  	  std::cout << std::endl;
//  }
  // initialize first column
  if(!aFull){
  		const uint32_t j = startb;
    if(0 == startb){
    		//back at 0, end gap
      for (uint32_t i = lastCursors.lena_; i < stopa; ++i) {
        parts.ScoreMatrix_[i][j].upInherit =
            parts.ScoreMatrix_[i - 1][j].upInherit -
            parts.gapScores_.gapLeftQueryExtend_;
        parts.ScoreMatrix_[i][j].upInheritPtr = 'U';

        parts.ScoreMatrix_[i][j].leftInheritPtr = '\0';
        parts.ScoreMatrix_[i][j].leftInherit = std::numeric_limits<int32_t>::lowest()/2;

        parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
        parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
      }
    } else {
    		//not an end gap
  		{
  			//start
  			const uint32_t i = lastCursors.lena_;
//  			if(print){
//  				std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//  				std::cout << "\t" << "i: " << i << std::endl;
//  				std::cout << "\t" << "j: " << j << std::endl;
//  				std::cout << "\t" << "lastCursors.lena_: " << lastCursors.lena_ << std::endl;
//  				std::cout << "\t" << "lastCursors.lenb_: " << lastCursors.lenb_ << std::endl;
//  				std::cout << "\t" << "lena             : " << lena << std::endl;
//  				std::cout << "\t" << "lenb             : " << lenb << std::endl;
//  				std::cout << "\t" << "startb == lastCursors.lastStartB_: " << njh::colorBool(startb == lastCursors.lastStartB_) << std::endl;
//  				std::cout << "\t" << "starta > lastCursors.lastStartA_ : " << njh::colorBool(starta > lastCursors.lastStartA_) << std::endl;
//
//  			}
  			char ptrFlag;
  			//regular up inherit
  			parts.ScoreMatrix_[i][j].upInherit =
  					needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
  														parts.gapScores_.gapExtend_,
  												parts.ScoreMatrix_[i - 1][j].leftInherit -
  														parts.gapScores_.gapOpen_,
  												parts.ScoreMatrix_[i - 1][j].diagInherit -
  														parts.gapScores_.gapOpen_,
  												ptrFlag);
  			parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
  			//if to where we filled to last time for b doesn't change but a does no reg diag
  			if(lenb == lastCursors.lenb_ && lena > lastCursors.lena_){
          parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
          parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
  			}else{
    			//regular diag inherit
    			int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
          parts.ScoreMatrix_[i][j].diagInherit =
              match +
              needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                            parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                            parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                            ptrFlag);
          parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  			}


        //first column so left inherit is impossible
        parts.ScoreMatrix_[i][j].leftInherit = std::numeric_limits<int32_t>::lowest()/2;
        parts.ScoreMatrix_[i][j].leftInheritPtr = '\0';
  		}
      for (uint32_t i = lastCursors.lena_ + 1; i < stopa; ++i) {
        parts.ScoreMatrix_[i][j].upInherit =
            parts.ScoreMatrix_[i - 1][j].upInherit -
            parts.gapScores_.gapExtend_;
        parts.ScoreMatrix_[i][j].upInheritPtr = 'U';

        //first column so left inherit and diag inherit is impossible
        parts.ScoreMatrix_[i][j].leftInherit = std::numeric_limits<int32_t>::lowest()/2;
        parts.ScoreMatrix_[i][j].leftInheritPtr = '\0';
        parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
        parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
      }
    }
  }

  // initialize first row:
  if(!bFull){
  		const uint32_t i = starta;
  		if(0 == starta){
  			//end gap
  		  for (uint32_t j = lastCursors.lenb_; j < stopb; ++j) {
  		    parts.ScoreMatrix_[i][j].leftInherit =
  		        parts.ScoreMatrix_[i][j - 1].leftInherit -
  		        parts.gapScores_.gapLeftRefExtend_;
  		    parts.ScoreMatrix_[i][j].leftInheritPtr = 'L';

  		    //first row so no up or diagonal inherit
  		    parts.ScoreMatrix_[i][j].upInheritPtr = '\0';
  		    parts.ScoreMatrix_[i][j].upInherit = std::numeric_limits<int32_t>::lowest()/2;

  		    parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
  		    parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;

  		  }
  		} else {
  			//not an end gap
  			{
    			//start
    			const uint32_t j = lastCursors.lenb_;
    			char ptrFlag;
//    			if(print){
//    				std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//    				std::cout << "\t" << "i: " << i << std::endl;
//    				std::cout << "\t" << "j: " << j << std::endl;
//
//    			}
    			//if to where we filled to last time for b doesn't change but a does no reg diag
    			if(lena == lastCursors.lena_ && lenb > lastCursors.lenb_){
    		    parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
    		    parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
    			}else{
      			//regular diag inherit
      			int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
  					parts.ScoreMatrix_[i][j].diagInherit =
  							match +
  							needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
  														parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
  														parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
  														ptrFlag);
  					parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    			}


					//first row so no up inherit
					parts.ScoreMatrix_[i][j].upInheritPtr = '\0';
					parts.ScoreMatrix_[i][j].upInherit = std::numeric_limits<int32_t>::lowest()/2;

					//regular left inherit
					parts.ScoreMatrix_[i][j].leftInherit =
							needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
																parts.gapScores_.gapOpen_,
														parts.ScoreMatrix_[i][j - 1].leftInherit -
																parts.gapScores_.gapExtend_,
														parts.ScoreMatrix_[i][j - 1].diagInherit -
																parts.gapScores_.gapOpen_,
														ptrFlag);
					parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  			}
  	    for (uint32_t j = lastCursors.lenb_ + 1; j < stopb; ++j) {

  		    parts.ScoreMatrix_[i][j].leftInherit =
  		        parts.ScoreMatrix_[i][j - 1].leftInherit -
  		        parts.gapScores_.gapExtend_;
  		    parts.ScoreMatrix_[i][j].leftInheritPtr = 'L';

  		    //first row so no up or diagonal inherit
  		    parts.ScoreMatrix_[i][j].upInheritPtr = '\0';
  		    parts.ScoreMatrix_[i][j].upInherit = std::numeric_limits<int32_t>::lowest()/2;

  		    parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
  		    parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
  	    }
  		}
  }

  //fill bottom
  for (uint32_t i = lastCursors.lena_; i < stopa ; ++i) {
    for (uint32_t j = startb + 1; j < lastCursors.lenb_; ++j) {
    	char ptrFlag;
      // std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
      //char ptrFlag;
      // first set the upInherit score. If we are in the first row, this will
      // be
      // the value of the above cell's leftInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max of the three scores in the cell above,
      // with the appropriate penalty applied (either a
      // parts.gapScores_.gapOpen or
      // parts.gapScores_.gapExtend_).
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

      // next set the leftInherit score. If we are in the first column, this
      // will be the value of the left cell's upInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max score of the three scores in the cell
      // to
      // the left, with the appropriate penalty applied.
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      //diag, match or mismatch
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }
  //fill right
  for (uint32_t i = starta + 1; i < lastCursors.lena_ ; ++i) {
    for (uint32_t j = lastCursors.lenb_; j < stopb; ++j) {
    	char ptrFlag;
      // std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
      //char ptrFlag;
      // first set the upInherit score. If we are in the first row, this will
      // be
      // the value of the above cell's leftInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max of the three scores in the cell above,
      // with the appropriate penalty applied (either a
      // parts.gapScores_.gapOpen or
      // parts.gapScores_.gapExtend_).
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

      // next set the leftInherit score. If we are in the first column, this
      // will be the value of the left cell's upInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max score of the three scores in the cell
      // to
      // the left, with the appropriate penalty applied.
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      //diag, match or mismatch
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }
  //fill rest
  for (uint32_t i = lastCursors.lena_; i < stopa ; ++i) {
    for (uint32_t j = lastCursors.lenb_; j < stopb; ++j) {
    	char ptrFlag;
      // std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
      //char ptrFlag;
      // first set the upInherit score. If we are in the first row, this will
      // be
      // the value of the above cell's leftInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max of the three scores in the cell above,
      // with the appropriate penalty applied (either a
      // parts.gapScores_.gapOpen or
      // parts.gapScores_.gapExtend_).
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;

      // next set the leftInherit score. If we are in the first column, this
      // will be the value of the left cell's upInherit score minus a gap open
      // penalty.
      // Otherwise, it will be the max score of the three scores in the cell
      // to
      // the left, with the appropriate penalty applied.
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      //diag, match or mismatch
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }


  //i = lena - 1
	//end gap
	if (stopa == lena - 1) {
		const int32_t i = lena - 1;
		{
			const uint32_t j = startb;
			if (0 == startb) {
				//end gap
				parts.ScoreMatrix_[i][j].upInherit =
						parts.ScoreMatrix_[i - 1][j].upInherit
								- parts.gapScores_.gapLeftQueryExtend_;
				parts.ScoreMatrix_[i][j].upInheritPtr = 'U';
			} else {
				char ptrFlag;
				//normal up
				parts.ScoreMatrix_[i][j].upInherit = needleMaximum(
						parts.ScoreMatrix_[i - 1][j].upInherit
								- parts.gapScores_.gapExtend_,
						parts.ScoreMatrix_[i - 1][j].leftInherit
								- parts.gapScores_.gapOpen_,
						parts.ScoreMatrix_[i - 1][j].diagInherit
								- parts.gapScores_.gapOpen_, ptrFlag);
				parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
			}

			if (i == lastCursors.icursor_ + 1) {
				//just one down from last time, so normal diagonal
				char ptrFlag;
				//normal diag
				int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
				parts.ScoreMatrix_[i][j].diagInherit = match
						+ needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
								parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
								parts.ScoreMatrix_[i - 1][j - 1].diagInherit, ptrFlag);
				parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
			} else {
				//diagonal is impossible
				parts.ScoreMatrix_[i][j].diagInherit =
						std::numeric_limits<int32_t>::lowest() / 2;
				parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
			}

			//first column so left inherit is impossible
			parts.ScoreMatrix_[i][j].leftInherit =
					std::numeric_limits<int32_t>::lowest() / 2;
			parts.ScoreMatrix_[i][j].leftInheritPtr = '\0';
		}

    for (uint32_t j = startb + 1; j < stopb; ++j) {
    	  char ptrFlag;
     	// normal up
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	// end left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapRightRefOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapRightRefExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapRightRefOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
     	// normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //j = lenb - 1
  //end gap
  if(stopb == lenb - 1 ){

   	const int32_t j = lenb - 1;

//   	std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << " putting in end gap " << std::endl;
//   	std::cout << "\tj:               " << j << std::endl;
//   	std::cout << "\tstopb:           " << stopb << std::endl;
//   	std::cout << "\tlenb:            " << lenb << std::endl;
//   	std::cout << "\tobjB.size() + 1: " << objB.size() + 1 << std::endl;
   	{
   		const uint32_t i = starta;
//     	if(print){
//     		   	std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << " putting in end gap " << std::endl;
//     		   	std::cout << "\ti:               " << i << std::endl;
//     		   	std::cout << "\tj:               " << j << std::endl;
//     		   	std::cout << "\tstopb:           " << stopb << std::endl;
//     		   	std::cout << "\tlenb:            " << lenb << std::endl;
//     		   	std::cout << "\tobjB.size() + 1: " << objB.size() + 1 << std::endl;
//     	}
   		if(0 == starta ){
   			//end gap
   			parts.ScoreMatrix_[i][j].leftInherit =
   			  		        parts.ScoreMatrix_[i][j - 1].leftInherit -
   			  		        parts.gapScores_.gapLeftRefExtend_;
   			parts.ScoreMatrix_[i][j].leftInheritPtr = 'L';
   		} else {
				//regular left
   			char ptrFlag;
				parts.ScoreMatrix_[i][j].leftInherit =
						needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
															parts.gapScores_.gapOpen_,
													parts.ScoreMatrix_[i][j - 1].leftInherit -
															parts.gapScores_.gapExtend_,
													parts.ScoreMatrix_[i][j - 1].diagInherit -
															parts.gapScores_.gapOpen_,
													ptrFlag);
				parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
   		}
			if(j == lastCursors.jcursor_ + 1){
				//just one over from last time, so normal diagonal
				char ptrFlag;
				//normal diag
				int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
				parts.ScoreMatrix_[i][j].diagInherit =
						match +
						needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
													parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
													parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
													ptrFlag);
				parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
			}else{
				//first row so no diag
		    parts.ScoreMatrix_[i][j].diagInheritPtr = '\0';
		    parts.ScoreMatrix_[i][j].diagInherit = std::numeric_limits<int32_t>::lowest()/2;
			}

	    //first row so no up
	    parts.ScoreMatrix_[i][j].upInheritPtr = '\0';
	    parts.ScoreMatrix_[i][j].upInherit = std::numeric_limits<int32_t>::lowest()/2;
   	}

    for (uint32_t i = starta + 1; i < stopa; ++i) {
    	  char ptrFlag;
    	  //end up
      parts.ScoreMatrix_[i][j].upInherit =
          needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapRightQueryExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	  //regular left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
      	//normal diag
      int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
      parts.ScoreMatrix_[i][j].diagInherit =
          match +
          needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                        parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
    }
  }

  //i = lena - 1, j = lenb -1

  //end
  if(stopa == lena - 1 && stopb == lenb - 1) {
  		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " putting in end gap " << std::endl;
		char ptrFlag;
		const uint32_t i = lena - 1;
		const uint32_t j = lenb - 1;
		//end up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapRightQueryExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    //end left
    parts.ScoreMatrix_[i][j].leftInherit =
        needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      parts.ScoreMatrix_[i][j - 1].leftInherit -
                          parts.gapScores_.gapRightRefExtend_,
                      parts.ScoreMatrix_[i][j - 1].diagInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  		//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  } else if(stopa == lena - 1 && stopb != lenb - 1) {
		char ptrFlag;
		const uint32_t i = stopa;
		const uint32_t j = stopb;
   	//normal up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    //end left
    parts.ScoreMatrix_[i][j].leftInherit =
        needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      parts.ScoreMatrix_[i][j - 1].leftInherit -
                          parts.gapScores_.gapRightRefExtend_,
                      parts.ScoreMatrix_[i][j - 1].diagInherit -
                          parts.gapScores_.gapRightRefOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  		//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  } else if(stopa != lena - 1 && stopb == lenb - 1) {
		char ptrFlag;
		const uint32_t i = stopa;
		const uint32_t j = stopb;
		//end up
    parts.ScoreMatrix_[i][j].upInherit =
        needleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapRightQueryExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
	  //regular left
		parts.ScoreMatrix_[i][j].leftInherit =
				needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
													parts.gapScores_.gapOpen_,
											parts.ScoreMatrix_[i][j - 1].leftInherit -
													parts.gapScores_.gapExtend_,
											parts.ScoreMatrix_[i][j - 1].diagInherit -
													parts.gapScores_.gapOpen_,
											ptrFlag);
		parts.ScoreMatrix_[i][j].leftInheritPtr = ptrFlag;
  		//normal diag
    int32_t match = parts.scoring_.mat_[objA[i - 1]][objB[j - 1]];
    parts.ScoreMatrix_[i][j].diagInherit =
        match +
        needleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].leftInherit,
                      parts.ScoreMatrix_[i - 1][j - 1].diagInherit,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].diagInheritPtr = ptrFlag;
  }

  int icursor = lena - 1;
  int jcursor = lenb - 1;
  if(lena == objA.size() + 1 && lenb == objB.size() + 1){
//  		std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//    std::cout << "\twithin icursor: " << icursor << "/" << lena << " within jcursor: " << jcursor<< "/" << lenb  << std::endl;
  	MatCursor ret{icursor, jcursor, lena, lenb};
  	ret.lastStartA_ = starta;
  	ret.lastStartB_ = startb;
  	return ret;
  } else {
    // tracerNext holds to where to go next in the matrix, will be (D) diagonal,
    // (U) up, or (L) left depending on the maximum score determined during the
    // matrix set up.
    char tracerNext = ' ';
    // get the alignment score from the  bottom right cell and set the tacer to
    // where to go next
    // keep tracing back until at the beginning of either sequence
    // Traceback algorithm follows. Score is the max of all three scores stored
    // in
    // the bottom right cell.
    // Alignments are constructed by following the correct pointer backwards at
    // each stage.
    needleMaximum(
        parts.ScoreMatrix_[icursor][jcursor].upInherit,
        parts.ScoreMatrix_[icursor][jcursor].leftInherit,
        parts.ScoreMatrix_[icursor][jcursor].diagInherit, tracerNext);
    // std::cout <<"rnv2" << std::endl;

//    std::cout << "within icursor: " << icursor << "/" << lena << " within jcursor: " << jcursor<< "/" << lenb << " tracerNext: " << tracerNext << std::endl;

    while (icursor != 0 || jcursor != 0) {
//      std::cout << "within icursor: " << icursor << "/" << lena << " within jcursor: " << jcursor<< "/" << lenb << " tracerNext: " << tracerNext << std::endl;
//				std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//				std::cout << "\twithin icursor: " << icursor << "/" << lena << " within jcursor: " << jcursor<< "/" << lenb  << std::endl;

      if (tracerNext == 'U') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
        --icursor;
      } else if (tracerNext == 'L') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
        --jcursor;
      } else if (tracerNext == 'D') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].diagInheritPtr;
        break;
      }
      // if ambigous traceback (can go either left or up), we give precedence
      // to an 'up' traceback.
      else if (tracerNext == 'B') {
        tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
        --icursor;
      }
//      if(tracerNext == '\0'){
//        std::cout << "within icursor: " << icursor << "/" << lena << " within jcursor: " << jcursor<< "/" << lenb << " tracerNext: " << tracerNext << std::endl;
//        exit(1);
//      }
//			if ('\0' == tracerNext) {
//				std::stringstream ss;
//				ss << __PRETTY_FUNCTION__ << ", error in aligning: " << "\n";
//				ss << objA << "\n";
//				ss << objB << "\n";
//				throw std::runtime_error { ss.str() };
//			}
    }
		if (0 == icursor && 0 == jcursor) {
			MatCursor ret{ static_cast<int32_t>(lena - 1),
					static_cast<int32_t>(lenb - 1), lena, lenb };
	  	ret.lastStartA_ = starta;
	  	ret.lastStartB_ = startb;
	  	return ret;
		} else {
			MatCursor ret{ icursor, jcursor, lena, lenb };
	  	ret.lastStartA_ = starta;
	  	ret.lastStartB_ = startb;
	  	return ret;
		}
  }
}

void alignCalc::runNeedleDiagonalSave(
									const std::string& objA,
									const std::string& objB,
									uint32_t inputAlignmentBlockSize,
									uint32_t inputAlignmentBlockWalkbackSize,
                   alnParts& parts) {
//  bool print = false;
//  if("AGTGAAGAAGCAATGGTGCCTTTGGAAAGTATGTAGCTCTTAGTTAAGCTTTTATCTGGAAGGTGCTGGAGAAGGTCTGGGCACAGAGCTGGTTGGCTCCATGGTTTTGTGACCTGTGCCATGATTTGAGAAGTGCCCATAGGTAGAAATCCAGGATGACTATGGGGCTCGTTTTGTATGTTTCCTGCCTCTCAAGGATTACAGACCTGTCTCTCTGTTGTCCAATGACTCTAAAGAATTGCTTCATAGATTTTGTTCAGTTTCACAGCTTTCTATAAGAGAAGTTAAATATGATACCTGTTACTCCTTCATGGCCAGAACTGGAAGCTAAGCATGGGTTTGATGTCACGCAGACCTGGGTTTGGACCCCAGCTCTGCCCCTTAATGGCTCTGTGACCTAAAGCAGGTTATTTCATCTCTTTGAGCCTCAGTGTCTTTACCTAAGAAACAGCAGGTTATTGTGAGGACCAAATAAAACAAAGCACACTTAGC" == objA &&
//  	 "TCAACATCACTAGGCAATAGGAATTTTTCAGCTCCATTTTAATCTTAAGGGACCACTGGTAATTTTAATTGGGTCTTGGATATTATATATTAAGATTTTGAGATATAATTTTAGGCTATGGGTGATTTTATTTTTTTTCCTGAGAAAGTGAAGTTTGCTTATTTTAGGCAACTGGGGTGAGGACATTAGTAATCCAAGGCCACTTTAAGCCAATCAGGTATTTATGTAATTTGAAGCTGGCCTTGATTCCCTACTATGGCTAGTTTGTTTTTAGTTTACTCTCACTCTTGGGGTATAGCCCTTCAAGTTCCCAGCCAAAAGCCTGGGGTGTTGGCCAAGGTTTGTTCTTCATTGAGGGCCCTGAGTACTAGTTTCTATCTGCCCAGATCCATAAGACTTGGTAAAAGCTCTCTTCAGCTCCTCAGCCTCTCAGCCCGGATGATGCATCAGACACACCTGGGACCCTCGGAGCTCCTGCTGTC" == objB){
//  	print = true;
//  }
//	bool print = false;
//	if ("CCCAGAGCTAGCACAGGCAGAACCCTCATGTGAATCAAAAAGGACACCACTCCGTCCCTGCTAACCCTGTCCTCTGTCTCTTCATTTTCCACCTAGACTGATATTGCCCCCAAACATTAAGCTTTCTTGTATTACCTTCATGAT"
//			== objB
//			&& "TATGTATACCTTCCTAAGATTAGCGGTGTCCTCTTGATTGGGGGATGGTTATAAATAGCTTTTAAGGTGATTTCCAGGCCTAACATTCTATTATTGACTTATAACTGCACTTTAAAGGCTTAAACACACACAATCCCATTGAATCCTCAAGATACTTCCATTGGCAGGCAAGTATTTTCCTTG"
//					== objA) {
//		print = true;
//	}
//	bool print = false;
//	if ("CCCAGAGCTAGCACAGGCAGAACCCTCATGTGAATCAAAAAGGACACCACTCCGTCCCTGCTAACCCTGTCCTCTGTCTCTTCATTTTCCACCTAGACTGATATTGCCCCCAAACATTAAGCTTTCTTGTATTACCTTCATGAT"
//			== objB
//			&& "TATGTATACCTTCCTAAGATTAGCGGTGTCCTCTTGATTGGGGGATGGTTATAAATAGCTTTTAAGGTGATTTCCAGGCCTAACATTCTATTATTGACTTATAACTGCACTTTAAAGGCTTAAACACACACAATCCCATTGAATCCTCAAGATACTTCCATTGGCAGGCAAGTATTTTCCTTG"
//					== objA) {
//		print = true;
//	}



  parts.gHolder_.gapInfos_.clear();
  parts.gHolder_.addFromFile_ = false;
//  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	auto cursors = runNeedleDiagonalSaveInit(objA, objB, inputAlignmentBlockSize, parts);
//  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	std::set<std::pair<uint32_t, uint32_t>> previousCursors;
	uint32_t internalAlignmentBlockSize  = inputAlignmentBlockSize;
	uint32_t internalAlignmentBlockWalkbackSize  = inputAlignmentBlockWalkbackSize;
	previousCursors.emplace(std::make_pair(cursors.icursor_, cursors.jcursor_));
//  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	while((    cursors.icursor_ != static_cast<int64_t>(objA.size()) && cursors.jcursor_ != static_cast<int64_t>(objB.size())) ||
			  (cursors.icursor_ == static_cast<int64_t>(objA.size()) && cursors.jcursor_ != static_cast<int64_t>(objB.size())) ||
			  (cursors.icursor_ != static_cast<int64_t>(objA.size()) && cursors.jcursor_ == static_cast<int64_t>(objB.size()))) {
//		std::cout << "icursor: " << cursors.icursor_ << " jcursor: "<< cursors.jcursor_<< " objA.size(): " << objA.size() << " objB.size(): "<< objB.size() << std::endl;

//		if(print){
//			std::cout << "icursor: " << cursors.icursor_ << " jcursor: "<< cursors.jcursor_<< " objA.size(): " << objA.size() << " objB.size(): "<< objB.size() << std::endl;
//		}
//	  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

		auto nextCursor = runNeedleDiagonalSaveStep(objA, objB, internalAlignmentBlockSize, internalAlignmentBlockWalkbackSize, cursors, parts);;
//	  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

		if(njh::in(std::make_pair(cursors.icursor_, cursors.jcursor_), previousCursors)){
			internalAlignmentBlockSize += internalAlignmentBlockSize/4;
			//internalAlignmentBlockWalkbackSize += internalAlignmentBlockWalkbackSize/4;
		}
//		if(nextCursor.icursor_ == cursors.icursor_ && nextCursor.jcursor_ == cursors.jcursor_){
//			nextCursor.icursor_ = std::min<int32_t>(objA.size(), nextCursor.icursor_ + internalAlignmentBlockSize);
//			nextCursor.jcursor_ = std::min<int32_t>(objB.size(), nextCursor.jcursor_ + internalAlignmentBlockSize);
//		}
//		if(print){
//			std::cout << "\tlastCursor.icursor: " << cursors.icursor_ << " lastCursor.jcursor: " << cursors.jcursor_  << std::endl;
//			std::cout << "\tnextCursor.icursor: " << nextCursor.icursor_ << " nextCursor.jcursor: " << nextCursor.jcursor_  << std::endl;
//		}
//		std::cout << "\tlastCursor.icursor: " << cursors.icursor_ << " lastCursor.jcursor: " << cursors.jcursor_  << std::endl;
//		std::cout << "\tnextCursor.icursor: " << nextCursor.icursor_ << " nextCursor.jcursor: " << nextCursor.jcursor_  << std::endl;
		cursors = nextCursor;
		previousCursors.emplace(std::make_pair(cursors.icursor_, cursors.jcursor_));
//		std::cout << "icursor: " << cursors.icursor_ << " jcursor: "<< cursors.jcursor_<< std::endl;



//		std::cout << std::endl;
	}
//  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

//	std::cout << "icursor: " << cursors.icursor_ << "/" << objA.size() << " jcursor: "<< cursors.jcursor_ << "/" << objB.size()<< std::endl;
//
//	std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  int icursor = cursors.icursor_;
  int jcursor = cursors.jcursor_;

  // tracerNext holds to where to go next in the matrix, will be (D) diagonal,
  // (U) up, or (L) left depending on the maximum score determined during the
  // matrix set up.
  char tracerNext = ' ';

  // get the alignment score from the  bottom right cell and set the tacer to
  // where to go next
  // keep tracing back until at the begining of either sequence
  // Traceback algorithm follows. Score is the max of all three scores stored
  // in
  // the bottom right cell.
  // Alignments are constructed by following the correct pointer backwards at
  // each stage.
  parts.score_ = needleMaximum(
      parts.ScoreMatrix_[icursor][jcursor].upInherit,
      parts.ScoreMatrix_[icursor][jcursor].leftInherit,
      parts.ScoreMatrix_[icursor][jcursor].diagInherit, tracerNext);

//  std::cout << "final icursor: " << icursor << "/" << objA.size() << " final jcursor: " << jcursor<< "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor][jcursor].upInherit: " << parts.ScoreMatrix_[icursor][jcursor].upInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor][jcursor].leftInherit: " << parts.ScoreMatrix_[icursor][jcursor].leftInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor][jcursor].diagInherit: " << parts.ScoreMatrix_[icursor][jcursor].diagInherit << std::endl;
//
//  std::cout << "final icursor up and over 1: " << icursor - 1 << "/" << objA.size() << " final jcursor: " << jcursor -1 << "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor - 1].upInherit: "   << parts.ScoreMatrix_[icursor - 1][jcursor - 1].upInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor - 1].leftInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor - 1].leftInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor - 1].diagInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor - 1].diagInherit << std::endl;
//
//  std::cout << "final icursor up 1: " << icursor << "/" << objA.size() << " final jcursor: " << jcursor -1 << "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor][jcursor - 1].upInherit: "   << parts.ScoreMatrix_[icursor][jcursor - 1].upInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor][jcursor - 1].leftInherit: " << parts.ScoreMatrix_[icursor][jcursor - 1].leftInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor][jcursor - 1].diagInherit: " << parts.ScoreMatrix_[icursor][jcursor - 1].diagInherit << std::endl;
//
//  std::cout << "final icursor over 1: " << icursor - 1 << "/" << objA.size() << " final jcursor: " << jcursor << "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor].upInherit: "   << parts.ScoreMatrix_[icursor - 1][jcursor].upInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor].leftInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor].leftInherit << std::endl;
//  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor].diagInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor].diagInherit << std::endl;
//
//  std::cout << std::endl;
//
//  if(objA == "ATATATTTGGACAGCAATGAAACATGGTGCAGAAATGAATAGTACTACGTGTAATTGTAGTGGTGATAGTAGTAGTGGTGAAAATCAAACTAATAGTTGTGATGATATGCCTACGATTGATTTGATCCCTCAATATTTACGGTTTTTGCAAGAATGGGTAGAACATTTTTGTAAACAACGTCAAGGAAAAGTAAAAGATGTGATAGAGAATTGTAAGTCGTGTAAGAATACATCGGGTGAAAGAATAATTGGAGGTACATGTAACGGTGAGTGTAAAACTGAATGTGAAAAAAAATGTAAAAATAAATGTGAAGCATACAAAACATTTATTGAAAAGTTTTGTACAGCTGATGGTGGTACTTCCGGATCTCCATGGAGCAAAAGGTGGTACCAAATATATATGAGGTATTCCAAATATA" &&
//  		objB == "ATATATTTGGTTAGCAATGAAACATGGTGCAGAAATGAATAGTACTATGTGTAATGCTGATGGTAGTGTCACTGGTAGTAGTGATAGTGGTAGTACTACGTGTAGTGGTGACAATGGTAGTATTAGTTGTGATGATATTCCTACGATTGATTTGATCCCTCAATATTTACGGTTTTTACAAGAATGGATAGAACATTTTTGTAAACAACGTCAAGAAAAAGTAAATGCTGTGATAGAGAATTGTAATTCGTGTAAAAATACCTCAAGTAAAACAAAAATTGGAGACACATGTGGTAATGATTGTAAAACTAAATGTAAAGTAGCGTGTGACGCATACAAAACATTTATTGAAAAGTGTACGGGTGGTGGTACTGGTACTGCCGGATCCTCATGGGTGAAAAGGTGGTACCAAATATATATGAGGTATTCCAAATATA"){
//    std::ofstream outScoresInfo("allScores_alnCalc.txt");
//    outScoresInfo << "icursor\tjcursor\tup\tleft\tdiag" << std::endl;
//    for(uint32_t i = 0; i <= icursor; ++i){
//    		for(uint32_t j = 0; j <= jcursor; ++j){
//    			outScoresInfo << i
//    					<< "\t" << j
//  					<< "\t" << parts.ScoreMatrix_[i][j].upInherit
//  					<< "\t" << parts.ScoreMatrix_[i][j].leftInherit
//  					<< "\t" << parts.ScoreMatrix_[i][j].diagInherit
//  					<< std::endl;
//    		}
//    }
//  }
//
//  if(print){
//  	  std::cout << "final icursor: " << icursor << "/" << objA.size() << " final jcursor: " << jcursor<< "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor][jcursor].upInherit: " << parts.ScoreMatrix_[icursor][jcursor].upInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor][jcursor].leftInherit: " << parts.ScoreMatrix_[icursor][jcursor].leftInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor][jcursor].diagInherit: " << parts.ScoreMatrix_[icursor][jcursor].diagInherit << std::endl;
//
//  	  std::cout << "final icursor up and over 1: " << icursor - 1 << "/" << objA.size() << " final jcursor: " << jcursor -1 << "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor - 1].upInherit: "   << parts.ScoreMatrix_[icursor - 1][jcursor - 1].upInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor - 1].leftInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor - 1].leftInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor - 1].diagInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor - 1].diagInherit << std::endl;
//
//  	  std::cout << "final icursor up 1: " << icursor << "/" << objA.size() << " final jcursor: " << jcursor -1 << "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor][jcursor - 1].upInherit: "   << parts.ScoreMatrix_[icursor][jcursor - 1].upInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor][jcursor - 1].leftInherit: " << parts.ScoreMatrix_[icursor][jcursor - 1].leftInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor][jcursor - 1].diagInherit: " << parts.ScoreMatrix_[icursor][jcursor - 1].diagInherit << std::endl;
//
//  	  std::cout << "final icursor over 1: " << icursor - 1 << "/" << objA.size() << " final jcursor: " << jcursor << "/" << objB.size() << " tracerNext: " << tracerNext << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor].upInherit: "   << parts.ScoreMatrix_[icursor - 1][jcursor].upInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor].leftInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor].leftInherit << std::endl;
//  	  std::cout << "parts.ScoreMatrix_[icursor - 1][jcursor].diagInherit: " << parts.ScoreMatrix_[icursor - 1][jcursor].diagInherit << std::endl;
//
//  	  std::cout << std::endl;
//
//  	  if(print){
//
//  	  	std::ofstream outScoresInfo("allScores_alnCalc.txt");
//  	    outScoresInfo << "icursor\tjcursor\tup\tleft\tdiag" << std::endl;
//  	    for(uint32_t i = 0; i <= icursor; ++i){
//  	    		for(uint32_t j = 0; j <= jcursor; ++j){
//  	    			outScoresInfo << i
//  	    					<< "\t" << j
//  	  					<< "\t" << parts.ScoreMatrix_[i][j].upInherit
//  	  					<< "\t" << parts.ScoreMatrix_[i][j].leftInherit
//  	  					<< "\t" << parts.ScoreMatrix_[i][j].diagInherit
//  	  					<< std::endl;
//  	    		}
//  	    }
//  	  }
//  }


  parts.gHolder_.score_ = parts.score_;
  uint32_t gapBSize = 0;
  uint32_t gapASize = 0;
//  std::cout << std::endl;
//  std::cout << objA << std::endl;
//  std::cout << objB << std::endl;

//  std::ofstream outInfo("alignmentPath_alnCalc.txt");
//  outInfo << "icursor\tjcursor\tscore" << std::endl;
  // std::cout <<"rnv2" << std::endl;
  //int32_t score = parts.score_;
  while (icursor != 0 || jcursor != 0) {
//  	std::cout << "icursor: " << icursor << " jcursor: " << jcursor << " score: " << score << " tracerNext: " << tracerNext << std::endl;
//  	if(print){
//  		std::cout << "icursor: " << icursor << " jcursor: " << jcursor << " score: " << score << " tracerNext: " << tracerNext << std::endl;
//  	}
  		//outInfo << icursor << '\t' << jcursor << "\t" << score << std::endl;
    if (tracerNext == 'U') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      //score = parts.ScoreMatrix_[icursor][jcursor].upInherit;
      if (tracerNext != 'U' && tracerNext != 'B') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    } else if (tracerNext == 'L') {
      ++gapASize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
      //score = parts.ScoreMatrix_[icursor][jcursor].leftInherit;
      if (tracerNext != 'L') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(icursor, gapASize, true));
        gapASize = 0;
      }
      --jcursor;
    } else if (tracerNext == 'D') {
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].diagInheritPtr;
      //score = parts.ScoreMatrix_[icursor][jcursor].diagInherit;
      --icursor;
      --jcursor;
    }
    // if ambigous traceback (can go either left or up), we give precedence
    // to an 'up' traceback.
    else if (tracerNext == 'B') {
      ++gapBSize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].upInheritPtr;
      //score = parts.ScoreMatrix_[icursor][jcursor].upInherit;
      if (tracerNext != 'U' && tracerNext != 'B') {
        parts.gHolder_.gapInfos_.emplace_back(
        		njhseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    }
    if(tracerNext == '\0'){
    	std::stringstream ss;
    	ss << __PRETTY_FUNCTION__ << ", error in aligning: " << "\n";
    	ss << objA << "\n";
    	ss << objB << "\n";
    	throw std::runtime_error{ss.str()};
    }
    if(icursor < 0 || jcursor < 0){
    		break;
    }
  }
  if ((tracerNext == 'U' || tracerNext == 'B') && gapBSize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(njhseq::gapInfo(jcursor, gapBSize, false));
  } else if (tracerNext == 'L' && gapASize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(njhseq::gapInfo(icursor, gapASize, true));
  }
  //std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
}


}  // namespace njhseq

