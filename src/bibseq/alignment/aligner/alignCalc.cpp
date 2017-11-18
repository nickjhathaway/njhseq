//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//

/*
 * alnCalc.cpp
 *
 *  Created on: Nov 29, 2015
 *      Author: nick
 */


#include "alignCalc.hpp"


namespace bibseq {

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
  for (uint32_t i = 1; i < parts.maxSize_; i++) {
    parts.ScoreMatrix_[i][0].upInherit = 0;
    parts.ScoreMatrix_[i][0].leftInherit = 0;
    parts.ScoreMatrix_[i][0].diagInherit = 0;
    parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
    parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
    parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
  }
  // initialize first row:
  for (uint32_t j = 1; j < parts.maxSize_; j++) {
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
      std::cout << "ERROR!!!!!" << std::endl;
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
		for (uint32_t j = 2; j < lenb - 1; j++) {
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
    for (uint32_t i = 2; i < lena - 1 ; i++) {
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
  for (uint32_t i = 2; i < lena - 1 ; i++) {
    for (uint32_t j = 2; j < lenb - 1 ; j++) {
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
    for (uint32_t i = 2; i < lena - 1; i++) {
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
        		bibseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    } else if (tracerNext == 'L') {
      ++gapASize;
      tracerNext = parts.ScoreMatrix_[icursor][jcursor].leftInheritPtr;
      if (tracerNext != 'L') {
        parts.gHolder_.gapInfos_.emplace_back(
        		bibseq::gapInfo(icursor, gapASize, true));
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
        		bibseq::gapInfo(jcursor, gapBSize, false));
        gapBSize = 0;
      }
      --icursor;
    }
  }
  if ((tracerNext == 'U' || tracerNext == 'B') && gapBSize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(bibseq::gapInfo(jcursor, gapBSize, false));
  } else if (tracerNext == 'L' && gapASize != 0) {
    parts.gHolder_.gapInfos_.emplace_back(bibseq::gapInfo(icursor, gapASize, true));
  }
}

}  // namespace bibseq

