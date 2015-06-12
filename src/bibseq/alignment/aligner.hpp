#pragma once
#include "bibseq/objects/seqObjects/baseReadObject.hpp"
#include "bibseq/objects/helperObjects/kmer.hpp"
#include "bibseq/objects/helperObjects/gaps.hpp"
#include "bibseq/objects/helperObjects/tandemRepeat.hpp"
#include "bibseq/objects/helperObjects/mismatch.hpp"
#include "bibseq/alignment/alignerUtils.h"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/alignment/alnInfoHolder.hpp"
#include "bibseq/programUtils/runningParameters.hpp"

namespace bibseq {

struct scoreMatrixCell {
  int32_t upInherit;
  int32_t leftInherit;
  int32_t diagInherit;
  // for traceback: 'U' = up, 'L' = Left, 'D' = diagonal, 'B' either up or left
  char upInheritPtr;
  char leftInheritPtr;
  char diagInheritPtr;
};


struct alnParts {

  alnParts(uint32_t maxSize, const gapScoringParameters& gapScores,
           const substituteMatrix& scoring)
      : maxSize_(maxSize + 10),
        gapScores_(gapScores),
        scoring_(scoring),
        ScoreMatrix_(std::vector<std::vector<scoreMatrixCell>>(
            maxSize_, std::vector<scoreMatrixCell>(maxSize_))) {}
  alnParts()
      : maxSize_(400),
        gapScores_(gapScoringParameters()),
        scoring_(substituteMatrix(2, -2)),
        ScoreMatrix_(std::vector<std::vector<scoreMatrixCell>>(
            maxSize_, std::vector<scoreMatrixCell>(maxSize_))) {}

  int32_t score_ = 0;
  uint32_t maxSize_;
  // gap scores
  gapScoringParameters gapScores_;
  substituteMatrix scoring_;
  alnInfoGlobal gHolder_;
  alnInfoLocal lHolder_;
  // the matrix
  std::vector<std::vector<scoreMatrixCell>> ScoreMatrix_;
};

class alignCalc {
 public:
  // score matrix maximum functions
  /*! \brief Needle Maximum
   *
   *
   *  Function is used to find the maximum of an alignment cell using needle
   *(global) style alignment
   */
  static int32_t needleMaximum(int32_t u, int32_t l, int32_t d, char& p) {
    if (u == l && u == d) {
      p = 'B';
      return u;
    } else if (u >= l && u >= d) {
      if (u == l) {
        p = 'B';
        return u;
      } else {
        p = 'U';
        return u;
      }
    } else if (l >= u && l >= d) {
      if (l == u) {
        p = 'B';
        return l;
      } else {
        p = 'L';
        return l;
      }
    } else {
      p = 'D';
      return d;
    }
  }
  /*! \brief Smith Maximum
   *
   *
   *  Function is used to find the maximum of an alignment cell using smith
   *(local) style alignment
   */
  // template <typename NUM>
  static int32_t smithMaximum(int32_t u, int32_t l, int32_t d, char& p) {
    if (u < 0 && l < 0 && d < 0) {
      p = '\0';
      return 0;
    }
    if (u == l && u == d) {
      p = 'B';
      return u;
    } else if (u >= l && u >= d) {
      if (u == l) {
        p = 'B';
        return u;
      } else {
        p = 'U';
        return u;
      }
    } else if (l >= u && l >= d) {
      if (l == u) {
        p = 'B';
        return l;
      } else {
        p = 'L';
        return l;
      }
    } else {
      p = 'D';
      return d;
    }
  }
  static void runNeedleSave(const std::string& objA, const std::string& objB,
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
    for (int32_t i = 1; i < parts.maxSize_; i++) {
      if (i == 1) {
        parts.ScoreMatrix_[i][0].upInherit = - parts.gapScores_.gapLeftOpen_;
        parts.ScoreMatrix_[i][0].leftInherit = 0;
        parts.ScoreMatrix_[i][0].diagInherit = 0;
        parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
        parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
        parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
      } else {
        parts.ScoreMatrix_[i][0].upInherit =
            parts.ScoreMatrix_[i - 1][0].upInherit -
            parts.gapScores_.gapLeftExtend_;
        parts.ScoreMatrix_[i][0].leftInherit = 0;
        parts.ScoreMatrix_[i][0].diagInherit = 0;
        parts.ScoreMatrix_[i][0].upInheritPtr = 'U';
        parts.ScoreMatrix_[i][0].leftInheritPtr = '\0';
        parts.ScoreMatrix_[i][0].diagInheritPtr = '\0';
      }
    }
    // initialize first row:
    for (int32_t j = 1; j < parts.maxSize_; j++) {
      if (j == 1) {
        parts.ScoreMatrix_[0][j].upInherit = 0;
        parts.ScoreMatrix_[0][j].leftInherit = - parts.gapScores_.gapLeftOpen_;
        parts.ScoreMatrix_[0][j].diagInherit = 0;
        parts.ScoreMatrix_[0][j].upInheritPtr = '\0';
        parts.ScoreMatrix_[0][j].leftInheritPtr = 'L';
        parts.ScoreMatrix_[0][j].diagInheritPtr = '\0';
      } else {
        parts.ScoreMatrix_[0][j].upInherit = 0;
        parts.ScoreMatrix_[0][j].leftInherit =
            parts.ScoreMatrix_[0][j - 1].leftInherit -
            parts.gapScores_.gapLeftExtend_;
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
                      parts.gapScores_.gapRightOpen_;
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
          parts.gapScores_.gapRightOpen_;
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
                              parts.gapScores_.gapRightOpen_,
                          parts.ScoreMatrix_[i][j - 1].leftInherit -
                              parts.gapScores_.gapRightExtend_,
                          parts.ScoreMatrix_[i][j - 1].diagInherit -
                              parts.gapScores_.gapRightOpen_,
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
                              parts.gapScores_.gapRightExtend_,
                          parts.ScoreMatrix_[i - 1][j].leftInherit -
                              parts.gapScores_.gapRightOpen_,
                          parts.ScoreMatrix_[i - 1][j].diagInherit -
                              parts.gapScores_.gapRightOpen_,
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
                            parts.gapScores_.gapRightExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapRightOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapRightOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//end left
      parts.ScoreMatrix_[i][j].leftInherit =
          needleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
                            parts.gapScores_.gapRightOpen_,
                        parts.ScoreMatrix_[i][j - 1].leftInherit -
                            parts.gapScores_.gapRightExtend_,
                        parts.ScoreMatrix_[i][j - 1].diagInherit -
                            parts.gapScores_.gapRightOpen_,
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

  static void runSmithSave(const std::string& objA, const std::string& objB,
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
    for (uint32_t i = 1; i < lena; i++) {
      for (uint32_t j = 1; j < lenb; j++) {
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

  template <typename T, typename F>
  static void rearrangeGlobal(T& readA, T& readB, const F& fill,
                              const alnInfoGlobal& gHolder) {
    for (const auto& g : gHolder.gapInfos_) {
      if (g.gapInA_) {
        readA.insert(readA.begin() + g.pos_, g.size_, fill);
      } else {
        readB.insert(readB.begin() + g.pos_, g.size_, fill);
      }
    }
  }
  template <typename T, typename F>
  static void rearrangeLocal(T& readA, T& readB, const F& fill,
                             const alnInfoLocal& lHolder) {
    readA = getSubVector(readA, lHolder.localAStart_, lHolder.localASize_);
    readB = getSubVector(readB, lHolder.localBStart_, lHolder.localBSize_);
    for (const auto& g : lHolder.gapInfos_) {
      if (g.gapInA_) {
        readA.insert(readA.begin() + g.pos_ - lHolder.localAStart_, g.size_,
                     fill);
      } else {
        readB.insert(readB.begin() + g.pos_ - lHolder.localBStart_, g.size_,
                     fill);
      }
    }
  }
};

/*! \brief Aligner Class
 *
 *
 *  This class can do local or global alignment using simple scoring or using a
 *  provide scoring matrix
 */
class aligner {

 public:
  // constructors
	/*
	 * \brief Default aligner, can handle alignments of up to 400 bps
	 */
  aligner() : parts_(alnParts()),
  alnHolder_(alnInfoMasterHolder(parts_.gapScores_, parts_.scoring_)){
  	setDefaultQualities();
  }

	aligner(uint32_t maxSize, const gapScoringParameters& gapPars,
          const substituteMatrix& scoreMatrix)
      : parts_(maxSize, gapPars, scoreMatrix),
        alnHolder_(alnInfoMasterHolder(gapPars, scoreMatrix)) {
		setDefaultQualities();
	}
	aligner(uint32_t maxSize, const gapScoringParameters& gapPars,
          const substituteMatrix& scoreMatrix, bool countEndGaps)
      : parts_(maxSize, gapPars, scoreMatrix),
        alnHolder_(alnInfoMasterHolder(gapPars, scoreMatrix)), countEndGaps_(countEndGaps) {
		setDefaultQualities();
	}

  /*! \brief Constructor with Scoring map
   *
   *
   *  This constructor sets up the scoring for a determined by a map of scores
   *allowing for more complex scoring set up
   */
  aligner(uint32_t maxSize, const gapScoringParameters & gapPars,
  				const substituteMatrix& subMatrix,
          const kmerMaps& inKmaps, uint32_t primaryQuality, uint32_t secondaryQuality,
          uint32_t qualThresholdWindow, bool countEndGaps
          ):alignObjectA_(baseReadObject()), alignObjectB_(baseReadObject()),
          		parts_(maxSize, gapPars, subMatrix),alnHolder_(gapPars, subMatrix),
          		primaryQual_(primaryQuality), secondaryQual_(secondaryQuality),
          		kMaps_(inKmaps),
          		countEndGaps_(countEndGaps){
  	setStartingParamters(inKmaps, primaryQuality, secondaryQuality, qualThresholdWindow, countEndGaps);
  }



  // to hold the sequence alignments
  baseReadObject alignObjectA_;
  baseReadObject alignObjectB_;

  alnParts parts_;
  alnInfoMasterHolder alnHolder_;

  uint32_t numberOfAlingmentsDone_ = 0;

  // Aligner
  void alignScore(const std::string& firstSeq,
  									 const std::string& secondSeq,
                    bool local) {
    if (local) {
    	alignCalc::runSmithSave(firstSeq, secondSeq, parts_);
    } else {
    	alignCalc::runNeedleSave(firstSeq, secondSeq, parts_);
    }
    ++numberOfAlingmentsDone_;
  }

  void alignScoreCache(const std::string& firstSeq, const std::string& secondSeq,
                bool local) {
    if (local) {
    	if(alnHolder_.localHolder_[parts_.gapScores_.uniqueIdentifer_]
            .getAlnInfo(firstSeq, secondSeq, parts_.lHolder_)){
    		parts_.score_ = parts_.lHolder_.score_;
    	}else{
    		alignCalc::runSmithSave(firstSeq, secondSeq, parts_);
    		alnHolder_.localHolder_[parts_.gapScores_.uniqueIdentifer_].addAlnInfo(firstSeq, secondSeq, parts_.lHolder_);
    		++numberOfAlingmentsDone_;
    	}
    } else {
    	if(alnHolder_.globalHolder_[parts_.gapScores_.uniqueIdentifer_]
    	            .getAlnInfo(firstSeq, secondSeq, parts_.gHolder_)){
    		parts_.score_ = parts_.gHolder_.score_;
			}else{
				alignCalc::runNeedleSave(firstSeq, secondSeq, parts_);
				alnHolder_.globalHolder_[parts_.gapScores_.uniqueIdentifer_].addAlnInfo(firstSeq, secondSeq, parts_.gHolder_);
				++numberOfAlingmentsDone_;
			}
    }
  }
  void alignVec(const baseReadObject & ref, const baseReadObject & read, bool local){
  	alignScoreCache(ref.seqBase_.seq_, read.seqBase_.seq_, local);
  	rearrangeObjs(ref.seqBase_, read.seqBase_, local);
  }
  void alignVec(const seqInfo & ref, const seqInfo & read, bool local){
  	alignScoreCache(ref.seq_, read.seq_, local);
  	rearrangeObjs(ref, read, local);
  }
  void alignReg(const baseReadObject & ref, const baseReadObject & read, bool local){
  	alignScore(ref.seqBase_.seq_, read.seqBase_.seq_, local);
  	rearrangeObjs(ref.seqBase_, read.seqBase_, local);
  }
  void alignReg(const seqInfo & ref, const seqInfo & read, bool local){
  	alignScore(ref.seq_, read.seq_, local);
  	rearrangeObjs(ref, read, local);
  }
  std::pair<uint32_t, uint32_t> findReversePrimer(const std::string& read,
                                          				const std::string& primer){
  	alignScoreCache(read, primer, true);
  	return {parts_.lHolder_.localAStart_,
  		parts_.lHolder_.localAStart_ + parts_.lHolder_.localASize_ - 1};
  }
  std::pair<uint32_t, uint32_t> findReversePrimer(const baseReadObject& read,
                                          				const baseReadObject& primer){
  	return findReversePrimer(read.seqBase_.seq_, primer.seqBase_.seq_);
  }

  void rearrangeSeq(const std::string& firstRead,
  		const std::string& secondRead, bool local) {
    alignObjectA_.seqBase_ = seqInfo("A", firstRead);
    alignObjectB_.seqBase_ = seqInfo("B", secondRead);;
    if(local){
    	alignCalc::rearrangeLocal(alignObjectA_.seqBase_.seq_, alignObjectB_.seqBase_.seq_, '-',
    	                                      parts_.lHolder_);
    	alignCalc::rearrangeLocal(alignObjectA_.seqBase_.qual_, alignObjectB_.seqBase_.qual_, 0,
    	                                      parts_.lHolder_);
    }else{
    	alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.seq_, alignObjectB_.seqBase_.seq_, '-',
    	                                      parts_.gHolder_);
    	alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.qual_, alignObjectB_.seqBase_.qual_, 0,
    	                                      parts_.gHolder_);
    }
  }

  void rearrangeObjs(const seqInfo& firstRead,
  		const seqInfo& secondRead, bool local) {
    alignObjectA_.seqBase_ = firstRead;
    alignObjectB_.seqBase_ = secondRead;
    if(local){
    	alignCalc::rearrangeLocal(alignObjectA_.seqBase_.seq_, alignObjectB_.seqBase_.seq_, '-',
    	                                      parts_.lHolder_);
    	alignCalc::rearrangeLocal(alignObjectA_.seqBase_.qual_, alignObjectB_.seqBase_.qual_, 0,
    	                                      parts_.lHolder_);
    }else{
    	alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.seq_, alignObjectB_.seqBase_.seq_, '-',
    	                                      parts_.gHolder_);
    	alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.qual_, alignObjectB_.seqBase_.qual_, 0,
    	                                      parts_.gHolder_);
    }
  }
  void rearrangeObjs(const baseReadObject& firstRead,
  		const baseReadObject& secondRead, bool local) {
  	rearrangeObjs(firstRead.seqBase_,secondRead.seqBase_, local );
  }




  //////////alignment
  // various scores

  comparison comp_;

  // tandem repeats

  std::map<uint32_t, mismatch> mismatches_;
  std::map<uint32_t, mismatch> lowKmerMismatches_;
  // the indels
  std::map<uint32_t, gap> alignmentGaps_;
  uint32_t primaryQual_;
  uint32_t secondaryQual_;

  uint32_t qualThresWindow_;
  kmerMaps kMaps_;
  bool countEndGaps_;


  void scoreAlignment(bool editTheSame);
  bool& CountEndGaps() { return countEndGaps_; }
  void setQual(int primaryQual, int secondayQual, int qualThresWindow);

  // kmer map setting
  void setKmerMpas(const kmerMaps& inKmerMaps) { kMaps_ = inKmerMaps; }
  kmerMaps getKmerMaps() { return kMaps_; }
  /*! \brief Find Reverse Primer
   *
   *
   *  Returns the start and end position of the highest matching local alignment
   *of the reverse primer to the sequence
   */
 // std::pair<int, int> findReversePrimer(const baseReadObject& read,
                                    //    const baseReadObject& primer);
  // Outputting
  void outPutParameterInfo(std::ostream& out) const;


  void handleGapCountingInA(gap& currentGap, bool weighHomopolymers);
  void handleGapCountingInB(gap& currentGap, bool weighHomopolymers);
  // profile with low kmer checking
  void profileAlignment(const baseReadObject& objectA,
                        const baseReadObject& objectB, int kLength,
                        bool kmersByPosition, bool checkKmer, bool usingQuality,
                        bool doingMatchQuality, bool weighHomopolyer,
                        uint32_t start = 0, uint32_t stop = 0);
  void profileAlignment(const seqInfo& objectA, const seqInfo& objectB,
                        int kLength, bool kmersByPosition, bool checkKmer,
                        bool usingQuality, bool doingMatchQuality,
                        bool weighHomopolyer, uint32_t start = 0,
                        uint32_t stop = 0);

  void profilePrimerAlignment(const baseReadObject& objectA,
                              const baseReadObject& objectB,
                              bool weighHomopolymers);
  void profilePrimerAlignment(const seqInfo& objectA,
                              const seqInfo& objectB,
                              bool weighHomopolymers);

  comparison compareAlignment(const baseReadObject& objectA,
                                         const baseReadObject& objectB,
                                         const runningParameters& runParams,
                                         bool checkKmers, bool kmersByPosition,
                                         bool weighHomopolymers);
  comparison compareAlignment(const seqInfo& objectA,
                                         const seqInfo& objectB,
                                         const runningParameters& runParams,
                                         bool checkKmers, bool kmersByPosition,
                                         bool weighHomopolymers);



  // constructor helper
  void setStartingParamters(const kmerMaps& inKmaps, uint32_t primaryQuality,
  		uint32_t secondaryQuality, uint32_t qualThresholdWindow,
                            bool countEndGaps);

  void setDefaultQualities() {
    primaryQual_ = 20;
    secondaryQual_ = 15;
    qualThresWindow_ = 5;
  };
  void resetCounts();
  void resetAlignmentInfo();
  bool checkForTandemRepeatGap();
  static bool checkTwoEqualSeqs(const std::string& seq1,
                                const std::string& seq2,
                                int allowableMismatches);
  // finding tandem repeats
  std::vector<tandemRepeat> findTandemRepeatsInSequence(
      const std::string& str, int match = 2, int mismatch = -2, int gap = -7,
      int minimumAlignScore = 50);
  tandemRepeat findTandemRepeatOfStrInSequence(std::string str,
                                               std::string tandem,
                                               int match = 2, int mismatch = -2,
                                               int gap = -7,
                                               int minimumAlignScore = 50);
  tandemRepeat findTandemRepeatOfStrInSequenceDegen(
      std::string str, std::string tandem, int match = 2, int mismatch = -2,
      int gap = -7, int minimumAlignScore = 50);
  static bool checkTwoStringsDegen(
      const std::string& str1, const std::string& str2, int allowableMismatches,
      const substituteMatrix& scoringArray);
  uint32_t getAlignPosForSeqAPos(uint32_t seqAPos);
  uint32_t getAlignPosForSeqBPos(uint32_t seqBPos);
 public:
  void setGeneralScorring(int32_t generalMatch, int32_t generalMismatch){
  	parts_.scoring_ = substituteMatrix(generalMatch, generalMismatch);
  }

  //void setScoringArrayWithMap(std::unordered_map<
    //  char, std::unordered_map<char, int>> matchMatrixScoresMap);
  //void setScoringArrayWithMap(const int matchMatrix[4][4]);
};



}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "aligner.cpp"
#endif
