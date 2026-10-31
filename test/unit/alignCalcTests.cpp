#include <catch_amalgamated.hpp>
#include "njhseq/alignment/aligner/alignCalc.hpp"
#include <random>

using namespace njhseq;

namespace {

//reference copies of needleMaximum and runNeedleSave in their straightforward form, alignCalc's versions must give exactly the same results
  int32_t referenceNeedleMaximum(int32_t u, int32_t l, int32_t d, char& p) {
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

void referenceRunNeedleSave(const std::string& objA, const std::string& objB,
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
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
        referenceNeedleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
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
        referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//end left
      parts.ScoreMatrix_[i][j].leftInherit =
          referenceNeedleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                            parts.gapScores_.gapRightQueryExtend_,
                        parts.ScoreMatrix_[i - 1][j].leftInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        parts.ScoreMatrix_[i - 1][j].diagInherit -
                            parts.gapScores_.gapRightQueryOpen_,
                        ptrFlag);
      parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
    	//regular left
      parts.ScoreMatrix_[i][j].leftInherit =
          referenceNeedleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
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
          referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
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
        referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j].upInherit -
                          parts.gapScores_.gapRightQueryExtend_,
                      parts.ScoreMatrix_[i - 1][j].leftInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      parts.ScoreMatrix_[i - 1][j].diagInherit -
                          parts.gapScores_.gapRightQueryOpen_,
                      ptrFlag);
    parts.ScoreMatrix_[i][j].upInheritPtr = ptrFlag;
  	//end left
    parts.ScoreMatrix_[i][j].leftInherit =
        referenceNeedleMaximum(parts.ScoreMatrix_[i][j - 1].upInherit -
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
        referenceNeedleMaximum(parts.ScoreMatrix_[i - 1][j - 1].upInherit,
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
  parts.score_ = referenceNeedleMaximum(
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

//fill the matrix with junk so a run that reads a cell it didn't write in that run gives a different result
void fillWithJunk(alnParts & parts, std::mt19937 & gen) {
	std::uniform_int_distribution<int32_t> scoreDist(-1000, 1000);
	const std::string ptrs("UDLB\0x", 6);
	std::uniform_int_distribution<size_t> ptrDist(0, ptrs.size() - 1);
	for (auto & row : parts.ScoreMatrix_) {
		for (auto & cell : row) {
			cell.upInherit = scoreDist(gen);
			cell.leftInherit = scoreDist(gen);
			cell.diagInherit = scoreDist(gen);
			cell.upInheritPtr = ptrs[ptrDist(gen)];
			cell.leftInheritPtr = ptrs[ptrDist(gen)];
			cell.diagInheritPtr = ptrs[ptrDist(gen)];
		}
	}
}

//empty string if the alignments match, otherwise a description of the first difference
std::string compareNeedle(const std::string & a, const std::string & b, const gapScoringParameters & gaps,
		const substituteMatrix & scoring, std::mt19937 & gen) {
	const uint64_t maxSize = std::max(a.size(), b.size());
	alnParts expected(maxSize, gaps, scoring);
	alnParts actual(maxSize, gaps, scoring);
	fillWithJunk(expected, gen);
	fillWithJunk(actual, gen);
	referenceRunNeedleSave(a, b, expected);
	alignCalc::runNeedleSave(a, b, actual);
	std::stringstream ss;
	ss << "a: " << a << " b: " << b << " gaps: " << gaps.getIdentifer() << " ";
	if (expected.score_ != actual.score_ || expected.gHolder_.score_ != actual.gHolder_.score_) {
		ss << "score " << expected.score_ << " vs " << actual.score_;
		return ss.str();
	}
	if (!(expected.gHolder_.gapInfos_ == actual.gHolder_.gapInfos_)) {
		ss << "gapInfos differ";
		return ss.str();
	}
	for (uint32_t i = 0; i <= a.size(); ++i) {
		for (uint32_t j = 0; j <= b.size(); ++j) {
			const auto & e = expected.ScoreMatrix_[i][j];
			const auto & o = actual.ScoreMatrix_[i][j];
			//row 0 and column 0 only have some fields set, so compare those that are
			if (0 == i && 0 == j) {
				continue;
			}
			bool same = true;
			if (0 == j) {
				same = e.upInherit == o.upInherit && e.upInheritPtr == o.upInheritPtr;
			} else if (0 == i) {
				same = e.leftInherit == o.leftInherit && e.leftInheritPtr == o.leftInheritPtr;
			} else {
				same = e.upInherit == o.upInherit && e.leftInherit == o.leftInherit && e.diagInherit == o.diagInherit &&
						e.upInheritPtr == o.upInheritPtr && e.leftInheritPtr == o.leftInheritPtr && e.diagInheritPtr == o.diagInheritPtr;
			}
			if (!same) {
				ss << "cell " << i << "," << j << " differs";
				return ss.str();
			}
		}
	}
	return "";
}

std::string randomSeq(std::mt19937 & gen, const std::string & alphabet, uint32_t len) {
	std::uniform_int_distribution<size_t> dist(0, alphabet.size() - 1);
	std::string ret;
	for (uint32_t pos = 0; pos < len; ++pos) {
		ret.push_back(alphabet[dist(gen)]);
	}
	return ret;
}

gapScoringParameters randomGaps(std::mt19937 & gen) {
	std::uniform_int_distribution<int32_t> dist(-2, 12);
	return {dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen)};
}

std::vector<substituteMatrix> testScorings(std::mt19937 & gen) {
	std::array<std::array<int32_t, 127>, 127> randomMat{};
	std::uniform_int_distribution<int32_t> dist(-6, 6);
	for (auto & row : randomMat) {
		for (auto & score : row) {
			score = dist(gen);
		}
	}
	return {substituteMatrix(2, -2),
					substituteMatrix::createDegenScoreMatrixNoNInRef(2, -2),
					substituteMatrix::createDegenScoreMatrix(1, -3),
					substituteMatrix::createDegenScoreMatrixCaseInsensitive(2, -2),
					substituteMatrix(randomMat)};
}

const std::vector<std::string> alphabets{"ACGT", "AC", "ACGTN", "ACGTRYKMSWBDHVN", "acgtACGT", "ACGT-"};

}  // namespace

TEST_CASE("needleMaximum matches the reference for every tie and order", "[alignment][alignCalc]") {
	uint32_t failures = 0;
	for (int32_t u = -6; u <= 6; ++u) {
		for (int32_t l = -6; l <= 6; ++l) {
			for (int32_t d = -6; d <= 6; ++d) {
				char expectedPtr = '\0';
				char actualPtr = '\0';
				auto expected = referenceNeedleMaximum(u, l, d, expectedPtr);
				auto actual = alignCalc::needleMaximum(u, l, d, actualPtr);
				failures += (expected != actual || expectedPtr != actualPtr);
			}
		}
	}
	REQUIRE(0 == failures);
	//and over a wide range of values
	std::mt19937 gen(17);
	std::uniform_int_distribution<int32_t> wide(-100000, 100000);
	std::uniform_int_distribution<int32_t> narrow(-3, 3);
	for (uint32_t iter = 0; iter < 200000; ++iter) {
		const int32_t base = wide(gen);
		//values near each other so ties stay common
		const int32_t u = base + narrow(gen);
		const int32_t l = base + narrow(gen);
		const int32_t d = iter % 2 ? base + narrow(gen) : wide(gen);
		char expectedPtr = '\0';
		char actualPtr = '\0';
		auto expected = referenceNeedleMaximum(u, l, d, expectedPtr);
		auto actual = alignCalc::needleMaximum(u, l, d, actualPtr);
		failures += (expected != actual || expectedPtr != actualPtr);
	}
	REQUIRE(0 == failures);
}

TEST_CASE("runNeedleSave matches the reference on every pair of short sequences", "[alignment][alignCalc]") {
	std::mt19937 gen(11);
	auto scorings = testScorings(gen);
	std::vector<gapScoringParameters> gapSets{gapScoringParameters(5, 1), gapScoringParameters(7, 1, 0, 0, 0, 0, 0, 0, 0, 0),
																						gapScoringParameters(3, 1, 0, 0, 3, 1, 0, 0, 3, 1), randomGaps(gen), randomGaps(gen)};
	//every string of length 2 to 5 over a two letter alphabet, against every other (a single base can make runNeedleSave trace back into row 0 and never finish)
	std::vector<std::string> seqs;
	for (uint32_t len = 2; len <= 5; ++len) {
		for (uint32_t bits = 0; bits < (1U << len); ++bits) {
			std::string seq;
			for (uint32_t pos = 0; pos < len; ++pos) {
				seq.push_back((bits >> pos) & 1U ? 'A' : 'C');
			}
			seqs.emplace_back(seq);
		}
	}
	std::string firstFailure;
	uint32_t failures = 0;
	for (const auto & scoring : scorings) {
		for (const auto & gaps : gapSets) {
			for (const auto & a : seqs) {
				for (const auto & b : seqs) {
					auto diff = compareNeedle(a, b, gaps, scoring, gen);
					if (!diff.empty()) {
						++failures;
						if (firstFailure.empty()) {
							firstFailure = diff;
						}
					}
				}
			}
		}
	}
	INFO(firstFailure);
	REQUIRE(0 == failures);
}

TEST_CASE("runNeedleSave matches the reference on random sequences, scorings and gap penalties", "[alignment][alignCalc]") {
	std::mt19937 gen(3);
	auto scorings = testScorings(gen);
	std::uniform_int_distribution<uint32_t> shortLen(2, 12);
	std::uniform_int_distribution<uint32_t> longLen(2, 120);
	std::uniform_int_distribution<size_t> alphabetDist(0, alphabets.size() - 1);
	std::uniform_int_distribution<size_t> scoringDist(0, scorings.size() - 1);
	std::uniform_real_distribution<double> unif(0, 1);
	std::string firstFailure;
	uint32_t failures = 0;
	for (uint32_t iter = 0; iter < 20000; ++iter) {
		const auto & alphabet = alphabets[alphabetDist(gen)];
		auto & lenDist = unif(gen) < 0.5 ? shortLen : longLen;
		auto a = randomSeq(gen, alphabet, lenDist(gen));
		std::string b;
		if (unif(gen) < 0.3) {
			//a mutated copy, so there are real alignments with gaps
			b = a;
			std::uniform_int_distribution<size_t> posDist(0, b.size() - 1);
			for (uint32_t mut = 0; mut < 1 + b.size() / 10; ++mut) {
				auto pos = posDist(gen);
				auto kind = unif(gen);
				if (kind < 0.4) {
					b[pos] = alphabet[std::uniform_int_distribution<size_t>(0, alphabet.size() - 1)(gen)];
				} else if (kind < 0.7 && b.size() > 2) {
					b.erase(pos, 1);
					posDist = std::uniform_int_distribution<size_t>(0, b.size() - 1);
				} else {
					b.insert(pos, 1, alphabet[std::uniform_int_distribution<size_t>(0, alphabet.size() - 1)(gen)]);
					posDist = std::uniform_int_distribution<size_t>(0, b.size() - 1);
				}
			}
		} else {
			b = randomSeq(gen, alphabet, lenDist(gen));
		}
		auto gaps = unif(gen) < 0.7 ? randomGaps(gen) : gapScoringParameters(5, 1, 0, 0, 0, 0, 0, 0, 0, 0);
		auto diff = compareNeedle(a, b, gaps, scorings[scoringDist(gen)], gen);
		if (!diff.empty()) {
			++failures;
			if (firstFailure.empty()) {
				firstFailure = diff;
			}
		}
	}
	INFO(firstFailure);
	REQUIRE(0 == failures);
}
