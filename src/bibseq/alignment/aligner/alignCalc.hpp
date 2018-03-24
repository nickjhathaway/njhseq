#pragma once
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
 * alnCalc.hpp
 *
 *  Created on: Nov 29, 2015
 *      Author: nick
 */

#include "bibseq/alignment/aligner/alnParts.hpp"


namespace bibseq {


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
                            alnParts& parts);
  static void runNeedleOnlyEndGapsSave(const std::string& objA, const std::string& objB,
                            alnParts& parts);

  static void runSmithSave(const std::string& objA, const std::string& objB,
                           alnParts& parts);





  struct MatCursor {
		MatCursor(int32_t icursor, int32_t jcursor, uint32_t lena, uint32_t lenb) :
				icursor_(icursor), jcursor_(jcursor), lena_(lena), lenb_(lenb) {

		}
		int32_t icursor_;
		int32_t jcursor_;
		uint32_t lena_;
		uint32_t lenb_;
	};

  static MatCursor runNeedleDiagonalSaveInit(
  									const std::string& objA,
  									const std::string& objB,
  									uint32_t alignmentBlockSize,
                     alnParts& parts) noexcept ;

  static MatCursor runNeedleDiagonalSaveStep(
  									const std::string& objA,
  									const std::string& objB,
  									uint32_t alignmentBlockSize,
  									uint32_t alignmentBlockWalkbackSize,
  									const MatCursor & lastCursors,
                     alnParts& parts) noexcept ;
  static void runNeedleDiagonalSave(
  									const std::string& objA,
  									const std::string& objB,
  									uint32_t alignmentBlockSize,
  									uint32_t alignmentBlockWalkbackSize,
                     alnParts& parts);

  template <typename T>
  static void rearrangeGlobal(T& readA, T& readB, const typename T::value_type fill,
                              const alnInfoGlobal& gHolder) {
    for (const auto& g : gHolder.gapInfos_) {
      if (g.gapInA_) {
        readA.insert(readA.begin() + g.pos_, g.size_, fill);
      } else {
        readB.insert(readB.begin() + g.pos_, g.size_, fill);
      }
    }
  }

  template <typename T>
  static void rearrangeGlobalQueryOnly(T& readB, const typename T::value_type fill,
                              const alnInfoGlobal& gHolder) {
    for (const auto& g : gHolder.gapInfos_) {
      if (!g.gapInA_) {
      	readB.insert(readB.begin() + g.pos_, g.size_, fill);
      }
    }
  }

  template <typename T>
  static void rearrangeLocal(T& readA, T& readB, const typename T::value_type fill,
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


}  // namespace bibseq

