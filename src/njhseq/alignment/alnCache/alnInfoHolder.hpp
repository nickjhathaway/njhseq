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

//
//  alnInfoHolder.hpp
//
//  Created by Nicholas Hathaway on 1/13/14.
//



#include "njhseq/alignment/alnCache/alnInfoHolderBase.hpp"


#if __APPLE__ == 1 && __cpp_lib_shared_timed_mutex < 201402L && __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ <= 101106
#include <sharedMutex.h>
#else
#include <shared_mutex>
#endif

namespace njhseq {


class alnInfoMasterHolder {

 public:
  // constructors
	alnInfoMasterHolder();
  alnInfoMasterHolder(const gapScoringParameters & gapPars,
  	  const substituteMatrix & scoringArray);
	alnInfoMasterHolder(const std::string &masterDirName, const gapScoringParameters & gapPars,
  	  const substituteMatrix & scoringArray, bool verbose = false);
  // members
  std::unordered_map<std::string, alnInfoHolderBase<alnInfoLocal>> localHolder_;
  std::unordered_map<std::string, alnInfoHolderBase<alnInfoGlobal>> globalHolder_;
  std::hash<std::string> strH_;

  void clearHolders();

  void addHolder(const gapScoringParameters & gapPars,
  	  const substituteMatrix & scoringArray);

  // reading
  void read(const std::string &masterDirName, bool verbose = false);

  // writing
  void write(const std::string &masterDirName, bool verbose = false);

  void mergeOtherHolder(const alnInfoMasterHolder & otherHolder);
};

namespace alignment {

static std::mutex alnCacheDirSearchLock;
static std::unordered_map<std::string, std::unique_ptr<std::shared_timed_mutex>> alnCacheDirLocks;

}  // namespace alignment


}  // namespace njhseq


