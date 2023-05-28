#pragma once
//
//  seqSetUp.hpp
//
//  Created by Nicholas Hathaway on 10/13/13.
//
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
#include "njhseq/programUtils/SeqSetUpPars.hpp"
#include <njhcpp/progutils.h>
#include "njhseq/objects/collapseObjects/opts/CollapseIterations.hpp"

namespace njhseq {

class seqSetUp : public njh::progutils::ProgramSetUp {
 public:
	using njh::progutils::ProgramSetUp::ProgramSetUp;

  SeqSetUpPars pars_;
  const static VecStr readInFormatsAvailable_;
	const static VecStr pairedReadInFormatsAvailable_;
	const static VecStr singleInFormatsAvailable_;

	void processQualityFiltering();
  bool processDefaultReader(bool readInNamesRequired = true);
  bool processDefaultReader(const VecStr & formats, bool readInNamesRequired = true);
  bool processReadInNames(bool required = true);
  bool processReadInNames(const VecStr & formats, bool required = true);

	bool processJustReadInNames(SeqIOOptions & opts, const VecStr & formats, bool required = true);


	void processGap();
  void processGapRef();
  void processQualThres();

  CollapseIterations processIteratorMap(const bfs::path& parametersFile);
  CollapseIterations processIteratorMapOnPerId(const bfs::path& parametersFile);

  void processKmerLenOptions();
  void processKmerProfilingOptions();
  void processScoringPars();
  void processAlignerDefualts();
  void processDirectoryOutputName(const std::string& defaultName,
                                  bool mustMakeDirectory);
  void processDirectoryOutputName(bool mustMakeDirectory);
  void processWritingOptions();
  void processWritingOptions(OutOptions & opts);
  bool processRefFilename(bool required = false);
  bool processSeq(bool required = false);
  bool processSeq(std::string& inputSeq, const std::string& flag,
                  const std::string& parName, bool required = false,
								 const std::string & flagGrouping = "Misc");
  bool processSeq(seqInfo& inputSeq, const std::string& flag,
                  const std::string& parName, bool required = false,
								 const std::string & flagGrouping = "Misc");
  bool processVerbose();
  bool processDebug();
  bool processQuiet();
  void processAlnInfoInput();

  void processClusteringOptions();
  void processSkipOnNucComp();
  void processAdjustHRuns();

  /**@brief the stub will be place before the flags
   *
   * @param comp the comparison object to set the errors for
   * @param stub the stub to place before the setting flags
   */
  void processComparison(comparison & comp, std::string stub = "");

  // usage prints
//  void printAdditionaInputUsage(std::ostream& out,
//                                const std::string& lowerRemove);
//
//  void printGapUsage(std::ostream & out)const;
//
//  void printKmerProfilingUsage(std::ostream& out);
//  void printQualThresUsage(std::ostream& out);
//  void printAlignmentUsage(std::ostream& out);
//  void printReferenceComparisonUsage(std::ostream& out);
//  void printFileWritingUsage(std::ostream& out, bool all);
//  void printAlnInfoDirUsage(std::ostream& out);
//  void printAdditionalClusteringUsage(std::ostream& out);

};
}  // namespace njhseq


