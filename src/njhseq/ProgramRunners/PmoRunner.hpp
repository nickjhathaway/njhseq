#pragma once
//
//  PmoRunner.hpp
//
//  Programs for working with whole PMO files, add as a ring to give a program these operations
//

#include "njhseq/programUtils/seqSetUp.hpp"
#include <njhcpp/progutils.h>

namespace njhseq {

class PmoRunner : public njh::progutils::ProgramRunner {
 public:
  PmoRunner();

  static int mergeMinimumPmos(const njh::progutils::CmdArgs & inputCommands);
  static int addSpecimenMetaToPmo(const njh::progutils::CmdArgs & inputCommands);
  static int validatePmo(const njh::progutils::CmdArgs & inputCommands);
};

}  // namespace njhseq
