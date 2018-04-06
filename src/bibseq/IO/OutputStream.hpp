#pragma once
/*
 * OutputFile.hpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/IO/IOOptions/OutOptions.hpp"
#include <bibcpp/files.h>

namespace bibseq {

class OutputStream : public std::ostream{
public:
	OutputStream(const OutOptions & outOpts);
	const OutOptions outOpts_;
	std::unique_ptr<bib::GZSTREAM::ogzstream> outFileGz_;
	std::unique_ptr<std::ofstream> outFile_;

	~OutputStream();
};




} /* namespace bibseq */


