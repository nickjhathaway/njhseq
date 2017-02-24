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
#include "IOUtils.hpp"

namespace bibseq {


void gzZipFile(const IoOptions & opts){
	if (opts.out_.outExists() && !opts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, file " << opts.out_.outName()
				<< " already exists, use --overWrite to overwrite it" << "\n";
		throw std::runtime_error { ss.str() };
	}
	bib::GZSTREAM::ogzstream outStream(opts.out_.outName());
	outStream << bib::files::get_file_contents(opts.in_.inFilename_, false);
}




}  // namespace bibseq
