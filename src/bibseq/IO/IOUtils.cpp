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
#include "bibseq/IO/OutputStream.hpp"
#include "bibseq/IO/InputStream.hpp"



namespace bibseq {

VecStr getInputValues(const std::string & valuesStr, const std::string & delim){
	VecStr ret;
	if("" == valuesStr){
		return ret;
	}
	if (bfs::path(valuesStr).filename().string().length() <= 255 && bfs::exists(valuesStr)) {
		InputStream infile{bfs::path(valuesStr)};
		std::string line = "";
		while(bib::files::crossPlatGetline(infile, line)){
			ret.emplace_back(line);
		}
	}else{
		ret = tokenizeString(valuesStr, delim);
	}
	return ret;
}


void gzZipFile(const IoOptions & opts){
	if (opts.out_.outExists() && !opts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, file " << opts.out_.outName()
				<< " already exists, use --overWrite to overwrite it" << "\n";
		throw std::runtime_error { ss.str() };
	}
//	bib::GZSTREAM::ogzstream outStream(opts.out_.outName());
//	outStream << bib::files::get_file_contents(opts.in_.inFilename_, false);
	//read in chunks so that the entire file doesn't have to be read in if it's very large
	/**@todo find an apprioprate chunkSize or */
	uint32_t chunkSize = 4096 * 10;
	bib::GZSTREAM::ogzstream outstream;
	opts.out_.openGzFile(outstream);
	std::ifstream infile(opts.in_.inFilename_.string(), std::ios::binary);
	std::vector<char> buffer(chunkSize);
	infile.read(buffer.data(), sizeof(char) * chunkSize);
	std::streamsize bytes = infile.gcount();

	while(bytes > 0){
		outstream.write(buffer.data(), bytes * sizeof(char));
		infile.read(buffer.data(), sizeof(char) * chunkSize);
		bytes = infile.gcount();
	}
}




}  // namespace bibseq
