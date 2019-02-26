#pragma once

/*
 * OutputStreamWrap.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: nicholashathaway
 */

#include "njhseq/IO/OutputStream.hpp"


namespace njhseq {

class OutputStreamWrap {
public:

	OutputStreamWrap(const OutOptions & opts);

	OutOptions opts_;
	std::unique_ptr<OutputStream> out_;
	std::mutex mut_;
	bool outOpen_ { false };

	void writeNoCheck(const std::string & line, bool flush = false);
	void write(const std::string & line, bool flush = false);

	bool outOpen() const;
	void openOut();

	void closeOut();
	void closeOutForReopening();

};


}  // namespace njhseq





