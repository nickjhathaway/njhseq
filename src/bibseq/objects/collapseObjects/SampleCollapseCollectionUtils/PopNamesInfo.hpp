#pragma once
/*
 * PopNamesInfo.hpp
 *
 *  Created on: Jul 31, 2016
 *      Author: nick
 */

#include "bibseq/utils.h"

namespace bibseq {

class PopNamesInfo {

public:
	PopNamesInfo(std::string populationName, std::set<std::string> samples);
	PopNamesInfo(std::string populationName, VecStr samples);

	std::string populationName_;
	std::set<std::string> samples_;

	bool hasSample(const std::string & sample)const;
private:
	void checkPopNameThrow()const;
public:

};


}  // namespace bibseq






