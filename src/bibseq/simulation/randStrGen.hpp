#pragma once
/*
 * randStrGen.hpp
 *
 *  Created on: Jul 27, 2014
 *      Author: nickhathaway
 */


#include <bibcpp/simulation/randObjGen.hpp>
#include "bibseq/simulation/simulationCommon.hpp"

namespace bibseq {

class randStrGen {
public:

	//constructor
	randStrGen(randomGenerator rGen,
			const std::vector<char> & letters):
				charGen_(bib::randObjectGen<char,uint32_t>(letters)),
				rGen_(rGen){}

	randStrGen(randomGenerator rGen,
			const std::vector<char> & letters,
			const std::vector<uint32_t> & counts):
					charGen_(bib::randObjectGen<char, uint32_t>(letters,counts)),
					rGen_(rGen){}

private:
	//members
	bib::randObjectGen<char, uint32_t> charGen_;
	randomGenerator rGen_;
public:
	//functions
	std::string rStr(uint64_t size);
	VecStr rStrs(uint64_t size, uint32_t num);
	std::string rStr(uint64_t minSize, uint64_t maxSize);
	VecStr rStrs(uint64_t minSize, uint64_t maxSize, uint32_t num);
};

} /* namespace bib */


