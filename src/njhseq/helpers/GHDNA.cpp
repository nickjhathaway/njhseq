/*
 * GHDNA.cpp
 *
 *  Created on: Jan 1, 2019
 *      Author: nicholashathaway
 */




#include "GHDNA.hpp"

namespace njhseq {

const std::array<int32_t, 256>  GHDNA::ntValNoN = {
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 3, 5, 7, 5, 5, 5, 11, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 7, 5, 5, 5, 11, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };


uint64_t GHDNA::gen(const std::string & str){
	uint32_t tau = (str.size() - (str.size() % 2));
	uint32_t beta = (tau/2) - 1;
	double R = 0;
	for(uint32_t u = tau; u < str.size() + 1; ++u){
		R += (ntValNoN[str[u - 1]]  - (str.size() - std::sqrt(u)))/ntValNoN[str[u - 1]];
	}
	double C = 0;
	for(uint32_t i = 1; i < beta + 1; ++i){
		uint32_t nuc1 = ntValNoN[str[(2*i -1 - 1)]];
		uint32_t nuc2 = ntValNoN[str[(2*i -0 - 1)]];
		uint32_t nuc3 = ntValNoN[str[(2*i +1 - 1)]];

		C += ((nuc1 - std::sqrt((i%2) + 1) ) * (nuc2 + + std::sqrt((i%3) + 1) ) )/(nuc3);
	}
	C -= R;
	double ID = str.size()/C;
	uint64_t prehash = std::abs(std::round(ID * std::pow(10, 14) )-str.size());
	std::string prehashAsStr = njh::leftPadNumStr<uint64_t>(prehash, std::pow(10, 13));

	std::string final = njh::pasteAsStr(prehashAsStr.substr(7,7), str.size()%10, prehashAsStr.substr(0,6));
	return njh::StrToNumConverter::stoToNum<uint64_t>(final);
}

}  // namespace njhseq
