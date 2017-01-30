#pragma once
//
//  vectorUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 7/18/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//
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

#include "bibseq/utils/stringUtils.hpp"
#include "bibseq/utils/vectorCalculations.hpp"
#include "bibseq/utils/numUtils.hpp"

namespace bibseq {
template<typename T>
void addAsStrToVec(VecStr& vec, const T& e) {
	vec.emplace_back(estd::to_string(e));
}

template<typename T>
void addAsStrToVec(VecStr& vec, const std::vector<T>& items) {
	for (const auto& e : items) {
		vec.emplace_back(estd::to_string(e));
	}
}

template<typename T>
void toVecStrAdd(VecStr & vec, const T& last) {
	addAsStrToVec(vec, last);
}

template<typename N, typename ... T>
void toVecStrAdd(VecStr& vec, const N& next, const T&... rest) {
	addAsStrToVec(vec, next);
	toVecStrAdd(vec, rest...);
}

template<typename ... T>
VecStr toVecStr(const T&... items) {
	VecStr ret;
	ret.reserve(sizeof...(items));
	toVecStrAdd(ret, items...);
	return ret;
}

template<typename T>
std::vector<T> vecStrToVecNum(const VecStr& vec) {
	std::vector<T> ans;
	for (const auto& s : vec) {
		ans.emplace_back(static_cast<T>(std::stof(s)));
	}
	return ans;
}

template<typename T>
void prependVec(std::vector<T>& vec, const std::vector<T>& otherVec) {
	vec.insert(vec.begin(), otherVec.begin(), otherVec.end());
}
template<typename T>
void prependVec(std::vector<T>& vec, const T& singleValue) {
	vec.insert(vec.begin(), singleValue);
}

template<typename T>
void printVector(const std::vector<T>& vec, const std::string& delim = " ",
		std::ostream& out = std::cout) {
	out << vectorToString(vec, delim) << "\n";

}

template<class T>
void addOtherVec(std::vector<T>& reads, const std::vector<T>& otherVec) {
	reads.reserve(reads.size() + otherVec.size());
	reads.insert(reads.end(), otherVec.begin(), otherVec.end());
}

template<typename T>
void increaseVectorByOne(std::vector<T>& vec) {
	std::for_each(vec.begin(), vec.end(), [](T& num) {++num;});
}
template<typename T>
T productOfVecElements(const std::vector<T>& vec) {
	T product = 1;
	for (const auto& iter : vec) {
		product *= (iter);
	}
	return product;
}

template<typename T>
std::vector<T> repeatVector(const std::vector<T>& vec,
		const std::vector<uint64_t>& repeatNumber) {
	std::vector<T> ans;
	if (repeatNumber.size() == 1) {
		ans.reserve(repeatNumber[0] * vec.size());
		for (uint64_t i = 0; i < repeatNumber[0]; ++i) {
			addOtherVec(ans, vec);
		}
	} else if (repeatNumber.size() == vec.size()) {
		uint32_t pos = 0;
		auto sum = vectorSum(repeatNumber);
		ans.reserve(sum);
		for (const auto& iter : vec) {
			addOtherVec(ans, std::vector<T>(repeatNumber[pos], iter));
			pos++;
		}
	} else {
		std::cout << "Repeat number vector needs to be either same size of the "
				"vector to be repeated or needs to be one number" << "\n";
	}
	return ans;
}

template<typename T>
std::vector<std::vector<T>> permuteVector(const std::vector<T>& vec,
		size_t numberOf) {

	std::vector<std::vector<T>> ans;
	std::vector<uint64_t> overAllRepeatVec;
	overAllRepeatVec.push_back(pow(vec.size(), numberOf));
	std::vector<uint64_t> repeatFactorVec;
	repeatFactorVec.push_back(1);
	std::vector<uint64_t> sizeOfVectorVec;
	sizeOfVectorVec.push_back(vec.size());
	for (size_t i = 0; i < numberOf; ++i) {
		overAllRepeatVec[0] = overAllRepeatVec[0] / vec.size();
		auto numberOfCurrentRepeats = repeatVector(repeatFactorVec,
				sizeOfVectorVec);
		auto subRepeat = repeatVector(vec, numberOfCurrentRepeats);
		auto bigRepeat = repeatVector(subRepeat, overAllRepeatVec);
		ans.push_back(bigRepeat);
		repeatFactorVec[0] = repeatFactorVec[0] * sizeOfVectorVec[0];
	}
	//std::cout << "RepeatThenNumbers:" << getDuration(startTime1) << "\n";
	std::reverse(ans.begin(), ans.end());
	size_t rowStop = ans.size();
	size_t colStop = ans[0].size();
	std::vector<std::vector<T>> realAns;
	for (size_t i = 0; i < colStop; ++i) {
		std::vector<T> tempVec;
		for (size_t j = 0; j < rowStop; ++j) {
			tempVec.push_back(ans[j][i]);
		}
		realAns.push_back(tempVec);
	}
	//std::cout << "Reorganize time: " << getDuration(startTime2) << "\n";
	return realAns;
}

template<typename T>
std::vector<std::vector<T>> fastPermuteVectorOneLength(std::vector<T> vec) {
	std::vector<std::vector<T>> ans;
	int numOfPermutes = Factorial((int) vec.size());
	ans.reserve(numOfPermutes);
	do {
		ans.push_back(vec);
	} while (std::next_permutation(vec.begin(), vec.end()));
	return ans;
}
const VecStr fastPermuteVectorOneLength(std::string vec);

template<typename T>
std::vector<std::vector<T>> findUniqueVectors(
		const std::vector<std::vector<T>>& vec) {
	std::vector<std::vector<T>> ans;

	for (const auto& iter : vec) {
		bool unique = true;
		for (auto i = 0; i < iter.size(); ++i) {
			for (auto j = i + 1; j < iter.size(); ++j) {
				if (iter[i] == iter[j]) {
					unique = false;
					break;
				}
			}
			if (!unique) {
				break;
			}
		}
		if (unique) {
			ans.push_back(iter);
		}
	}
	return ans;
}

std::vector<std::vector<int>> permuteLengthN(int N);

template<typename T>
void eraseTrailingZeros(std::vector<T>& vec) {
	while (vec.size() > 0 && vec[vec.size() - 1] == 0) {
		vec.erase(vec.end() - 1);
	}
	return;
}

template<typename T>
bool vectorContains(const std::vector<T>& vec, const T& search) {
	for (const auto& iter : vec) {
		if (iter == search) {
			return true;
		}
	}
	return false;
}
/*
 template <typename T>
 std::vector<T> getUnique(const std::vector<T>& vec) {
 std::vector<T> ans;
 for (const auto& element : vec) {
 if (!vectorContains(ans, element)) {
 ans.push_back(element);
 }
 }
 return ans;
 }*/
/*
 template<typename T>
 void removeElements;
 */
template<class T>
std::vector<std::vector<T>> collapseUniqueVectors(
		const std::vector<std::vector<T>>& vec) {
	std::vector<std::vector<T>> ans;
	int count = 0;
	for (const auto& iter : vec) {
		if (count == 0) {
			ans.push_back(iter);
		} else {
			if (!vectorContains(ans, iter)) {
				ans.push_back(iter);
			}
		}
		++count;
	}
	return ans;
}
template<typename T>
void removeDuplicates(std::vector<T>& vec) {
	std::sort(vec.begin(), vec.end());
	auto last = std::unique(vec.begin(), vec.end());
	vec.erase(last, vec.end());
	return;
}

template<typename T>
void outputVectorOfVectors(const std::vector<std::vector<T>>& vec,
		const std::string& delim, std::ostream& out) {
	for (const auto& iter : vec) {
		printVector(iter, delim, out);
	}
}

//catenateVectors
template<class T>
const std::vector<T> concatVecs(const std::vector<T>& vec1,
		const std::vector<T>& vec2) {
	std::vector<T> ans;
	ans.reserve(vec1.size() + vec2.size());
	ans.insert(ans.end(), vec1.begin(), vec1.end());
	ans.insert(ans.end(), vec2.begin(), vec2.end());
	return ans;
}

// getting positions of targets in vectors
template<class T>
std::vector<uint32_t> getPositionsOfTarget(const std::vector<T>& vec,
		const T& target) {
	uint32_t pos = 0;
	std::vector<uint32_t> positions;
	for (const auto& iter : vec) {
		if (iter == target) {
			positions.push_back(pos);
		}
		++pos;
	}
	return positions;
}
template<class T>
uint32_t getFirstPositionOfTarget(const std::vector<T>& vec, const T& target) {
	uint32_t pos = 0;
	for (const auto& iter : vec) {
		if (iter == target) {
			return pos;
		}
		++pos;
	}
	return UINT32_MAX;
}

std::vector<uint32_t> getPositionsOfTarget(const VecStr& vec,
		const std::string& target);

std::vector<uint32_t> getPositionsOfSubStrTarget(const VecStr& vec,
		const std::string& target);

std::vector<uint32_t> getPositionsOfTargetStartsWith(const VecStr& vec,
		const std::string& target);

template<class T>
std::vector<T> getTargetsAtPositions(const std::vector<T>& vec,
		const std::vector<uint32_t>& positions) {
	std::vector<T> ans;
	for (const auto& pIter : positions) {
		if (pIter < vec.size()) {
			ans.push_back(vec[pIter]);
		}
	}
	return ans;
}
template<class T>
std::vector<T> getTargetsNotAtPositions(const std::vector<T>& vec,
		const std::vector<uint32_t>& positions) {
	std::vector<T> ans;
	uint32_t currentPos = 0;
	for (const auto& element : vec) {
		if (!vectorContains(positions, currentPos)) {
			ans.push_back(element);
		}
		++currentPos;
	}
	return ans;
}

template<class T>
std::vector<uint32_t> getPositionsMultipleTargets(const std::vector<T>& vec,
		const std::vector<T>& targets) {
	std::vector<uint32_t> positions;
	for (const auto& iter : targets) {
		addOtherVec(positions, getPositionsOfTarget(vec, iter));
	}
	return positions;
}

VecStr getUniqueStrings(const VecStr& vec);

VecStr getStringsContains(const VecStr& vec, const std::string& contains);

template<typename T>
void removeElement(std::vector<T>& vec, const T& element) {
	vec.erase(std::remove(vec.begin(), vec.end(), element), vec.end());
}
template<typename T>
void removeElements(std::vector<T>& vec, const std::vector<T>& elements) {
	for (const auto& element : elements) {
		removeElement(vec, element);
	}
}

template<typename T>
std::vector<std::vector<T>> findClumps(const std::vector<T>& positions,
		int windowSize, int appearsAtLeast, int extraLength = 0,
		bool findAllClumps = false) {
	std::vector<std::vector<T>> ans;
	for (auto i : iter::range(positions.size() - appearsAtLeast + 1)) {
		std::vector<T> currentClump = { positions[i] };
		for (auto j : iter::range(i + appearsAtLeast - 1, positions.size())) {
			if (positions[j] - positions[i] + extraLength <= windowSize) {
				if (j == i + appearsAtLeast - 1) {
					for (auto k : iter::range(1, appearsAtLeast)) {
						currentClump.push_back(positions[i + k]);
					}
				} else {
					currentClump.push_back(positions[j]);
				}
			} else {
				break;
			}
		}
		if (currentClump.size() >= appearsAtLeast) {
			ans.push_back(currentClump);
			if (!findAllClumps) {
				return ans;
			}
		}
	}
	return ans;
}
template<typename T>
std::vector<T> getSubVector(const std::vector<T>& vec, uint32_t start,
		uint32_t size) {
	std::vector<T> ans(vec.begin() + start, vec.begin() + size + start);
	return ans;
}
template<typename T>
std::vector<T> getSubVector(const std::vector<T>& vec, uint32_t start) {
	std::vector<T> ans(vec.begin() + start, vec.end());
	return ans;
}

std::string getSubVector(const std::string& vec, uint32_t start, uint32_t size);

std::string getSubVector(const std::string& vec, uint32_t start);

template<typename T>
std::vector<T> convertStringToVector(const std::string& str,
		const std::string& delim) {
	std::vector<T> ans;
	VecStr toks = tokenizeString(str, delim);
	for (const auto& t : toks) {
		std::stringstream tempStream;
		tempStream << t;
		T tempObject;
		tempStream >> tempObject;
		ans.push_back(tempObject);
	}
	return ans;
}

template<typename T>
VecStr numVecToVecStr(const std::vector<T>& nums) {
	VecStr ans;
	ans.reserve(nums.size());
	// std::generate_n(ans.begin(), nums.size(), [](const T & num) {return
	// estd::to_string(num);});
	for (const auto& num : nums) {
		ans.emplace_back(estd::to_string(num));
	}
	return ans;
}

template<typename T>
std::vector<T> getRange(const T& start, const T& stop) {
	std::vector<T> ans(stop - start + 1);
	bib::iota(ans, start);
	return ans;
}

double getMeanFromVecStr(const VecStr & strNums);

template<typename T>
T getSumOfVecStr(const VecStr & vec) {
	return vectorSum(bib::lexical_cast_con<VecStr, std::vector<T>>(vec));
}

}  // namespace bibseq

