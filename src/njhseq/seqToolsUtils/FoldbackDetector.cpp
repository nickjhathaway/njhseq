/*
 * FoldbackDetector.cpp
 *
 *  Created on: Nov 5, 2026
 *      Author: nicholashathaway
 */
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

#include "FoldbackDetector.hpp"

namespace njhseq {

namespace {

constexpr int8_t notABase = -1;

int8_t baseCode(const char base){
	switch(base){
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
		default: return notABase;
	}
}

struct Vote{
	int64_t sum_;
	uint32_t left_;
	uint32_t right_;
};

/**@brief Edit distance between a stretch of the left arm and the matching stretch of the right arm, read back from its end
 *
 */
uint32_t armEditDistance(const std::vector<int8_t> & codes, const uint32_t leftStart, const uint32_t leftLen, const uint32_t rightEnd, const uint32_t rightLen){
	std::vector<uint32_t> previous(rightLen + 1);
	std::vector<uint32_t> current(rightLen + 1);
	std::iota(previous.begin(), previous.end(), 0);
	for(uint32_t pos = 0; pos < leftLen; ++pos){
		current[0] = pos + 1;
		for(uint32_t other = 0; other < rightLen; ++other){
			const bool pairs = notABase != codes[leftStart + pos] && 3 - codes[leftStart + pos] == codes[rightEnd - 1 - other];
			current[other + 1] = std::min({previous[other] + (pairs ? 0 : 1), previous[other + 1] + 1, current[other] + 1});
		}
		std::swap(previous, current);
	}
	return previous[rightLen];
}

/**@brief Extend the pairing of the two arms, with errors allowed, as far as it stays worth it
 *
 * @param codes the read's bases as codes
 * @param left the first position on the left arm to pair
 * @param leftStep which way the left arm goes
 * @param right the first position on the right arm to pair
 * @param rightStep which way the right arm goes
 * @param leftAvailable how many bases the left arm has to use
 * @param rightAvailable how many bases the right arm has to use
 * @param available how many bases can be used between both arms, going inward the arms can't use the same bases
 * @param pars the settings
 * @return how many bases of the left arm and of the right arm the best extension used
 */
std::pair<uint32_t, uint32_t> extendArms(const std::vector<int8_t> & codes, const int64_t left, const int32_t leftStep,
		const int64_t right, const int32_t rightStep, const uint64_t leftAvailable, const uint64_t rightAvailable, const uint64_t available, const FoldbackDetector::Pars & pars){
	constexpr int32_t match = 1;
	//harsh enough that random sequence past the end of an arm doesn't pair, two arms at 10% error each still gain from each base
	constexpr int32_t penalty = -3;
	constexpr int32_t unreachable = std::numeric_limits<int32_t>::min() / 2;
	const auto band = static_cast<int64_t>(pars.extensionBand_);
	//columns are how many more of the right arm than of the left arm have been used, plus the band
	std::vector<int32_t> previous(2 * band + 1, unreachable);
	std::vector<int32_t> current(2 * band + 1, unreachable);
	int32_t best = 0;
	std::pair<uint32_t, uint32_t> bestUsed{0, 0};
	for(int64_t leftUsed = 0; static_cast<uint64_t>(leftUsed) <= leftAvailable; ++leftUsed){
		std::fill(current.begin(), current.end(), unreachable);
		int32_t rowBest = unreachable;
		for(int64_t column = 0; column <= 2 * band; ++column){
			const int64_t rightUsed = leftUsed + column - band;
			if(rightUsed < 0 || static_cast<uint64_t>(rightUsed) > rightAvailable || static_cast<uint64_t>(leftUsed + rightUsed) > available){
				continue;
			}
			int32_t score = 0 == leftUsed && 0 == rightUsed ? 0 : unreachable;
			if(leftUsed > 0 && rightUsed > 0 && unreachable != previous[column]){
				const auto leftCode = codes[left + leftStep * (leftUsed - 1)];
				const bool pairs = notABase != leftCode && 3 - leftCode == codes[right + rightStep * (rightUsed - 1)];
				score = std::max(score, previous[column] + (pairs ? match : penalty));
			}
			//a base used on only the left arm comes from the previous row one column over
			if(leftUsed > 0 && column + 1 <= 2 * band && unreachable != previous[column + 1]){
				score = std::max(score, previous[column + 1] + penalty);
			}
			if(rightUsed > 0 && column > 0 && unreachable != current[column - 1]){
				score = std::max(score, current[column - 1] + penalty);
			}
			current[column] = score;
			rowBest = std::max(rowBest, score);
			if(score > best){
				best = score;
				bestUsed = {static_cast<uint32_t>(leftUsed), static_cast<uint32_t>(rightUsed)};
			}
		}
		if(rowBest < best - pars.extensionXDrop_){
			break;
		}
		std::swap(previous, current);
	}
	return bestUsed;
}

}  // namespace

uint32_t FoldbackDetector::Fold::leftArmLen() const{
	return leftEnd_ - leftStart_;
}

uint32_t FoldbackDetector::Fold::rightArmLen() const{
	return rightEnd_ - rightStart_;
}

uint32_t FoldbackDetector::Fold::shorterArmLen() const{
	return std::min(leftArmLen(), rightArmLen());
}

uint32_t FoldbackDetector::Fold::loopLen() const{
	return rightStart_ - leftEnd_;
}

double FoldbackDetector::Fold::anchorsPerBase() const{
	return 0 == anchoredLen_ ? 0.0 : static_cast<double>(anchors_) / anchoredLen_;
}

bool FoldbackDetector::Fold::hasCandidate() const{
	return anchors_ > 0;
}

FoldbackDetector::FoldbackDetector() = default;

FoldbackDetector::FoldbackDetector(Pars pars) : pars_(pars){
}

FoldbackDetector::Fold FoldbackDetector::detect(const std::string & seq) const{
	Fold ret;
	const uint32_t kmerLen = pars_.kmerLen_;
	if(seq.size() < 2 * static_cast<uint64_t>(kmerLen)){
		return ret;
	}
	std::vector<int8_t> codes(seq.size());
	std::transform(seq.begin(), seq.end(), codes.begin(), baseCode);

	//each kmer and its position
	const uint64_t kmerMask = (uint64_t{1} << (2 * kmerLen)) - 1;
	std::vector<std::pair<uint64_t, uint32_t>> kmers;
	kmers.reserve(seq.size());
	uint64_t kmer = 0;
	uint32_t basesInARow = 0;
	for(uint32_t pos = 0; pos < codes.size(); ++pos){
		if(notABase == codes[pos]){
			basesInARow = 0;
			continue;
		}
		kmer = ((kmer << 2) | static_cast<uint64_t>(codes[pos])) & kmerMask;
		++basesInARow;
		if(basesInARow >= kmerLen){
			kmers.emplace_back(kmer, pos + 1 - kmerLen);
		}
	}
	std::sort(kmers.begin(), kmers.end());

	const auto revCompOf = [kmerLen](uint64_t forward){
		uint64_t revComp = 0;
		for(uint32_t pos = 0; pos < kmerLen; ++pos){
			revComp = (revComp << 2) | (3 - (forward & 3));
			forward >>= 2;
		}
		return revComp;
	};
	const auto lowerBound = [&kmers](const uint64_t lookingFor){
		return std::lower_bound(kmers.begin(), kmers.end(), std::make_pair(lookingFor, uint32_t{0}));
	};

	std::vector<Vote> votes;
	auto groupStart = kmers.begin();
	while(kmers.end() != groupStart){
		const auto groupEnd = std::upper_bound(groupStart, kmers.end(), std::make_pair(groupStart->first, std::numeric_limits<uint32_t>::max()));
		const auto groupRevComp = revCompOf(groupStart->first);
		//each pair is seen from both of its kmers so only take it from the kmer on the left which is checked with the positions below
		if(groupRevComp != groupStart->first && static_cast<uint64_t>(std::distance(groupStart, groupEnd)) <= pars_.maxKmerOccurrences_){
			const auto otherStart = lowerBound(groupRevComp);
			auto otherEnd = otherStart;
			while(kmers.end() != otherEnd && groupRevComp == otherEnd->first){
				++otherEnd;
			}
			if(static_cast<uint64_t>(std::distance(otherStart, otherEnd)) <= pars_.maxKmerOccurrences_){
				for(auto leftKmer = groupStart; leftKmer != groupEnd; ++leftKmer){
					for(auto rightKmer = otherStart; rightKmer != otherEnd; ++rightKmer){
						if(rightKmer->second >= leftKmer->second + kmerLen){
							votes.emplace_back(Vote{static_cast<int64_t>(leftKmer->second) + rightKmer->second, leftKmer->second, rightKmer->second});
						}
					}
				}
			}
		}
		groupStart = groupEnd;
	}
	if(votes.empty()){
		return ret;
	}

	//the strongest range of sums
	std::sort(votes.begin(), votes.end(), [](const Vote & vote1, const Vote & vote2){
		return vote1.sum_ == vote2.sum_ ? vote1.left_ < vote2.left_ : vote1.sum_ < vote2.sum_;
	});
	size_t peakStart = 0;
	size_t peakSize = 0;
	size_t windowStart = 0;
	for(size_t pos = 0; pos < votes.size(); ++pos){
		while(votes[pos].sum_ - votes[windowStart].sum_ > static_cast<int64_t>(pars_.seedBand_)){
			++windowStart;
		}
		if(pos + 1 - windowStart > peakSize){
			peakSize = pos + 1 - windowStart;
			peakStart = windowStart;
		}
	}
	const std::vector<Vote> peak(votes.begin() + static_cast<int64_t>(peakStart), votes.begin() + static_cast<int64_t>(peakStart + peakSize));
	const auto seed = [&peak](){
		std::vector<uint32_t> lefts;
		lefts.reserve(peak.size());
		for(const auto & vote : peak){
			lefts.emplace_back(vote.left_);
		}
		std::nth_element(lefts.begin(), lefts.begin() + static_cast<int64_t>(lefts.size() / 2), lefts.end());
		const auto middle = lefts[lefts.size() / 2];
		return *std::find_if(peak.begin(), peak.end(), [middle](const Vote & vote){
			return vote.left_ == middle;
		});
	}();

	//chain out from the middle of the peak in both directions, the sum is allowed to drift more the further the next anchor is
	//but only so far, anchors that are far apart with this level of error aren't from the same arms
	std::sort(votes.begin(), votes.end(), [](const Vote & vote1, const Vote & vote2){
		return vote1.left_ == vote2.left_ ? vote1.right_ > vote2.right_ : vote1.left_ < vote2.left_;
	});
	const auto allowedDrift = [](const uint32_t distance){
		return static_cast<int64_t>(8 + distance / 10);
	};
	const uint32_t maxChainGap = pars_.maxChainGap_;
	std::vector<Vote> chain{seed};
	{
		Vote last = seed;
		for(auto vote = votes.rbegin(); votes.rend() != vote; ++vote){
			if(vote->left_ < last.left_ && vote->right_ > last.right_ && last.left_ - vote->left_ <= maxChainGap && std::abs(vote->sum_ - last.sum_) <= allowedDrift(last.left_ - vote->left_)){
				chain.emplace_back(*vote);
				last = *vote;
			}
		}
		std::reverse(chain.begin(), chain.end());
		last = seed;
		for(const auto & vote : votes){
			if(vote.left_ > last.left_ && vote.right_ < last.right_ && vote.left_ - last.left_ <= maxChainGap && std::abs(vote.sum_ - last.sum_) <= allowedDrift(vote.left_ - last.left_)){
				chain.emplace_back(vote);
				last = vote;
			}
		}
	}
	const auto & outer = chain.front();
	const auto & inner = chain.back();
	ret.anchors_ = static_cast<uint32_t>(chain.size());
	ret.anchoredLen_ = inner.left_ + kmerLen - outer.left_;
	ret.otherVotes_ = static_cast<uint32_t>(votes.size() - chain.size());

	//how well the arms pair between the anchors, going from the outside in
	uint64_t paired = kmerLen;
	uint64_t edits = 0;
	for(size_t pos = 1; pos < chain.size(); ++pos){
		const auto & previous = chain[pos - 1];
		const auto & next = chain[pos];
		const int64_t leftGap = static_cast<int64_t>(next.left_) - (previous.left_ + kmerLen);
		const int64_t rightGap = static_cast<int64_t>(previous.right_) - (next.right_ + kmerLen);
		if(leftGap >= 0 && rightGap >= 0){
			constexpr int64_t longestToAlign = 300;
			edits += leftGap > longestToAlign || rightGap > longestToAlign ?
					static_cast<uint64_t>(std::max(leftGap, rightGap)) :
					armEditDistance(codes, previous.left_ + kmerLen, static_cast<uint32_t>(leftGap), previous.right_, static_cast<uint32_t>(rightGap));
			paired += kmerLen + static_cast<uint64_t>(std::min(leftGap, rightGap));
		}else{
			//overlapping anchors, what they share has to be taken off and if they don't overlap by the same amount there was an insertion or deletion
			edits += static_cast<uint64_t>(std::abs(leftGap - rightGap));
			paired += static_cast<uint64_t>(std::max<int64_t>(0, kmerLen + std::min(leftGap, rightGap)));
		}
	}
	ret.identity_ = static_cast<double>(paired > edits ? paired - edits : 0) / static_cast<double>(paired);

	//find the ends of the arms
	{
		const uint64_t innerLeft = inner.left_ + kmerLen;
		const uint64_t innerRight = inner.right_;
		const uint64_t between = innerRight > innerLeft ? innerRight - innerLeft : 0;
		const auto usedInner = extendArms(codes, static_cast<int64_t>(innerLeft), 1, static_cast<int64_t>(innerRight) - 1, -1, between, between, between, pars_);
		ret.leftEnd_ = static_cast<uint32_t>(innerLeft + usedInner.first);
		ret.rightStart_ = static_cast<uint32_t>(innerRight - usedInner.second);

		const uint64_t leftAvailable = outer.left_;
		const uint64_t rightAvailable = seq.size() - (outer.right_ + kmerLen);
		const auto usedOuter = extendArms(codes, static_cast<int64_t>(outer.left_) - 1, -1, static_cast<int64_t>(outer.right_ + kmerLen), 1,
				leftAvailable, rightAvailable, leftAvailable + rightAvailable, pars_);
		ret.leftStart_ = outer.left_ - usedOuter.first;
		ret.rightEnd_ = static_cast<uint32_t>(outer.right_ + kmerLen + usedOuter.second);
	}

	const double maxLoop = pars_.maxLoopPerArm_ * std::max(ret.leftArmLen(), ret.rightArmLen()) + pars_.maxLoopFlank_;
	ret.called_ = ret.shorterArmLen() >= pars_.minArmLen_ && ret.anchors_ >= pars_.minAnchors_ && ret.loopLen() <= maxLoop &&
			ret.identity_ >= pars_.minIdentity_ && ret.anchorsPerBase() >= pars_.minAnchorsPerBase_;
	return ret;
}

}  // namespace njhseq
