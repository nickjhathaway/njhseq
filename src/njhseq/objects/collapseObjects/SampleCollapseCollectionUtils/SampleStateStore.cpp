/*
 * SampleStateStore.cpp
 *
 *  Created on: Oct 22, 2026
 *      Author: nick
 */
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
#include "SampleStateStore.hpp"

#include <array>
#include <cstring>
#include <zlib.h>

namespace njhseq {
namespace collapse {

namespace {
//the texts in a fixed order so a state can be serialized as one buffer
std::array<std::optional<std::string> *, 6> stateTexts(SampleState & state) {
	return { &state.final_, &state.finalInputClusters_, &state.excluded_,
			&state.excludedInputClusters_, &state.refCompInfos_, &state.refCompInfosExcluded_ };
}
std::array<const std::optional<std::string> *, 6> stateTexts(const SampleState & state) {
	return { &state.final_, &state.finalInputClusters_, &state.excluded_,
			&state.excludedInputClusters_, &state.refCompInfos_, &state.refCompInfosExcluded_ };
}
}  // namespace

size_t SampleState::sizeInBytes() const {
	size_t ret = 0;
	for (const auto * text : stateTexts(*this)) {
		if (text->has_value()) {
			ret += (*text)->size();
		}
	}
	return ret;
}

SampleStateStore::CompressedState SampleStateStore::compress(
		const SampleState & state) {
	//each text is a presence byte then, if present, its length and its bytes
	std::string raw;
	raw.reserve(state.sizeInBytes() + 6 * (1 + sizeof(uint64_t)));
	for (const auto * text : stateTexts(state)) {
		raw.push_back(text->has_value() ? 1 : 0);
		if (text->has_value()) {
			uint64_t textSize = (*text)->size();
			raw.append(reinterpret_cast<const char *>(&textSize), sizeof(textSize));
			raw.append(**text);
		}
	}
	CompressedState ret;
	ret.uncompressedSize_ = raw.size();
	uLongf compressedSize = compressBound(raw.size());
	ret.data_.resize(compressedSize);
	auto status = compress2(reinterpret_cast<Bytef *>(ret.data_.data()), &compressedSize,
			reinterpret_cast<const Bytef *>(raw.data()), raw.size(), Z_BEST_SPEED);
	if (Z_OK != status) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, zlib compress2 failed with status " << status << "\n";
		throw std::runtime_error { ss.str() };
	}
	ret.data_.resize(compressedSize);
	ret.data_.shrink_to_fit();
	return ret;
}

SampleState SampleStateStore::decompress(const CompressedState & compressed) {
	std::string raw(compressed.uncompressedSize_, '\0');
	uLongf rawSize = raw.size();
	auto status = uncompress(reinterpret_cast<Bytef *>(raw.data()), &rawSize,
			reinterpret_cast<const Bytef *>(compressed.data_.data()), compressed.data_.size());
	if (Z_OK != status || rawSize != raw.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, zlib uncompress failed with status " << status << "\n";
		throw std::runtime_error { ss.str() };
	}
	SampleState state;
	size_t pos = 0;
	for (auto * text : stateTexts(state)) {
		bool present = 1 == raw[pos];
		++pos;
		if (present) {
			uint64_t textSize = 0;
			std::memcpy(&textSize, raw.data() + pos, sizeof(textSize));
			pos += sizeof(textSize);
			*text = raw.substr(pos, textSize);
			pos += textSize;
		}
	}
	return state;
}

void SampleStateStore::set(const std::string & sampleName, SampleState state) {
	auto compressed = compress(state);
	std::lock_guard<std::mutex> lock(mut_);
	states_[sampleName] = std::move(compressed);
}

std::optional<SampleState> SampleStateStore::get(
		const std::string & sampleName) const {
	CompressedState compressed;
	{
		std::lock_guard<std::mutex> lock(mut_);
		auto search = states_.find(sampleName);
		if (states_.end() == search) {
			return std::nullopt;
		}
		compressed = search->second;
	}
	return decompress(compressed);
}

bool SampleStateStore::has(const std::string & sampleName) const {
	std::lock_guard<std::mutex> lock(mut_);
	return states_.end() != states_.find(sampleName);
}

void SampleStateStore::erase(const std::string & sampleName) {
	std::lock_guard<std::mutex> lock(mut_);
	states_.erase(sampleName);
}

VecStr SampleStateStore::samples() const {
	VecStr ret;
	{
		std::lock_guard<std::mutex> lock(mut_);
		for (const auto & state : states_) {
			ret.emplace_back(state.first);
		}
	}
	njh::sort(ret);
	return ret;
}

size_t SampleStateStore::sizeInBytes() const {
	std::lock_guard<std::mutex> lock(mut_);
	size_t ret = 0;
	for (const auto & state : states_) {
		ret += state.second.data_.size();
	}
	return ret;
}

size_t SampleStateStore::uncompressedSizeInBytes() const {
	std::lock_guard<std::mutex> lock(mut_);
	size_t ret = 0;
	for (const auto & state : states_) {
		ret += state.second.uncompressedSize_;
	}
	return ret;
}

}  // namespace collapse
}  // namespace njhseq
