/*
 * MetaDataInName.cpp
 *
 *  Created on: Jan 27, 2017
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

#include "MetaDataInName.hpp"
#include "njhseq/utils/stringUtils.hpp"

namespace njhseq {


MetaDataInName::MetaDataInName() {

}

MetaDataInName::MetaDataInName(const std::string & str) {
	processNameForMeta(str, false);
}


void MetaDataInName::removeMeta(const std::string & metaField){
	containsMetaThrow(metaField, __PRETTY_FUNCTION__);
	meta_.erase(metaField);
}

void MetaDataInName::addMeta(const MetaDataInName & otherMeta, bool replace) {
	for(const auto & meta : otherMeta.meta_){
		addMeta(meta.first, meta.second, replace);
	}
}


void MetaDataInName::processNameForMeta(const std::string & name, bool replace){

	auto firstBracket = name.find("[");
	if(std::string::npos == firstBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", could not find [ in " << name << std::endl;
		throw std::runtime_error{ss.str()};
	}
	auto secondBracket = name.rfind("]");
	if(std::string::npos == secondBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", could not find ] in " << name
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	if(firstBracket > secondBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", [ must come before ] "
				<< "\n";
		ss << "[ pos: " << firstBracket << "; ] pos: " << secondBracket << "\n"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	auto openingBrackets = findOccurences(name, "[");
	auto closingBrackets = findOccurences(name, "]");
	if(openingBrackets.size() > 1){
		/**@todo still some corner cases that need to be handled below, like if there is two opening brackets before the close e.g. [text[text] ]*/
		if(openingBrackets.size() != closingBrackets.size()){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__
					<< " unequal number of opening and closing brackets" << " from name: " << name
					<< ", opening:" << openingBrackets.size() << ", closing:" <<closingBrackets.size()
					<< std::endl;
			throw std::runtime_error{ss.str()};
		}

		auto nextBracketPos = name.find("[", firstBracket + 1);
		std::vector<std::pair<size_t, size_t>> bracketsPositions;
		while(std::string::npos != nextBracketPos ){
			auto nextCloseBracket = name.find("]", nextBracketPos);
			uint32_t openBracketCount = 0;
			auto theNextOpeningBracket =  name.find("[", nextBracketPos + 1);
			while(theNextOpeningBracket < nextCloseBracket){
				++openBracketCount;
				theNextOpeningBracket =  name.find("[", theNextOpeningBracket + 1);
			}
			for(uint32_t bracketCount = 0; bracketCount < openBracketCount; ++bracketCount){
				nextCloseBracket = name.find("]", nextCloseBracket + 1);
			}
			if(std::string::npos == nextCloseBracket){
				std::stringstream ss;
				ss << "Error in : " << __PRETTY_FUNCTION__
						<< "could not find the next closing bracket after open bracketPosition: "
						<< nextBracketPos << " from name: " << name << std::endl;
				throw std::runtime_error { ss.str() };
			}
			bracketsPositions.emplace_back(std::make_pair(nextBracketPos, nextCloseBracket));
			nextBracketPos = name.find("[", nextCloseBracket);
		}
		std::unordered_map<std::string, std::string> extractKeyToks;
		uint32_t tokNum = 0;
		auto modName = name;
		for(const auto & position : iter::reversed(bracketsPositions)){
			std::string replacementTokenName = "TOKEN_TO_BE_REPLACED." + estd::to_string(tokNum);
			extractKeyToks[replacementTokenName] = name.substr(position.first, position.second + 1 - position.first);
			modName.erase(modName.begin() + position.first, modName.begin() + 1 + position.second);
			modName.insert(modName.begin() + position.first, replacementTokenName.begin(), replacementTokenName.end());
			++tokNum;
		}
		auto modName_firstBracket = modName.find("[");
		if(std::string::npos == modName_firstBracket){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__
					<< ", could not find [ in " << modName << std::endl;
			throw std::runtime_error{ss.str()};
		}
		auto modName_secondBracket = modName.rfind("]");
		if(std::string::npos == modName_secondBracket){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__
					<< ", could not find ] in " << modName  << " after " << modName_firstBracket
					<< std::endl;
			throw std::runtime_error{ss.str()};
		}
		if(modName_firstBracket > modName_secondBracket){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__ << ", [ must come before ] "
					<< "\n";
			ss << "[ pos: " << modName_firstBracket << "; ] pos: " << modName_secondBracket << "\n"
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		auto toks = njh::tokenizeString(modName.substr(modName_firstBracket + 1, modName_secondBracket - modName_firstBracket - 1), ";");
		for(const auto & tok : toks){
//			auto subToks = njh::tokenizeString(tok, "=");
//			if(2 != subToks.size()){
//				std::stringstream ss;
//				ss << "Error in : " << __PRETTY_FUNCTION__
//						<< "values should be separated by one =, not " << subToks.size()
//						<< " for tok: " << tok << " from modName: " << modName
//						<< ", originalName: " << name << std::endl;
//				throw std::runtime_error { ss.str() };
//			}else{
//				auto metaFieldKey = subToks[0];
//				auto metaFieldValue = subToks[1];
//				for(const auto & nameReplace : extractKeyToks){
//					metaFieldKey = njh::replaceString(metaFieldKey, nameReplace.first, nameReplace.second);
//					metaFieldValue = njh::replaceString(metaFieldValue, nameReplace.first, nameReplace.second);
//				}
//				addMeta(metaFieldKey, metaFieldValue, replace);
//			}
			auto firstEqual = tok.find("=");
			if(std::string::npos == firstEqual){
				std::stringstream ss;
				ss << "Error in : " << __PRETTY_FUNCTION__
						<< "values should be separated by one =, no equal sign found in  " << "tok: " << tok << " from modName: " << modName
						<< ", originalName: " << name << std::endl;
				throw std::runtime_error { ss.str() };
			}else{
				auto metaFieldKey = tok.substr(0, tok.find("="));
				auto metaFieldValue = firstEqual == tok.size()? std::string("") : tok.substr(tok.find("=") + 1);
				for(const auto & nameReplace : extractKeyToks){
					metaFieldKey = njh::replaceString(metaFieldKey, nameReplace.first, nameReplace.second);
					metaFieldValue = njh::replaceString(metaFieldValue, nameReplace.first, nameReplace.second);
				}
				addMeta(metaFieldKey, metaFieldValue, replace);
			}
		}
	}else{
		auto toks = njh::tokenizeString(name.substr(firstBracket + 1, secondBracket - firstBracket - 1), ";");
		for(const auto & tok : toks){
//			auto subToks = njh::tokenizeString(tok, "=");
//			if(2 != subToks.size()){
//				std::stringstream ss;
//				ss << "Error in : " << __PRETTY_FUNCTION__
//						<< "values should be separated by one =, not " << subToks.size() << " for tok: " << tok << " from name: " << name
//						<< std::endl;
//				throw std::runtime_error{ss.str()};
//			}else{
//				addMeta(subToks[0], subToks[1], replace);
//			}
			auto firstEqual = tok.find("=");
			if(std::string::npos == firstEqual){

				std::stringstream ss;
				ss << "Error in : " << __PRETTY_FUNCTION__
						<< "values should be separated by one =, no = found in tok: " << tok << " from name: " << name
						<< std::endl;
				throw std::runtime_error{ss.str()};
			}else{
				addMeta(tok.substr(0, tok.find("=")), firstEqual == tok.size() ? std::string("") : tok.substr(tok.find("=") + 1), replace);
			}
		}
	}
}

bool MetaDataInName::containsMeta(const std::string & key) const {
	return meta_.find(key) != meta_.end();
}

void MetaDataInName::containsMetaThrow(const std::string & key, const std::string & funcName) const{
	if(!containsMeta(key)){
		std::stringstream ss;
		ss << funcName << ", error no meta field " << key << "\n";
		ss << "Options are: " << njh::conToStr(njh::getVecOfMapKeys(meta_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
}

std::string MetaDataInName::getMeta(const std::string & key) const {
	auto search = meta_.find(key);
	if (search != meta_.end()) {
		return search->second;
	}else{
		std::stringstream ss;
		ss << __FILE__ << " - " << __LINE__ << " : " << __PRETTY_FUNCTION__ << ", error no meta field " << key << "\n";
		ss << "Options are: " << njh::conToStr(njh::getVecOfMapKeys(meta_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return "";
}

std::string MetaDataInName::pasteLevels(const std::string & sep) const{
	return pasteLevels(getVectorOfMapKeys(meta_), sep);
}

std::string MetaDataInName::pasteLevels(const VecStr & metalevels,
		const std::string & sep) const {
	VecStr levelVals;
	for(const auto & lev : metalevels){
		levelVals.emplace_back(getMeta(lev));
	}
	if ("" != sep) {
		return njh::conToStr(levelVals, sep);
	} else {
		return njh::pasteAsStr(levelVals);
	}
}

std::string MetaDataInName::createMetaName() const {
	std::string newMeta = "[";
	auto metaKeys = njh::getVecOfMapKeys(meta_);
	//sort integer by their actual numerical values if all keys numbers
	if(std::all_of(metaKeys.begin(), metaKeys.end(), [](const std::string & str){ return njh::strAllDigits(str);}) ){
		njh::sort(metaKeys, [](const std::string & str1, const std::string & str2){
			return njh::StrToNumConverter::stoToNum<uint32_t>(str1) < njh::StrToNumConverter::stoToNum<uint32_t>(str2);
		});
	}else{
		njh::sort(metaKeys);
	}
	for (const auto & metaKey : metaKeys) {
		const auto & meta = meta_.at(metaKey);
		if ("[" != newMeta) {
			newMeta.append(";" + metaKey + "=" + meta);
		} else {
			newMeta.append(metaKey + "=" + meta);
		}
	}
	newMeta.append("]");
	return newMeta;
}

std::string MetaDataInName::createMetaName(const std::function<bool(const std::string &, const std::string &)> & metaKeyPredSorter) const{
	std::string newMeta = "[";
	auto metaKeys = njh::getVecOfMapKeys(meta_);
	njh::sort(metaKeys, metaKeyPredSorter);
	for (const auto & metaKey : metaKeys) {
		const auto & meta = meta_.at(metaKey);
		if ("[" != newMeta) {
			newMeta.append(";" + metaKey + "=" + meta);
		} else {
			newMeta.append(metaKey + "=" + meta);
		}
	}
	newMeta.append("]");
	return newMeta;
}

void MetaDataInName::resetMetaInName(std::string & name,
		size_t pos) const{
	auto firstOpeningBracket = name.find("[");
	auto secondBracket = name.rfind("]");
	if(firstOpeningBracket > secondBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", [ must come before ] "
				<< "\n";
		ss << "[ pos: " << firstOpeningBracket << "; ] pos: " << secondBracket << "\n"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::string newMeta = createMetaName();
	if (std::string::npos != firstOpeningBracket
			&& std::string::npos != secondBracket) {
		name = name.substr(0, firstOpeningBracket) + newMeta
				+ name.substr(secondBracket + 1);
	} else {
		if (std::numeric_limits<size_t>::max() != pos && pos < name.size()) {
			name.insert(name.begin() + pos, newMeta.begin(), newMeta.end());
		} else {
			name += newMeta;
		}
	}
}

void MetaDataInName::removeMetaDataInName(std::string & name){
	if(nameHasMetaData(name)){
		auto firstBracket = name.find("[");
		auto secondBracket = name.rfind("]");
		if(firstBracket > secondBracket){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__ << ", [ must come before ] "
					<< "\n";
			ss << "[ pos: " << firstBracket << "; ] pos: " << secondBracket << "\n"
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		if (std::string::npos != firstBracket
				&& std::string::npos != secondBracket) {
			name = name.substr(0, firstBracket) + name.substr(secondBracket + 1);
		}
	}
}

bool MetaDataInName::nameHasMetaData(const std::string & name) {
	auto firstBracket = name.find("[");
	if (std::string::npos == firstBracket) {
		return false;
	}
	auto secondBracket = name.rfind("]");
	if (std::string::npos == secondBracket) {
		return false;
	}
	if (firstBracket > secondBracket && secondBracket != (firstBracket + 1)) {
		return false;
	}
	return true;
}


MetaDataInName MetaDataInName::genMetaFromJson(const Json::Value & val){
	MetaDataInName ret;
	njh::json::MemberChecker memChecker(val);
	memChecker.failMemberCheckThrow(VecStr{"meta_"}, __PRETTY_FUNCTION__);
	for(const auto & metaName : val["meta_"].getMemberNames()){
		ret.addMeta(metaName, val["meta_"][metaName]);
	}
	return ret;
}


Json::Value MetaDataInName::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["meta_"] = njh::json::toJson(meta_);
	return ret;
}

}  // namespace njhseq


