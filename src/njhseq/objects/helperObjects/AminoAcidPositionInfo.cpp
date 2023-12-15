//
// Created by Nicholas Hathaway on 12/15/23.
//

#include "AminoAcidPositionInfo.hpp"



namespace njhseq {

void AminoAcidPositionInfo::intialize(){

	VecStr oldColNames = infoTab_.columnNames_;
	njh::for_each(infoTab_.columnNames_, [](std::string&col) {
		njh::strToLower(col);
	});
	infoTab_.setColNamePositions();

	std::string idColName = "id";
	if (njh::in(std::string("transcriptid"), infoTab_.columnNames_)) {
		idColName = "transcriptid";
	}
	std::string aaPosColName = "aaposition";
	bool byRange = false;
	VecStr requiredColumns{idColName, aaPosColName};
	if (njh::in(std::string("aastart"), infoTab_.columnNames_) && njh::in(std::string("aastop"), infoTab_.columnNames_)) {
		requiredColumns = VecStr{idColName, "aastart", "aastop"};
		byRange = true;
	}
	infoTab_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);

	auto idsVec = infoTab_.getColumnLevels(idColName);
	ids_ = std::set<std::string>(idsVec.begin(), idsVec.end());

	VecStr additionalColumns;
	for (const auto&col: infoTab_.columnNames_) {
		if (!njh::in(col, requiredColumns)) {
			additionalColumns.emplace_back(col);
		}
	}
	std::vector<uint32_t> addColPositions;
	for (const auto&col: additionalColumns) {
		addColPositions.emplace_back(infoTab_.getColPos(col));
	}
	if (byRange) {
		for (const auto&row: infoTab_.content_) {
			auto aaStart =
					zeroBased_
						? njh::StrToNumConverter::stoToNum<uint32_t>(
							row[infoTab_.getColPos("aastart")])
						: njh::StrToNumConverter::stoToNum<uint32_t>(
							  row[infoTab_.getColPos("aastart")]) - 1;
			auto aastop =
					zeroBased_
						? njh::StrToNumConverter::stoToNum<uint32_t>(
							row[infoTab_.getColPos("aastop")])
						: njh::StrToNumConverter::stoToNum<uint32_t>(
							  row[infoTab_.getColPos("aastop")]) + 1;

			auto idName = row[infoTab_.getColPos(idColName)];

			if (aastop < aaStart) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "error for " << idName << " stop: " << aastop <<
						" can't be lower than start: " << aaStart << "\n";
				throw std::runtime_error{ss.str()};
			}

			for (const auto&aaPos: iter::range(aaStart, aastop)) {
				aminoPositionsPerId_[idName].emplace(aaPos);
				//add any additional column as meta data
				MetaDataInName meta;
				for (const auto&pos: addColPositions) {
					meta.addMeta(oldColNames[pos], row[pos]);
				}
				metaDataForAAPos_[idName][aaPos] = meta;
			}
		}
	}
	else {
		for (const auto&row: infoTab_.content_) {
			auto aaPos =
					zeroBased_
						? njh::StrToNumConverter::stoToNum<uint32_t>(
							row[infoTab_.getColPos(aaPosColName)])
						: njh::StrToNumConverter::stoToNum<uint32_t>(
							  row[infoTab_.getColPos(aaPosColName)]) - 1;
			auto idName = row[infoTab_.getColPos(idColName)];
			aminoPositionsPerId_[idName].emplace(aaPos);
			//add any additional column as meta data
			MetaDataInName meta;
			for (const auto&pos: addColPositions) {
				meta.addMeta(oldColNames[pos], row[pos]);
			}
			metaDataForAAPos_[idName][aaPos] = meta;
		}
	}
}


AminoAcidPositionInfo::AminoAcidPositionInfo(table infoTab,
                                             bool zeroBased) :
                                                               infoTab_(std::move(infoTab)),
                                                               zeroBased_(zeroBased) {
	intialize();
}

AminoAcidPositionInfo::AminoAcidPositionInfo(const bfs::path&inputInfoFnp,
                                             bool zeroBased) : inputInfoFnp_(inputInfoFnp),
                                                               infoTab_(inputInfoFnp, "\t", true),
                                                               zeroBased_(zeroBased) {
	intialize();
}
bool AminoAcidPositionInfo::byRange() const {
  return njh::in(std::string("aastart"), infoTab_.columnNames_) && njh::in(std::string("aastop"), infoTab_.columnNames_);
}

} //namespace njhseq


