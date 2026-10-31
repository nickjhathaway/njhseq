/*
 * KmerCollectionIO.cpp
 *
 *  Created on: Oct 31, 2026
 *      Author: nick
 */

#include "KmerCollectionIO.hpp"
#include "njhseq/objects/dataContainers/tables/TableReader.hpp"

namespace njhseq {

namespace {
void checkTwoColumns(const TableReader & reader, const std::string & funcName) {
	if (reader.header_.nCol() < 2) {
		std::stringstream ss;
		ss << funcName << ", error " << "need to have 2 columns" << "\n";
		throw std::runtime_error{ss.str()};
	}
}
}  // namespace

uint32_t KmerCollectionIO::getKmerLength(const bfs::path & tableFnp) {
	auto toks = tokenizeString(njh::files::getFirstLine(tableFnp), "\t");
	return toks.at(1).size();
}

KmerCollectionIO::HashedCollections KmerCollectionIO::readHashedCollections(const bfs::path & tableFnp) {
	uint32_t kmerLength = 0;
	return readHashedCollections(tableFnp, kmerLength);
}

KmerCollectionIO::HashedCollections KmerCollectionIO::readHashedCollections(const bfs::path & tableFnp, uint32_t & kmerLength) {
	HashedCollections ret;
	kmerLength = 0;
	SimpleKmerHash hasher;
	TableReader reader(TableIOOpts::genTabFileIn(tableFnp, false));
	checkTwoColumns(reader, __PRETTY_FUNCTION__);
	VecStr row;
	while (reader.getNextRow(row)) {
		kmerLength = row[1].size();
		ret[row[0]].emplace(hasher.hash(row[1]));
	}
	return ret;
}

std::unordered_set<uint64_t> KmerCollectionIO::readHashedKmers(const bfs::path & tableFnp) {
	std::unordered_set<uint64_t> ret;
	SimpleKmerHash hasher;
	TableReader reader(TableIOOpts::genTabFileIn(tableFnp, false));
	checkTwoColumns(reader, __PRETTY_FUNCTION__);
	VecStr row;
	while (reader.getNextRow(row)) {
		ret.emplace(hasher.hash(row[1]));
	}
	return ret;
}

}  // namespace njhseq
