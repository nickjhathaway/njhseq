#include <catch_amalgamated.hpp>
#include "njhseq/GenomeUtils/GenomeMapping/MultiGenomeMapper.hpp"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"

using namespace njhseq;

//nothing else uses these so without this they aren't compiled and stop building without anyone noticing
template std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>>
MultiGenomeMapper::addIntersectingGeneInfosToLocs<Bed6RecordCore>(
	std::vector<Bed6RecordCore>&, const std::string&, const std::string&, bool);

template std::unordered_map<uint32_t, std::vector<MultiGenomeMapper::IntersectedProteinInfo>>
MultiGenomeMapper::addIntersectingGeneInfosToLocs<Bed6RecordCore>(
	std::vector<Bed6RecordCore>&, const intersectBedLocsWtihGffRecordsPars&, const bfs::path&, bool);

namespace {

const uint32_t noAA = std::numeric_limits<uint32_t>::max();

/**
 * \brief make the positions of a transcript from the genomic positions of its coding bases, in coding order
 */
std::unordered_map<uint32_t, GeneSeqInfo::GenePosInfo> makePosInfos(const std::vector<uint32_t>& codingPositions) {
	std::unordered_map<uint32_t, GeneSeqInfo::GenePosInfo> ret;
	uint32_t cDNAPos = 0;
	for (const auto pos : codingPositions) {
		GeneSeqInfo::GenePosInfo info;
		info.gDNAPos_ = pos;
		info.cDNAPos_ = cDNAPos;
		info.aaPos_ = cDNAPos / 3;
		info.codonPos_ = cDNAPos % 3;
		ret[pos] = info;
		++cDNAPos;
	}
	return ret;
}

/**
 * \brief the coding positions [start, end), walked backwards for a reverse strand gene
 */
std::vector<uint32_t> codingRange(const uint32_t start, const uint32_t end, const bool reverse = false) {
	std::vector<uint32_t> ret;
	for (uint32_t pos = start; pos < end; ++pos) {
		ret.emplace_back(pos);
	}
	if (reverse) {
		std::reverse(ret.begin(), ret.end());
	}
	return ret;
}

}  // namespace

TEST_CASE("getAARangeForRegion any base versus full codon", "[getAARangeForRegion]") {
	//12 bases, 4 codons: 100-102 AA1, 103-105 AA2, 106-108 AA3, 109-111 AA4
	const auto posInfos = makePosInfos(codingRange(100, 112));
	using range = std::pair<uint32_t, uint32_t>;

	SECTION("region on codon boundaries is the same either way") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 109, false) == range{2, 3});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 109, true) == range{2, 3});
	}
	SECTION("partial codons at both ends only count with any base") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 102, 107, false) == range{1, 3});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 102, 107, true) == range{2, 2});
	}
	SECTION("1 base short of a codon") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 105, false) == range{2, 2});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 105, true) == range{noAA, noAA});
	}
	SECTION("defaults to any base") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 102, 107) == range{1, 3});
	}
	SECTION("region past the end of the transcript") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 110, 200, false) == range{4, 4});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 110, 200, true) == range{noAA, noAA});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 109, 200, true) == range{4, 4});
	}
	SECTION("region outside of the transcript") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 10, 50, false) == range{noAA, noAA});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 10, 50, true) == range{noAA, noAA});
	}
}

TEST_CASE("getAARangeForRegion reverse strand", "[getAARangeForRegion]") {
	//AA1 is at the highest genomic positions: 109-111 AA1, 106-108 AA2, 103-105 AA3, 100-102 AA4
	const auto posInfos = makePosInfos(codingRange(100, 112, true));
	using range = std::pair<uint32_t, uint32_t>;

	REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 102, 107, false) == range{2, 4});
	REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 102, 107, true) == range{3, 3});
	REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 109, true) == range{2, 3});
}

TEST_CASE("getAARangeForRegion with introns", "[getAARangeForRegion]") {
	//exon 100-103, intron 104-199, exon 200-207, the intron splits AA2: 100-102 AA1, 103,200,201 AA2, 202-204 AA3, 205-207 AA4
	auto coding = codingRange(100, 104);
	const auto secondExon = codingRange(200, 208);
	coding.insert(coding.end(), secondExon.begin(), secondExon.end());
	auto posInfos = makePosInfos(coding);
	//introns are in the real infos but have no amino acid
	for (uint32_t pos = 104; pos < 200; ++pos) {
		posInfos[pos].gDNAPos_ = pos;
	}
	using range = std::pair<uint32_t, uint32_t>;

	SECTION("codon split by the intron is complete when the region spans it") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 202, true) == range{2, 2});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 103, 202, false) == range{2, 2});
	}
	SECTION("codon split by the intron is partial when the region stops in the intron") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 100, 150, true) == range{1, 1});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 100, 150, false) == range{1, 2});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 150, 205, true) == range{3, 3});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 150, 205, false) == range{2, 3});
	}
	SECTION("completely within the intron") {
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 110, 190, false) == range{noAA, noAA});
		REQUIRE(MultiGenomeMapper::getAARangeForRegion(posInfos, 110, 190, true) == range{noAA, noAA});
	}
}
