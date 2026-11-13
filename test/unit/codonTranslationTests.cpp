#include <catch_amalgamated.hpp>
#include "njhseq/helpers/seqUtil.hpp"

using namespace njhseq;

TEST_CASE("standard genetic code", "[seqUtil][translation]") {
	const std::string bases = "TCAG";
	const std::string standardCode = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	for (uint32_t first = 0; first < 4; ++first) {
		for (uint32_t second = 0; second < 4; ++second) {
			for (uint32_t third = 0; third < 4; ++third) {
				const std::string codon{bases[first], bases[second], bases[third]};
				const std::string expected(1, standardCode[16 * first + 4 * second + third]);
				INFO(codon);
				CHECK(seqUtil::convertToProtein(codon) == expected);
				CHECK(seqUtil::convertOneCodon(codon) == expected);
				std::string rna = codon;
				std::replace(rna.begin(), rna.end(), 'T', 'U');
				CHECK(seqUtil::convertToProtein(rna) == expected);
				std::string lower = codon;
				std::transform(lower.begin(), lower.end(), lower.begin(), [](const char c) { return static_cast<char>(tolower(c)); });
				CHECK(seqUtil::convertToProtein(lower) == expected);
			}
		}
	}
}

TEST_CASE("ambiguous bases translate only when every expansion agrees", "[seqUtil][translation]") {
	//every expansion gives the same amino acid
	CHECK(seqUtil::convertToProtein("GCN") == "A");
	CHECK(seqUtil::convertToProtein("AAR") == "K");
	CHECK(seqUtil::convertToProtein("AAY") == "N");
	CHECK(seqUtil::convertToProtein("TAR") == "*");
	CHECK(seqUtil::convertToProtein("YTR") == "L");
	CHECK(seqUtil::convertToProtein("MGR") == "R");
	CHECK(seqUtil::convertToProtein("ATH") == "I");
	CHECK(seqUtil::convertToProtein("gcn") == "A");
	//expansions disagree
	CHECK(seqUtil::convertToProtein("NNN") == "X");
	CHECK(seqUtil::convertToProtein("ANN") == "X");
	CHECK(seqUtil::convertToProtein("AAN") == "X");
	CHECK(seqUtil::convertToProtein("TAN") == "X");
	CHECK(seqUtil::convertToProtein("TGN") == "X");
	CHECK(seqUtil::convertToProtein("RAT") == "X");
	CHECK(seqUtil::convertOneCodon("NNN") == "X");
}

TEST_CASE("gaps and non nucleotide characters give X", "[seqUtil][translation]") {
	CHECK(seqUtil::convertToProtein("A-T") == "X");
	CHECK(seqUtil::convertToProtein("---") == "X");
	CHECK(seqUtil::convertToProtein("AXG") == "X");
	CHECK(seqUtil::convertToProtein("A.G") == "X");
	CHECK(seqUtil::convertToProtein("ATG---GCT") == "MXA");
}

TEST_CASE("frame offset, trailing partial codon, and short input", "[seqUtil][translation]") {
	CHECK(seqUtil::convertToProtein("ATGGCTTAA") == "MA*");
	CHECK(seqUtil::convertToProtein("ATGGCTTA") == "MA");
	CHECK(seqUtil::convertToProtein("CATGGCTTAA", 1) == "MA*");
	CHECK(seqUtil::convertToProtein("CCATGGCTTAA", 2) == "MA*");
	CHECK(seqUtil::convertToProtein("AT").empty());
	CHECK(seqUtil::convertToProtein("").empty());
	CHECK(seqUtil::convertToProtein("ATGG", 2).empty());
	CHECK(seqUtil::convertToProtein("ATG", 5).empty());
	CHECK(seqUtil::convertOneCodon("AT").empty());
}

TEST_CASE("forceStartM only changes a known start codon at the start", "[seqUtil][translation]") {
	CHECK(seqUtil::convertToProtein("GTGGTG", 0, true) == "MV");
	CHECK(seqUtil::convertToProtein("TTGTTG", 0, true) == "ML");
	CHECK(seqUtil::convertToProtein("ATTATT", 0, true) == "MI");
	CHECK(seqUtil::convertToProtein("CTGCTG", 0, true) == "ML");
	CHECK(seqUtil::convertToProtein("CGTGGTG", 1, true) == "MV");
	CHECK(seqUtil::convertToProtein("GTGGTG", 0, false) == "VV");
	CHECK(seqUtil::convertToProtein("GCTGTG", 0, true) == "AV");
	CHECK(seqUtil::convertToProtein("NTGGTG", 0, true) == "XV");
}
