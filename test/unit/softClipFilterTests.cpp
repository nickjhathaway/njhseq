#include <catch_amalgamated.hpp>
#include "njhseq/BamToolsUtils/SoftClipFilter.hpp"

using namespace njhseq;

namespace {

//Pf3D7_04_v3:283313-284016 (LISP2), has a 12 bp tandem repeat (ATTGTGATAATC-like units) at about 210-245
const std::string lisp2Ref =
		"CAATAAAATCTTCCTCTTCTTTAAATGATTCGTCCTGACCATCCTGCACCTCATGATAATCATCAACATTATATTCATAATTAATAACGTCTTGATTATT"
		"TTTATTTTCATGTTTAATATTATAAAATTTTAAATTATTAAAAGAGATATTGTATTTATTCTTAAGAAATGTATCATGTATTTTGTTTGTGTGGTTTCCT"
		"TTACTAGAATGATTATTGTGATAATCGTTGTGATAATTATTATTGGGATGATTATTTCTATTGTAAAAATTATTTATATCATTTATTTTAAAATTTTCAT"
		"TATAATTATTAACATAATTTTTATTAAATGACATATCATTTTCGTCTTCATCTTTCTTTTCTTTAATTATATGATTCATTTTTTTTTCGTGATATATTTC"
		"CTTTTTATTATTTATAATATCATTAAATTCAGTTTGTAATTTTTTGGTTTTATATATTTTCTCTTCTTCAACTTTAATAATAATTTTACCTGTATAAACA"
		"ATATTAGGTTCTGTAAAATAATTACAAGTACAAGAGATTGTTTTATCTTGGTTCGTCTCATTTAATAAAATATATCTGGTCTGATTAACATCTAAATCAT"
		"GATGTATCAAATCAAATCTAGGAAATGGAACCGCTCCTGGTAATATATCTTCCACATTTATGTTATCATTTACAAATTCAAAACAATTAATTGGTTCTGT"
		"TTG";

std::string reversed(const std::string & seq) {
	return std::string(seq.rbegin(), seq.rend());
}

BamTools::BamAlignment makeAln(const std::string & seq, int32_t pos, const std::vector<BamTools::CigarOp> & cigar) {
	BamTools::BamAlignment aln;
	aln.Name = "read";
	aln.QueryBases = seq;
	aln.Qualities = std::string(seq.size(), 'I');
	aln.RefID = 0;
	aln.Position = pos;
	aln.CigarData = cigar;
	aln.MapQuality = 60;
	aln.AlignmentFlag = 0;
	return aln;
}

std::string randomSeq(uint32_t length, uint32_t seed) {
	std::mt19937 gen(seed);
	std::uniform_int_distribution<int> base(0, 3);
	std::string ret(length, 'A');
	for (auto & letter : ret) {
		letter = "ACGT"[base(gen)];
	}
	return ret;
}

}  // namespace

TEST_CASE("a soft clipped LISP2 repeat expansion that runs into the far flank is a tandem extension", "[SoftClipTandemExtension]") {
	//1S1-1C1-10000, aligned 60S85M at 283529, four extra 12 bp units plus the left flank were soft clipped
	const std::string read = "CTAGAATGATTATTGTGATAATCATTGTGATAGTTGTTGTGATAGTTGTTGTGATAATTGTTGTGATAATTGTTGTGATAATTATTATTGGGATGATTATTTCTATTGTAAAAATTATTTATATCATTTATTTTAAAATTTTCAT";
	SoftClipTandemExtension::Pars pars;
	auto res = SoftClipTandemExtension::classifyClip(read.substr(0, 60), read.substr(60), lisp2Ref, pars);
	REQUIRE(res.pass_);
	REQUIRE(res.period_ == 12);
	//the flank's own ATTATTGTGATAAT carries the period, only 5 bases are left past the tract
	REQUIRE(res.tractLen_ == 55);
	REQUIRE_FALSE(res.usedFlank_);

	//with more of the far flank in the clip what's past the tract has to match the reference
	const std::string moreFlank = "TTTGTGTGGTTTCCTTTA";
	auto withFlank = SoftClipTandemExtension::classifyClip(moreFlank + read.substr(0, 60), read.substr(60), lisp2Ref, pars);
	REQUIRE(withFlank.pass_);
	REQUIRE(withFlank.usedFlank_);
	auto wrongFlank = SoftClipTandemExtension::classifyClip("GCGCGCATATGCGCGCATAT" + read.substr(0, 60), read.substr(60), lisp2Ref, pars);
	REQUIRE_FALSE(wrongFlank.pass_);
}

TEST_CASE("a trailing clip that is the reverse complement of the flank (a foldback) is not a tandem extension", "[SoftClipTandemExtension]") {
	//1S1-1C1-10000, aligned 57M61S at 283595
	const std::string read = "TTTATTTTAAAATTTTCATTATAATTATTAACATAATTTTTATTAAATGACATATCACAACAACTATCACAATGATTATCACAATAATCATTCTAGTAAAGGAAACCACCCAAACAAA";
	SoftClipTandemExtension::Pars pars;
	const auto rRead = reversed(read);
	auto res = SoftClipTandemExtension::classifyClip(rRead.substr(0, 61), rRead.substr(61), reversed(lisp2Ref), pars);
	REQUIRE_FALSE(res.pass_);
}

TEST_CASE("a clip from somewhere else is not a tandem extension", "[SoftClipTandemExtension]") {
	SoftClipTandemExtension::Pars pars;
	const std::string anchor = lisp2Ref.substr(220, 90);
	auto res = SoftClipTandemExtension::classifyClip(randomSeq(55, 7), anchor, lisp2Ref, pars);
	REQUIRE_FALSE(res.pass_);
}

TEST_CASE("a low complexity clip next to a low complexity anchor is held to the entropy cut off", "[SoftClipTandemExtension]") {
	SoftClipTandemExtension::Pars pars;
	const std::string clip(40, 'T');
	const std::string anchor = std::string(30, 'T') + "CCCTTGTCGACCTTAACAGATGGCTCACG";
	auto res = SoftClipTandemExtension::classifyClip(clip, anchor, std::string(600, 'T'), pars);
	REQUIRE_FALSE(res.pass_);
	pars.minEntropy_ = 0;
	REQUIRE(SoftClipTandemExtension::classifyClip(clip, anchor, std::string(600, 'T'), pars).pass_);
}

TEST_CASE("trimming to a region keeps soft clipped bases that sit inside the region", "[trimAlnSeqToRegion]") {
	BamTools::RefVector refData{BamTools::RefData("chr", 10000)};
	GenomicRegion region("reg", "chr", 100, 400, false);
	const std::string seq = randomSeq(100, 11);
	using BamTools::CigarOp;

	SECTION("leading clip well inside the region is kept whole") {
		auto aln = makeAln(seq, 200, {CigarOp('S', 20), CigarOp('M', 80)});
		REQUIRE(trimAlnSeqToRegion(aln, region, refData, true).seq_ == seq);
		REQUIRE(trimAlnSeqToRegion(aln, region, refData, false).seq_ == seq.substr(20));
	}
	SECTION("leading clip near the region start keeps only what would be inside") {
		auto aln = makeAln(seq, 105, {CigarOp('S', 20), CigarOp('M', 80)});
		REQUIRE(trimAlnSeqToRegion(aln, region, refData, true).seq_ == seq.substr(15));
	}
	SECTION("trailing clip inside the region is kept, and qualities follow") {
		auto aln = makeAln(seq, 200, {CigarOp('M', 70), CigarOp('S', 30)});
		auto out = trimAlnSeqToRegion(aln, region, refData, true);
		REQUIRE(out.seq_ == seq);
		REQUIRE(out.qual_.size() == seq.size());
		REQUIRE(out.qual_.back() == 40);
	}
	SECTION("an alignment crossing the region start is still trimmed there") {
		auto aln = makeAln(seq, 50, {CigarOp('M', 100)});
		REQUIRE(trimAlnSeqToRegion(aln, region, refData, true).seq_ == seq.substr(50));
	}
	SECTION("each clip is kept or dropped by the keep function") {
		auto aln = makeAln(seq, 200, {CigarOp('S', 20), CigarOp('M', 60), CigarOp('S', 20)});
		std::function<bool(const BamTools::BamAlignment &, bool)> leadingOnly = [](const BamTools::BamAlignment &, bool leading) { return leading; };
		REQUIRE(trimAlnSeqToRegion(aln, region, refData, leadingOnly).seq_ == seq.substr(0, 80));
		REQUIRE(trimAlnSeqToRegion(aln, region, refData, std::function<bool(const BamTools::BamAlignment &, bool)>{}).seq_ == seq.substr(20, 60));
	}
	SECTION("explained clips without a genome are all dropped") {
		SoftClipFilterPars pars;
		pars.trimClips_ = SoftClipFilterPars::TrimClips::explained;
		SoftClipFilter filter(pars, 30, 1, refData);
		auto aln = makeAln(seq, 200, {CigarOp('S', 20), CigarOp('M', 80)});
		REQUIRE_FALSE(filter.keepClipWhenTrimming(aln, true));
	}
	SECTION("clip off the far side of an alignment crossing the edge is dropped") {
		auto aln = makeAln(seq, 380, {CigarOp('S', 10), CigarOp('M', 60), CigarOp('S', 30)});
		auto out = trimAlnSeqToRegion(aln, region, refData, true);
		REQUIRE(out.seq_ == seq.substr(0, 10 + 20));
	}
}

TEST_CASE("soft clip filter modes", "[SoftClipFilter]") {
	BamTools::RefVector refData{BamTools::RefData("chr", 100000), BamTools::RefData("chr2", 100000)};
	GenomicRegion region("reg", "chr", 1000, 2000, false);
	using BamTools::CigarOp;
	const std::string seq = randomSeq(150, 3);
	auto clipped = makeAln(seq, 1500, {CigarOp('S', 60), CigarOp('M', 90)});
	auto unclipped = makeAln(seq, 1500, {CigarOp('S', 10), CigarOp('M', 140)});

	SoftClipFilterPars pars;
	SECTION("hardCutOff filters on amount only") {
		SoftClipFilter filter(pars, 30, 1, refData);
		REQUIRE(filter.passes(unclipped, region));
		REQUIRE_FALSE(filter.passes(clipped, region));
	}
	SECTION("saElsewhere keeps clips without a supplementary alignment or with one in the region") {
		pars.mode_ = SoftClipFilterPars::Mode::saElsewhere;
		SoftClipFilter filter(pars, 30, 1, refData);
		REQUIRE(filter.passes(clipped, region));
		auto sameLocus = clipped;
		sameLocus.AddTag("SA", "Z", std::string("chr,1450,-,90S60M,60,0;"));
		REQUIRE(filter.passes(sameLocus, region));
		auto elsewhere = clipped;
		elsewhere.AddTag("SA", "Z", std::string("chr2,50000,+,90S60M,60,0;"));
		REQUIRE_FALSE(filter.passes(elsewhere, region));
		auto farOnSameChrom = clipped;
		farOnSameChrom.AddTag("SA", "Z", std::string("chr,9000,+,90S60M,60,0;"));
		REQUIRE_FALSE(filter.passes(farOnSameChrom, region));
	}
	SECTION("saElsewhereOrInverted also filters an opposite strand supplementary alignment in the region") {
		pars.mode_ = SoftClipFilterPars::Mode::saElsewhereOrInverted;
		SoftClipFilter filter(pars, 30, 1, refData);
		auto inverted = clipped;
		inverted.AddTag("SA", "Z", std::string("chr,1450,-,90S60M,60,0;"));
		REQUIRE_FALSE(filter.passes(inverted, region));
		auto sameStrand = clipped;
		sameStrand.AddTag("SA", "Z", std::string("chr,1450,+,90S60M,60,0;"));
		REQUIRE(filter.passes(sameStrand, region));
	}
	SECTION("tandem needs a genome") {
		pars.mode_ = SoftClipFilterPars::Mode::tandem;
		REQUIRE_THROWS(SoftClipFilter(pars, 30, 1, refData));
	}
}
