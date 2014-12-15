#include <catch.hpp>

#include "../src/bibseq/objects/counters/charCounter.hpp"
using namespace bibseq;

TEST_CASE("Basic tests for charCounterArray "){
  SECTION("correct count and fractions set for DNA library"){
  	std::string testStr = std::string(10,'A') + std::string(20, 'C') + std::string(30,'G') + std::string(40, 'T');
  	charCounterArray counter(std::vector<char> {'A', 'C', 'G', 'T'});
  	counter.increaseCountByString(testStr);
  	counter.setFractions(counter.alphabet_);

  	//check the counts
  	REQUIRE(10 == counter.chars_['A']);
  	REQUIRE(20 == counter.chars_['C']);
  	REQUIRE(30 == counter.chars_['G']);
  	REQUIRE(40 == counter.chars_['T']);
  	//check the fractions
  	REQUIRE((10/100.0) == counter.fractions_['A']);
  	REQUIRE((20/100.0) == counter.fractions_['C']);
  	REQUIRE((30/100.0) == counter.fractions_['G']);
  	REQUIRE((40/100.0) == counter.fractions_['T']);

  	//
  	counter.reset();
  	counter.increaseCountByString(testStr);
  	std::vector<char> all(127);
  	std::iota(all.begin(), all.end(), 0);
  	counter.setFractions(all);
  	//check the counts again to make sure reset works
  	REQUIRE(10 == counter.chars_['A']);
  	REQUIRE(20 == counter.chars_['C']);
  	REQUIRE(30 == counter.chars_['G']);
  	REQUIRE(40 == counter.chars_['T']);
  	//check the fractions after setting by all
  	REQUIRE((10/100.0) == counter.fractions_['A']);
  	REQUIRE((20/100.0) == counter.fractions_['C']);
  	REQUIRE((30/100.0) == counter.fractions_['G']);
  	REQUIRE((40/100.0) == counter.fractions_['T']);


  }
  SECTION("correct GC content DNA library"){
  	std::string testStr = std::string(10,'A') + std::string(20, 'C') + std::string(30,'G') + std::string(40, 'T');
  	charCounterArray counter(std::vector<char> {'A', 'C', 'G', 'T'});
  	counter.increaseCountByString(testStr);
  	counter.setFractions(counter.alphabet_);
  	counter.calcGcContent();
  	//check the counts
  	REQUIRE((50/100.0) == counter.gcContent);


  }

}

