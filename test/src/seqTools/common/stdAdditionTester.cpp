#include <catch.hpp>

#include "../src/bibseq/common/stdAddition.hpp"
using namespace bibseq;

TEST_CASE("Basic tests for stdAddition", "[stdAddition]" ){
    SECTION("contains/has/in"){
        std::vector<int> v = {0,2,4,6,8};
        REQUIRE(contains(v, 2));
        REQUIRE(!contains(v, 3));
        REQUIRE(has(v, 2));
        REQUIRE(!has(v, 3));
        REQUIRE(in(2, v));
        REQUIRE(!in(3, v));
    }

    SECTION("std::vec"){
        REQUIRE(contains({0,2,4,6,8}, 2));
        REQUIRE(!contains({0,2,4,6,8}, 3));
        REQUIRE(has({0,2,4,6,8}, 2));
        REQUIRE(!has({0,2,4,6,8}, 3));
        REQUIRE(in(2, {0,2,4,6,8}));
        REQUIRE(!in(3, {0,2,4,6,8}));
    }
}

