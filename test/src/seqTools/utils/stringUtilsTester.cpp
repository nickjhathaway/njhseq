#include <catch.hpp>

#include "../src/bibseq/utils/stringUtils.hpp"
using namespace bibseq;

TEST_CASE( "String utils", "[stringUtils]" ) {
    SECTION("getStringFromSubstrings"){
        std::vector<size_t> pos = {0,2,4,6,8};
        auto ret = getStringFromSubstrings("abcdefghijklmnopqrstuvwxyz",
                                           pos,
                                           1);
        REQUIRE(ret == "acegi");
    }

    SECTION("trimStringAtFirstOccurence"){
        std::string str = "aaaaaaaaaabbbbbbbbbbb";
        trimStringAtFirstOccurence(str, "ab");
        REQUIRE(str == "aaaaaaaaa");
    }

    SECTION("trimStringsAtFirstOccurence"){
        std::vector<std::string> v = { "ab",
                                       "aaab",
                                       "aaaaaaaaaabbbbbbbbbbb"};
        trimStringsAtFirstOccurence(v, "ab");
        std::vector<std::string> r = {"", "aa", "aaaaaaaaa"};
        REQUIRE(v == r);
    }

    SECTION("containsSubString"){
        REQUIRE(!containsSubString("", "bc"));
        REQUIRE(containsSubString("aaaaaaabc", "bc"));
    }

    SECTION("stringStartsWith"){
        REQUIRE(!beginsWith("", "aaaaa"));
        REQUIRE(beginsWith("aaaaabbbbb", "aaaaa"));
    }

    SECTION("repeatString"){
        REQUIRE(repeatString("a", 0) == std::string(""));
        REQUIRE(repeatString("a", 4) == std::string("aaaa"));
    }

    SECTION("trimEndWhiteSpace"){
        std::string in = "  a  ";
        std::string out = "a";
        auto c = in;
        trimEndWhiteSpace(c);
        REQUIRE(c == out);
        REQUIRE(trimEndWhiteSpaceReturn(in) == out);
        REQUIRE(trimEndWhiteSpaceReturn("        ") == "");

        std::string whiteSpaceMiddle("ab       cd");
        REQUIRE(trimEndWhiteSpaceReturn(whiteSpaceMiddle) == whiteSpaceMiddle);

        std::string whiteSpaceMiddle2(" ab       cd ");
        REQUIRE(trimEndWhiteSpaceReturn(whiteSpaceMiddle2) == whiteSpaceMiddle);
    }
}
