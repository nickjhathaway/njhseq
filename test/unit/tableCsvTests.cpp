#include <catch_amalgamated.hpp>
#include "njhseq/objects/dataContainers/tables/table.hpp"

using namespace njhseq;

TEST_CASE("joinRow quotes CSV fields that tokenizeLine unquotes", "[table][csv]") {
	VecStr row{"S1", "Kenya, Kisumu", "say \"hi\"", "", "plain"};
	auto line = table::joinRow(row, ",");
	REQUIRE(line == "S1,\"Kenya, Kisumu\",\"say \"\"hi\"\"\",,plain");
	REQUIRE(table::tokenizeLine(line, ",") == row);
}

TEST_CASE("a quoted CSV row read and written again keeps its columns", "[table][csv]") {
	std::string line = "S1,\"Kenya, Kisumu\",3";
	auto row = table::tokenizeLine(line, ",");
	REQUIRE(row.size() == 3);
	REQUIRE(row[1] == "Kenya, Kisumu");
	REQUIRE(table::tokenizeLine(table::joinRow(row, ","), ",") == row);
}

TEST_CASE("joinRow leaves tab delimited rows unquoted", "[table][csv]") {
	VecStr row{"a,b", "c\"d", "e"};
	REQUIRE(table::joinRow(row, "\t") == "a,b\tc\"d\te");
}
