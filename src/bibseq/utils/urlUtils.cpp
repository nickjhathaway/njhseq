#include "urlUtils.hpp"
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//

namespace bibseq {

std::string urldecode(char const *begin, char const *end) {
	//from cppcms 1.05
	std::string result;
	result.reserve(end - begin);
	for (; begin < end; begin++) {
		char c = *begin;
		switch (c) {
		case '+':
			result += ' ';
			break;
		case '%':
			if (end - begin >= 3 && xdigit(begin[1]) && xdigit(begin[2])) {
				char buf[3] = { begin[1], begin[2], 0 };
				int value;
				sscanf(buf, "%x", &value);
				result += char(value);
				begin += 2;
			}
			break;
		default:
			result += c;
		}
	}
	return result;
}

std::string urldecode(std::string const &s) {
	//from cppcms 1.05
	return urldecode(s.c_str(), s.c_str() + s.size());
}

void urlencode(char const *b, char const *e, std::ostream &out) {
	//from cppcms 1.05
	std::ostream_iterator<char> it(out);
	urlencode_impl(b, e, it);
}

std::string urlencode(const std::string &s) {
	//from cppcms 1.05
	std::string content;
	content.reserve(3 * s.size());
	std::back_insert_iterator<std::string> out(content);
	urlencode_impl(s.c_str(), s.c_str() + s.size(), out);
	return content;
}

}  // namespace bib
