#include "urlUtils.hpp"
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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

size_t WriteCallback(char *contents, size_t size, size_t nmemb,
                     std::ostream *stream) {
  stream->write(contents, size * nmemb);
  return size * nmemb;
}

std::string GetURL(const std::string url) {
  CURL *curl;
  CURLcode res;
  std::stringstream ss;
  curl = curl_easy_init();
  if (curl) {
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &ss);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "libcurl-agent/1.0");
    res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
      fprintf(stderr, "curl_easy_perform() failed: %s\n",
              curl_easy_strerror(res));
    }
    curl_easy_cleanup(curl);
  }
  return ss.str();
}
void GetURLStream(const std::string url, std::ostream & out){
  CURL *curl;
  CURLcode res;
  //std::stringstream ss;
  curl = curl_easy_init();
  if (curl) {
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &out);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "libcurl-agent/1.0");
    res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
      fprintf(stderr, "curl_easy_perform() failed: %s\n",
              curl_easy_strerror(res));
    }
    curl_easy_cleanup(curl);
  }
  return;
}


}  // namespace bib
