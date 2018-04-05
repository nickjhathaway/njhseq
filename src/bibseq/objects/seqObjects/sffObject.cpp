#include "sffObject.hpp"
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

/*for toJson()
 * 	ret["flowValues"] = bib::json::toJson(flowValues);
	ret["processedFlowValues"] = bib::json::toJson(processedFlowValues);
	ret["numberOfFlows"] = bib::json::toJson(numberOfFlows);
 *
 */

uint32_t fromBase36(std::string base36) {
	try {
		uint32_t num = 0;

		std::unordered_map<char, uint32_t> converts;
		converts['A'] = 0;
		converts['a'] = 0;
		converts['B'] = 1;
		converts['b'] = 1;
		converts['C'] = 2;
		converts['c'] = 2;
		converts['D'] = 3;
		converts['d'] = 3;
		converts['E'] = 4;
		converts['e'] = 4;
		converts['F'] = 5;
		converts['f'] = 5;
		converts['G'] = 6;
		converts['g'] = 6;
		converts['H'] = 7;
		converts['h'] = 7;
		converts['I'] = 8;
		converts['i'] = 8;
		converts['J'] = 9;
		converts['j'] = 9;
		converts['K'] = 10;
		converts['k'] = 10;
		converts['L'] = 11;
		converts['l'] = 11;
		converts['M'] = 12;
		converts['m'] = 12;
		converts['N'] = 13;
		converts['n'] = 13;
		converts['O'] = 14;
		converts['o'] = 14;
		converts['P'] = 15;
		converts['p'] = 15;
		converts['Q'] = 16;
		converts['q'] = 16;
		converts['R'] = 17;
		converts['r'] = 17;
		converts['S'] = 18;
		converts['s'] = 18;
		converts['T'] = 19;
		converts['t'] = 19;
		converts['U'] = 20;
		converts['u'] = 20;
		converts['V'] = 21;
		converts['v'] = 21;
		converts['W'] = 22;
		converts['w'] = 22;
		converts['X'] = 23;
		converts['x'] = 23;
		converts['Y'] = 24;
		converts['y'] = 24;
		converts['Z'] = 25;
		converts['z'] = 25;
		converts['0'] = 26;
		converts['1'] = 27;
		converts['2'] = 28;
		converts['3'] = 29;
		converts['4'] = 30;
		converts['5'] = 31;
		converts['6'] = 32;
		converts['7'] = 33;
		converts['8'] = 34;
		converts['9'] = 35;

		uint32_t i = 0;
		while (i < base36.length()) {
			char c = base36[i];
			num = 36 * num + converts[c];
			++i;
		}

		return num;

	} catch (std::exception& e) {

		throw;
	}
}

void sffObject::sffAdditionalInitialation() {
	baseFlowIndex_.clear();
	flowIndex_.clear();
	selectFlowValues.clear();
	numBases_ = 0;
	clipQualLeft_ = 0;
	clipQualRight_ = 0;
	clipAdapterLeft_ = 0;
	clipAdapterRight_ = 0;

  flowValues.clear();
  numberOfFlows = 0;

}

template<typename T>
std::vector<T> stringToVectorTest(const std::string& strToConvert) {
	std::vector<T> ret;
	std::istringstream in(strToConvert);
	std::copy(std::istream_iterator<T>(in), std::istream_iterator<T>(),
			std::back_inserter(ret));
	return ret;
}

void sffObject::addSffInfo(std::unordered_map<std::string, std::string> &info) {
	// add various info and the sequence
	seqBase_.name_ = info["Name"];
	seqBase_.seq_ = info["Bases"];
	numBases_ = std::stoi(info["# of Bases"]);
	clipAdapterLeft_ = std::stoi(info["Clip Adap Left"]);
	clipAdapterRight_ = std::stoi(info["Clip Adap Right"]);
	clipQualLeft_ = std::stoi(info["Clip Qual Left"]);
	clipQualRight_ = std::stoi(info["Clip Qual Right"]);
	// add qual info
	seqBase_.qual_ = stringToVectorTest<uint32_t>(info["Quality Scores"]);
	// add flow indexes
	flowIndex_ = stringToVectorTest<uint32_t>(info["Flow Indexes"]);
	// add flowgram data
	flowValues = stringToVectorTest<double>(info["Flowgram"]);
	averageErrorRate = getAverageErrorRate();
	numberOfFlows = flowIndex_[clipQualRight_ - 1];
}

void sffObject::setSelectFlowValues() {
	flowValues.clear();
	for (const auto & flowIn : flowIndex_) {
		selectFlowValues.emplace_back(flowValues[flowIn - 1]);
		// std::cout<<(*flowIter)-1<<std::endl;
	}
	// exit(1);
}

void sffObject::setSffClip(size_t left, size_t right) {
	seqBase_.seq_ = seqBase_.seq_.substr(left - 1, right - left + 1);
	seqBase_.qual_.erase(seqBase_.qual_.begin() + right, seqBase_.qual_.end());
	seqBase_.qual_.erase(seqBase_.qual_.begin(),
			seqBase_.qual_.begin() + left - 1);
	averageErrorRate = getAverageErrorRate();
}

std::string sffObject::getFlowIndexesString() const {
	return vectorToString(flowIndex_);
}

std::string sffObject::getFlowValuesString() const {
	return vectorToString(flowValues);
}

void sffObject::outPutSelectFlowvalues(std::ostream &out) const {
	out << ">" << seqBase_.name_ << std::endl;
	out << vectorToString(selectFlowValues) << std::endl;
}

void sffObject::outPutFlowIndexes(std::ostream &out) const {
	out << ">" << seqBase_.name_ << std::endl;
	out << vectorToString(flowIndex_) << std::endl;
}

void sffObject::outPutFlowValuesPryo(std::ostream &flowValuesPryoFile) const {
	flowValuesPryoFile << seqBase_.name_ << " " << flowIndex_[clipQualRight_ - 1]
			<< " " << vectorToString(flowValues) << std::endl;
}

void sffObject::decodeName() {
	try {
		if (seqBase_.name_.length() >= 6) {
			std::string time = seqBase_.name_.substr(0, 6);
			unsigned int timeNum = fromBase36(time);

			int q1 = timeNum / 60;
			int sec = timeNum - 60 * q1;
			int q2 = q1 / 60;
			int minute = q1 - 60 * q2;
			int q3 = q2 / 24;
			int hr = q2 - 24 * q3;
			int q4 = q3 / 32;
			int day = q3 - 32 * q4;
			int q5 = q4 / 13;
			int mon = q4 - 13 * q5;
			int year = 2000 + q5;

			timestamp_ = estd::to_string(year) + "_" + estd::to_string(mon) + "_"
					+ estd::to_string(day) + "_" + estd::to_string(hr) + "_"
					+ estd::to_string(minute) + "_" + estd::to_string(sec);
		}

		if (seqBase_.name_.length() >= 9) {
			region_ = seqBase_.name_.substr(7, 2);

			std::string xyNum = seqBase_.name_.substr(9);
			unsigned int myXy = fromBase36(xyNum);
			int x = myXy >> 12;
			int y = myXy & 4095;

			xy_ = estd::to_string(x) + "_" + estd::to_string(y);
		}
	} catch (std::exception& e) {
		throw;
	}
}

void sffObject::printHeaderSffTxtStyle(std::ofstream& out) {
	out << ">" << seqBase_.name_ << std::endl;
	out << "Run Prefix: " << timestamp_ << std::endl;
	out << "Region #:  " << region_ << std::endl;
	out << "XY Location: " << xy_ << std::endl << std::endl;

	out << "Run Name:  " << std::endl;
	out << "Analysis Name:  " << std::endl;
	out << "Full Path: " << std::endl << std::endl;

	out << "Read Header Len: " << headerLength_ << std::endl;
	out << "Name Length: " << nameLength_ << std::endl;
	out << "# of Bases: " << numBases_ << std::endl;
	out << "Clip Qual Left: " << clipQualLeft_ << std::endl;
	out << "Clip Qual Right: " << clipQualRight_ << std::endl;
	out << "Clip Adap Left: " << clipAdapterLeft_ << std::endl;
	out << "Clip Adap Right: " << clipAdapterRight_ << std::endl << std::endl;
}
void sffObject::printSffTxtSeqData(std::ostream& out) {
	out << "Flowgram: ";
	//for (int i = 0; i < flowgram.size(); ++i) { out << std::setprecision(2) << (flowgram[i]/(float)100) << '\t';  }
	//for (int i = 0; i < flowgram.size(); ++i) { out << std::setprecision(3) << roundDecPlaces((flowgram[i]/(float)100), 2) << '\t';  }
	for (uint32_t i = 0; i < flowValues.size(); ++i) {
		//out << std::setprecision(3) << (flowgram[i]/(float)100) << '\t';
		out << std::setprecision(3) << flowValues[i] << '\t';
	}
	out << std::endl << "Flow Indexes: ";
	int sum = 0;
	for (uint32_t i = 0; i < baseFlowIndex_.size(); ++i) {
		sum += baseFlowIndex_[i];
		out << sum << '\t';
	}

	//make the bases you want to clip lower case and the bases you want to keep upper case
	int endValue = clipQualRight_;
	if (endValue == 0) {
		endValue = seqBase_.seq_.length();
	}
	for (int i = 0; i < (clipQualLeft_ - 1); ++i) {
		seqBase_.seq_[i] = tolower(seqBase_.seq_[i]);
	}
	for (int i = (clipQualLeft_ - 1); i < (endValue); ++i) {
		seqBase_.seq_[i] = toupper(seqBase_.seq_[i]);
	}
	for (uint32_t i = (endValue); i < seqBase_.seq_.length(); ++i) {
		seqBase_.seq_[i] = tolower(seqBase_.seq_[i]);
	}
	out << std::endl << "Bases: " << seqBase_.seq_ << std::endl
			<< "Quality Scores: ";
	for (uint32_t i = 0; i < seqBase_.qual_.size(); ++i) {
		out << seqBase_.qual_[i] << '\t';
	}
	out << std::endl << std::endl;
}
bool sffObject::sanityCheck() {
	bool okay = true;
	std::string message = "[WARNING]: Your sff file may be corrupted! Sequence: "
			+ seqBase_.name_ + "\n";

	if (clipQualLeft_ > seqBase_.seq_.length()) {
		okay = false;
		message += "Clip Qual Left = " + estd::to_string(clipQualLeft_)
				+ ", but we only read " + estd::to_string(seqBase_.seq_.length())
				+ " bases.\n";
	}
	if (clipQualRight_ > seqBase_.seq_.length()) {
		okay = false;
		message += "Clip Qual Right = " + estd::to_string(clipQualRight_)
				+ ", but we only read " + estd::to_string(seqBase_.seq_.length())
				+ " bases.\n";
	}
	if (clipQualLeft_ > seqBase_.qual_.size()) {
		okay = false;
		message += "Clip Qual Left = " + estd::to_string(clipQualLeft_)
				+ ", but we only read " + estd::to_string(seqBase_.qual_.size())
				+ " quality scores.\n";
	}
	if (clipQualRight_ > seqBase_.qual_.size()) {
		okay = false;
		message += "Clip Qual Right = " + estd::to_string(clipQualRight_)
				+ ", but we only read " + estd::to_string(seqBase_.qual_.size())
				+ " quality scores.\n";
	}

	if (okay == false) {
		//m->mothurOut(message); m->mothurOutEndLine();
	}

	return okay;
}
int sffBinaryHeader::printSffTxtStyle(std::ostream & out) {
	out << "Common Header: \nMagic Number: " << magicNumber << std::endl;
	out << "Version: " << version << std::endl;
	out << "Index Offset: " << indexOffset << std::endl;
	out << "Index Length: " << indexLength << std::endl;
	out << "Number of Reads: " << numReads << std::endl;
	out << "Header Length: " << headerLength << std::endl;
	out << "Key Length: " << keyLength << std::endl;
	out << "Number of Flows: " << numFlowsPerRead << std::endl;
	out << "Format Code: " << flogramFormatCode << std::endl;
	out << "Flow Chars: " << flowChars << std::endl;
	out << "Key Sequence: " << keySequence << std::endl << std::endl;
	return 0;
}

void sffBinaryHeader::printDescription(std::ostream & out, bool deep) const {
	out << "CommonHeader{" << std::endl;
	out << "magicNumber: " << magicNumber << std::endl;
	out << "version: " << version << std::endl;
	out << "indexOffset: " << indexOffset << std::endl;
	out << "indexLength: " << indexLength << std::endl;
	out << "numReads: " << numReads << std::endl;
	out << "headerLength: " << headerLength << std::endl;
	out << "keyLength: " << keyLength << std::endl;
	out << "numFlowsPerRead: " << numFlowsPerRead << std::endl;
	out << "flogramFormatCode: " << flogramFormatCode << std::endl;
	out << "flowChars: " << flowChars << std::endl;
	out << "keySequence: " << keySequence << std::endl;
	out << "}" << std::endl;
}


void sffObject::outPutFlows(std::ostream& flowsFile) const {
  flowsFile << ">" << seqBase_.name_ << std::endl;
  flowsFile << vectorToString(flowValues) << std::endl;
}

void sffObject::outPutFlowsRaw(std::ostream& out) const {
  out << ">" << seqBase_.name_ << std::endl;
  out << numberOfFlows << " " << vectorToString(flowValues) << std::endl;
}

void sffObject::outPutPyroData(std::ostream& pyroNoiseFile) const {
  pyroNoiseFile << seqBase_.name_ << " " << numberOfFlows << " "
                << vectorToString(flowValues) << std::endl;
}

uint32_t sffObject::getNumberOfBasesFromFlow(
    const std::vector<double>& inFlows) {
  uint32_t numOfBases = 0;
  for (auto fIter = inFlows.begin(); fIter != inFlows.end(); ++fIter) {
    if (fIter == inFlows.begin() && *fIter > 20) {
      continue;
    }
    if ((*fIter) > 0.7) {
      // non noisy signal increase bases count by nearest integer
      numOfBases += std::round(*fIter);
    }
  }
  return numOfBases;
}
int sffObject::clipToNinetyPercentOfFlows(size_t cutOff) {
  processedFlowValues = flowValues;
  if (processedFlowValues.size() > cutOff) {
    processedFlowValues.erase(processedFlowValues.begin() + cutOff,
                              processedFlowValues.end());
  }
  std::cout << "Clipped to : " << processedFlowValues.size() << std::endl;
  int clipTo = getNumberOfBasesFromFlow(flowValues);
  // setClip(0, clipTo);
  //exit(
  return clipTo;
}

size_t sffObject::findFirstNoisyFlow(const std::vector<double>& inFlows) {
  int noise = 0;
  int signal = 0;
  int basesSinceSignal = 0;
  size_t pos = 0;
  for (std::vector<double>::const_iterator fIter = inFlows.begin();
       fIter != inFlows.end(); ++fIter) {

    if (*fIter > 0.5) {
      if (*fIter < 0.7) {
        ++noise;
      } else {
        ++signal;
        basesSinceSignal = 0;
      }
    } else {
      ++basesSinceSignal;
    }
    if (basesSinceSignal >= 4 || noise >= 1) {
      break;
    }
    ++pos;
  }
  return pos;
}

bool sffObject::flowNoiseProcess(size_t cutoff) {
  size_t firstNoiseReadPos = findFirstNoisyFlow(flowValues);
  if (firstNoiseReadPos < 360) {
    return false;
  }
  if (firstNoiseReadPos < cutoff) {
    cutoff = firstNoiseReadPos;
  }
  processedFlowValues = flowValues;
  processedFlowValues.erase(processedFlowValues.begin() + cutoff,
                            processedFlowValues.end());

  auto numberOfBases = getNumberOfBasesFromFlow(processedFlowValues);
  if (numberOfBases < seqBase_.seq_.size()) {
    setClip(0, numberOfBases);
  }
  return true;
}



}  // namespace bib
