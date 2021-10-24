//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "seqUtil.hpp"
#include "njhseq/IO/fileUtils.hpp"

namespace njhseq {

double seqUtil::calculateWeightOfProteinDouble(const std::string &protein) {
  double weight = 0;
  for (const auto &c : protein) {
    weight += aminoAcidInfo::infos::allInfo.at(c).weight_;
  }
  return weight;
}
int seqUtil::calculateWeightOfProteinInt(const std::string &protein) {
  int weight = 0;
  for (const auto &c : protein) {
    weight += std::round(aminoAcidInfo::infos::allInfo.at(c).weight_);
  }
  return weight;
}

std::string seqUtil::reverseComplement(const std::string &seq,
                                       const std::string &seqType) {
  size_t numChar = seq.size();
  std::string outSeq("");
  outSeq.resize(numChar);
  char thisBase, outBase;
  for (size_t i = 0; i < numChar; ++i) {
    thisBase = seq[i];
    switch (thisBase) {
      case 'a':
        if (seqType == "DNA") {
          outBase = 't';
        } else  // RNA
        {
          outBase = 'u';
        }
        break;
      case 'A':
        if (seqType == "DNA") {
          outBase = 'T';
        } else  // RNA
        {
          outBase = 'U';
        }
        break;
      case 'c':
        outBase = 'g';
        break;
      case 'C':
        outBase = 'G';
        break;
      case 'g':
        outBase = 'c';
        break;
      case 'G':
        outBase = 'C';
        break;
      case 't':
        outBase = 'a';
        break;
      case 'T':
        outBase = 'A';
        break;
      case 'u':
        outBase = 'a';
        break;
      case 'U':
        outBase = 'A';
        break;
      case 'r':
        outBase = 'y';
        break;  // R = A or G
      case 'R':
        outBase = 'Y';
        break;
      case 'y':
        outBase = 'r';
        break;  // Y = C or T/U
      case 'Y':
        outBase = 'R';
        break;
      case 's':
        outBase = 's';
        break;  // S = C or G
      case 'S':
        outBase = 'S';
        break;
      case 'w':
        outBase = 'w';
        break;  // W = A or T/U
      case 'W':
        outBase = 'W';
        break;
      case 'k':
        outBase = 'm';
        break;  // K = G or T/U
      case 'K':
        outBase = 'M';
        break;
      case 'm':
        outBase = 'k';
        break;  // M = A or C
      case 'M':
        outBase = 'K';
        break;
      case 'b':
        outBase = 'v';
        break;  // B = C, G, or T/U
      case 'B':
        outBase = 'V';
        break;
      case 'd':
        outBase = 'h';
        break;  // D = A, G, or T/U
      case 'D':
        outBase = 'H';
        break;
      case 'h':
        outBase = 'd';
        break;  // H = A, C, or T/U
      case 'H':
        outBase = 'D';
        break;
      case 'v':
        outBase = 'b';
        break;  // V = A, C, or G
      case 'V':
        outBase = 'B';
        break;
      case 'n':
        outBase = 'n';
        break;  // N = A, C, G, or T/U
      case 'N':
        outBase = 'N';
        break;
      case '.':
        outBase = '.';
        break;  // . and - denote a gap
      case '-':
        outBase = '-';
        break;
      case '*':
        outBase = '*';
        break;  // * is a padding symbol from
                // the SAM/CIGAR format.
      default:  // Should never get here.
        // Unrecognized base.
      	std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << " : Error, Unknown base character = "
                  << thisBase << ".  Exiting." << "\n";
        ss << "Whole string: " << seq << "\n";
        throw std::runtime_error{ss.str()};  // Stop right away.
        break;
    }
    outSeq[numChar - 1 - i] = outBase;
  }
  return outSeq;
}

bool seqUtil::reversePalindrome(const std::string &seq,
                                const std::string &seqType) {
  if (seq.size() % 2 != 0 || seq.size() < 2) {
    return false;
  } else {
    std::string reverse = reverseComplement(seq, seqType);
    if (reverse == seq) {
      return true;
    } else {
      return false;
    }
  }
}

std::map<size_t, size_t> seqUtil::findReversePalindromes(
    const std::string &seq, const std::string &seqType, size_t lowerSizeLimit,
    size_t upperSizeLimit) {
  std::map<size_t, size_t> locationsLengths;
  for (size_t size = lowerSizeLimit; size <= upperSizeLimit; size += 2) {
    size_t pos = 0;
    while (pos + size <= seq.size()) {
      if (reversePalindrome(seq.substr(pos, size), seqType)) {
        locationsLengths.insert(std::make_pair(pos, size));
      }
      ++pos;
    }
  }
  return locationsLengths;
}

void seqUtil::printOutReversePalindromes(
    const std::string &seq, const std::string &seqType, std::ostream &out,
    bool multipleAtSameSite, size_t lowerSizeLimit, size_t upperSizeLimit) {
  std::map<size_t, size_t> locationsLengths = seqUtil::findReversePalindromes(
      seq, seqType, lowerSizeLimit, upperSizeLimit);
  size_t previousPosition = std::string::npos;
  for (const auto &iter : locationsLengths) {
    if (!multipleAtSameSite && iter.first == previousPosition) {

    } else {
      out << iter.first + 1 << " " << iter.second << std::endl;
    }
    previousPosition = iter.first;
  }
}

// forceStartM means if the string starts AUG, GUG, UUG, AUU, or CUG, then
// put M (start codon) at the beginning.
std::string seqUtil::convertOneCodon(const std::string &codon) {
  if (codon.size() != 3) {
    std::cout << "codon needs to be three bases long" << std::endl;
    return "";
  }
  size_t numChar = codon.size();
  std::string outSeq("");
  outSeq.resize(numChar / 3);  // numChar is exactly divisible by 3.
  // In below, cB is currentBase.
  char cB[3], newBase;
  std::string cBstring = "   ";
  int i = 0;
  cB[0] = static_cast<char>(toupper(codon[i]));
  cB[1] = static_cast<char>(toupper(codon[i + 1]));
  cB[2] = static_cast<char>(toupper(codon[i + 2]));
  if (cB[0] == 'T') {
    cB[0] = 'U';
  }
  if (cB[1] == 'T') {
    cB[1] = 'U';
  }
  if (cB[2] == 'T') {
    cB[2] = 'U';
  }
  newBase = '0';  // We'll designate * as Stop,
  // and 1 as error.
  if (cB[0] == 'U') {
    if (cB[1] == 'U') {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'F';
      } else {
        newBase = 'L';
      }
    } else if (cB[1] == 'C') {
      newBase = 'S';
    } else if (cB[1] == 'A') {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'Y';
      } else {
        newBase = '*';
      }
    } else  // cB[1] == 'G'
    {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'C';
      } else if (cB[2] == 'A') {
        newBase = '*';
      } else {
        newBase = 'W';
      }
    }
  } else if (cB[0] == 'C') {
    if (cB[1] == 'U') {
      newBase = 'L';
    } else if (cB[1] == 'C') {
      newBase = 'P';
    } else if (cB[1] == 'A') {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'H';
      } else if (cB[2] == 'A') {
        newBase = 'Q';
      } else {
        newBase = 'Q';
      }
    } else  // cB[1] == 'G'
    {
      newBase = 'R';
    }
  } else if (cB[0] == 'A') {
    if (cB[1] == 'U') {
      if (cB[2] == 'G') {
        newBase = 'M';
      } else {
        newBase = 'I';
      }
    } else if (cB[1] == 'C') {
      newBase = 'T';
    } else if (cB[1] == 'A') {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'N';
      } else {
        newBase = 'K';
      }
    } else  // cB[1] == 'G'
    {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'S';
      } else {
        newBase = 'R';
      }
    }
  } else  // cb[0] == 'G'
  {
    if (cB[1] == 'U') {
      newBase = 'V';
    } else if (cB[1] == 'C') {
      newBase = 'A';
    } else if (cB[1] == 'A') {
      if (cB[2] == 'U' || cB[2] == 'C') {
        newBase = 'D';
      } else {
        newBase = 'E';
      }
    } else  // cB[1] == 'G'
    {
      newBase = 'G';
    }
  }

  outSeq[0] = newBase;
  return outSeq;
}

std::string seqUtil::convertToProtein(const std::string &seq, size_t start,
                                      bool forceStartM) {

  size_t numChar = seq.size() > start ? seq.size() - start : seq.size();

  std::string outSeq("");
  if(numChar < 3){
  		return outSeq;
  }
  outSeq.resize(numChar / 3);  // numChar is exactly divisible by 3.
  // In below, cB is currentBase.
  char cB[3];
	char newBase;

//  std::cout << "seq: " << seq << std::endl;

  for (size_t i = start; i < seq.size() - 2; i += 3) {
//  	std::cout << "i: " << i << ", outSeq: " << outSeq << std::endl;
    cB[0] = static_cast<char>(toupper(seq[i]));
    cB[1] = static_cast<char>(toupper(seq[i + 1]));
    cB[2] = static_cast<char>(toupper(seq[i + 2]));
//    std::cout << "\tcB:" << cB << std::endl;
    if (cB[0] == 'T') {
      cB[0] = 'U';
    }
    if (cB[1] == 'T') {
      cB[1] = 'U';
    }
    if (cB[2] == 'T') {
      cB[2] = 'U';
    }
    newBase = '0';  // We'll designate * as Stop,
    // and 1 as error.
    if(cB[0] == '-' || cB[1] == '-' || cB[2] == '-'){
    	newBase = 'X';
    } else if (cB[0] == 'U') {
      if (cB[1] == 'U') {
        if (cB[2] == 'U' || cB[2] == 'C') {
          newBase = 'F';
        } else {
          newBase = 'L';
        }
      } else if (cB[1] == 'C') {
        newBase = 'S';
      } else if (cB[1] == 'A') {
        if (cB[2] == 'U' || cB[2] == 'C') {
          newBase = 'Y';
        } else {
          newBase = '*';
        }
      } else  // cB[1] == 'G'
      {
        if (cB[2] == 'U' || cB[2] == 'C') {
          newBase = 'C';
        } else if (cB[2] == 'A') {
          newBase = '*';
        } else {
          newBase = 'W';
        }
      }
    } else if (cB[0] == 'C') {
      if (cB[1] == 'U') {
        newBase = 'L';
      } else if (cB[1] == 'C') {
        newBase = 'P';
      } else if (cB[1] == 'A') {
        if (cB[2] == 'U' || cB[2] == 'C') {
          newBase = 'H';
        } else if (cB[2] == 'A') {
          newBase = 'Q';
        } else {
          newBase = 'Q';
        }
      } else  // cB[1] == 'G'
      {
        newBase = 'R';
      }
    } else if (cB[0] == 'A') {
      if (cB[1] == 'U') {
        if (cB[2] == 'G') {
          newBase = 'M';
        } else {
          newBase = 'I';
        }
      } else if (cB[1] == 'C') {
        newBase = 'T';
      } else if (cB[1] == 'A') {
        if (cB[2] == 'U' || cB[2] == 'C') {
          newBase = 'N';
        } else {
          newBase = 'K';
        }
      } else  // cB[1] == 'G'
      {
        if (cB[2] == 'U' || cB[2] == 'C') {
          newBase = 'S';
        } else {
          newBase = 'R';
        }
      }
		} else  // cb[0] == 'G'
		{
			if (cB[1] == 'U') {
				newBase = 'V';
			} else if (cB[1] == 'C') {
				newBase = 'A';
			} else if (cB[1] == 'A') {
				if (cB[2] == 'U' || cB[2] == 'C') {
					newBase = 'D';
				} else {
					newBase = 'E';
				}
			} else  // cB[1] == 'G'
			{
				newBase = 'G';
			}
		}

    if (forceStartM && (i == start)) {
    	std::string cBstring = "   ";
      cBstring[0] = cB[0];
      cBstring[1] = cB[1];
      cBstring[2] = cB[2];
      // Below are all the known start codons.
      if ((cBstring == "AUG") || (cBstring == "GUG") || (cBstring == "UUG") ||
          (cBstring == "AUU") || (cBstring == "CUG")) {
        newBase = 'M';
      }
    }
    outSeq[(i - start) / 3] = newBase;
  }

//	std::cout << "final outSeq: " << outSeq << std::endl;

  return outSeq;
}

void seqUtil::transcribe(std::string &theWord) {

  for (uint32_t i = 0; i < theWord.length(); i++) {
    if (theWord[i] == 'T') {
      theWord[i] = 'A';
    } else if (theWord[i] == 'A') {
      theWord[i] = 'T';
    } else if (theWord[i] == 'G') {
      theWord[i] = 'C';
    } else if (theWord[i] == 'C') {
      theWord[i] = 'G';
    } else if (theWord[i] == 't') {
      theWord[i] = 'a';
    } else if (theWord[i] == 'a') {
      theWord[i] = 't';
    } else if (theWord[i] == 'g') {
      theWord[i] = 'c';
    } else if (theWord[i] == 'c') {
      theWord[i] = 'g';
    }
  }
}

void seqUtil::changeToRNA(std::string &theWord) {
  for (uint32_t i = 0; i < theWord.length(); i++) {
    if (theWord[i] == 'T') {
      theWord[i] = 'U';
    } else if (theWord[i] == 't') {
      theWord[i] = 'u';
    }
  }
}
std::string seqUtil::changeToRNAOut(const std::string &theWord) {
  std::string out = theWord;
  changeToRNA(out);
  return out;
}

void seqUtil::transcribeToRNA(std::string &theWord) {
  transcribe(theWord);
  changeToRNA(theWord);
  return;
}

void seqUtil::convertToProteinFromcDNA(std::string &theWord, size_t start,
                                       bool forceStartM) {
  transcribeToRNA(theWord);
  theWord = convertToProtein(theWord, start, forceStartM);
}
std::string seqUtil::convertToProteinFromcDNAReturn(const std::string &theWord,
                                                    size_t start,
                                                    bool forceStartM) {
  std::string hold = theWord;
  transcribeToRNA(hold);
  return convertToProtein(hold, start, forceStartM);
}

char seqUtil::translatePhred(int qual) {
  if (qual <= 93) {
    return ((char)(qual + 33));
  } else {
    return ((char)(93 + 33));
  }
}

int seqUtil::translatePhred(char qual) {
  int quality = qual;
  return (quality - 33);
}

int seqUtil::computeHammingDistance(const std::string &firstSeq,
                                    const std::string &secondSeq) {
  int distance = 0;
  if (firstSeq.size() != secondSeq.size()) {
    std::cout << "strings need to be the same length" << std::endl;
    return 0;
  }
  if (firstSeq.find("-") != std::string::npos ||
      secondSeq.find("-") != std::string::npos) {
    std::cout << "strings can't contain gaps" << std::endl;
    return 0;
  }
  for (uint32_t i = 0; i < firstSeq.size(); ++i) {
    if (firstSeq[i] != secondSeq[i]) {
      ++distance;
    }
  }
  return distance;
}
/*
std::string seqUtil::consensusFromVectorOfStrings(const VecStr &strings,
                                                  std::ostream &out) {
  std::vector<letterCounter> counts;
  std::vector<size_t> sizes;
  for (VecStr::const_iterator strIter = strings.begin();
       strIter != strings.end(); ++strIter) {
    sizes.push_back(strIter->size());
  }
  std::vector<char> dnaLetters={'A','C','G','T'};
  size_t maxSize = vectorMaximum(sizes);
  for (size_t i = 0; i < maxSize; ++i) {
    counts.push_back(letterCounter(dnaLetters));
  }
  for (VecStr::const_iterator strIter = strings.begin();
       strIter != strings.end(); ++strIter) {
    for (uint32_t i = 0; i < strIter->size(); ++i) {
      counts[i].increaseCountOfBase(strIter->substr(i, 1));
    }
  }
  std::stringstream consensus;
  for (std::vector<letterCounter>::iterator countIter = counts.begin();
       countIter != counts.end(); ++countIter) {
    consensus << countIter->outputBestLetter();
    countIter->setFractions();
  }
  out << consensus.str() << std::endl;

  out << "A:";
  for (std::vector<letterCounter>::iterator countIter = counts.begin();
       countIter != counts.end(); ++countIter) {
    if (countIter->letters.find("A") == countIter->letters.end()) {
      out << " " << 0;
    } else {
      out << " " << countIter->letters["A"];
    }
  }
  out << std::endl;

  out << "C:";
  for (std::vector<letterCounter>::iterator countIter = counts.begin();
       countIter != counts.end(); ++countIter) {
    if (countIter->letters.find("C") == countIter->letters.end()) {
      out << " " << 0;
    } else {
      out << " " << countIter->letters["C"];
    }
  }
  out << std::endl;

  out << "G:";
  for (std::vector<letterCounter>::iterator countIter = counts.begin();
       countIter != counts.end(); ++countIter) {
    if (countIter->letters.find("G") == countIter->letters.end()) {
      out << " " << 0;
    } else {
      out << " " << countIter->letters["G"];
    }
  }
  out << std::endl;

  out << "T:";
  for (std::vector<letterCounter>::iterator countIter = counts.begin();
       countIter != counts.end(); ++countIter) {
    if (countIter->letters.find("T") == countIter->letters.end()) {
      out << " " << 0;
    } else {
      out << " " << countIter->letters["T"];
    }
  }
  out << std::endl;
  std::map<std::string, std::vector<double>> fractions;
  for(const auto & counter : counts){
    for(const auto & letterCount : counter.fractions){
      fractions[letterCount.first].push_back(letterCount.second);
    }
  }
  VecStr dnaLettersStrings={"A", "C", "G", "T"};
  for(const auto & letter : dnaLettersStrings){
    std::cout<<letter<<": ";
    printVector(fractions[letter]);
  }
  double sumOfEntrophy=0;
  std::cout<<"E:";
  for(auto & counter : counts){
    double currentEntrophy=counter.computeEntrophy();
    sumOfEntrophy+=currentEntrophy;
    std::cout<<" "<<currentEntrophy;
  }
  std::cout<<std::endl;
  std::cout<<"Score: "<<sumOfEntrophy<<std::endl;;
  return consensus.str();
}*/

// read scoring matrix
std::unordered_map<char, std::unordered_map<char, int>>
seqUtil::readScoringMatrix(const std::string &fileName) {
	std::ifstream inFile(fileName);
  return(readScoringMatrix(inFile));
}
std::unordered_map<char, std::unordered_map<char, int>> seqUtil::readScoringMatrix(std::istream& in){
	table inTab(in, " ", false);
  std::unordered_map<char, std::unordered_map<char, int>> matchMap;
  std::vector<char> baseOrder;
  int count = 0;
  for (const auto &fileIter : inTab.content_) {
    if (count == 0) {
      for (const auto &strIter : fileIter) {
        baseOrder.push_back((strIter)[0]);
      }
    } else {
      char currentBase = ' ';
      int basePos = 0;
      int secondCount = 0;
      for (const auto &strIter : fileIter) {
        if (secondCount == 0) {
          currentBase = (strIter)[0];
        } else {
          matchMap[baseOrder[basePos]][currentBase] = std::stoi(strIter);
          basePos++;
        }
        ++secondCount;
      }
    }
    ++count;
  }
  return matchMap;
}


VecStr seqUtil::getAllProteins(const std::string &seq) {
  VecStr ans;
  if (seq.size() < 6) {
    std::cout << "sequence needs to be at least 6 bases long" << std::endl;
    return ans;
  }

  for (int i = 0; i < 3; ++i) {
    ans.push_back(convertToProtein(seq.substr(i)));
    std::string hold = reverseComplement(seq, "DNA");
    ans.push_back(convertToProtein(hold.substr(i)));
  }
  return ans;
}

VecStr seqUtil::findAllOpenFrames(const VecStr &proteins) {
  VecStr ans;
  for (VecStr::const_iterator siter = proteins.begin(); siter != proteins.end();
       ++siter) {
    if (siter->find("M") == std::string::npos ||
        siter->find("*") == std::string::npos) {
      continue;
    }
    std::vector<size_t> mPositions = findOccurences(*siter, "M");
    std::vector<size_t> starPositions = findOccurences(*siter, "*");
    std::sort(starPositions.begin(), starPositions.end());
    std::sort(mPositions.begin(), mPositions.end());
    for (std::vector<size_t>::iterator miter = mPositions.begin();
         miter != mPositions.end(); ++miter) {
      for (std::vector<size_t>::iterator starIter = starPositions.begin();
           starIter != starPositions.end(); ++starIter) {
        if (*starIter > *miter) {
          if (vectorContains(ans,
                             siter->substr(*miter, (*starIter - *miter)))) {

          } else {
            ans.push_back(siter->substr(*miter, (*starIter - *miter)));
          }
          break;
        }
      }
    }
  }
  return ans;
}

VecStr seqUtil::findOpenFramesFromSeq(const std::string &seq) {
  VecStr proteins = getAllProteins(seq);
  return findAllOpenFrames(proteins);
}



table seqUtil::readBarcodes(const std::string &idFileName,
                            const std::string &fileDelim,
														int &barcodeSize,
                            bool forceRead) {
	table inTab(idFileName, fileDelim);
	inTab.removeEmpty(false);
  table ans;
  ans.columnNames_ = {"id", "barcode"};
  ans.hasHeader_ = true;
  bool readingGene = false;
  bool readingBarcode = forceRead;
  for (const auto &row : inTab.content_) {
    if (stringToLowerReturn(row[0]) == "gene") {
      readingGene = true;
      readingBarcode = false;
      continue;
    } else if (stringToLowerReturn(row[0]) == "id") {
      readingGene = false;
      readingBarcode = true;
      continue;
    } else {
      if (readingGene) {

      } else if (readingBarcode) {
        barcodeSize = row[1].size();
        ans.content_.emplace_back(row);
      }
    }
  }
  ans.setColNamePositions();
  return ans;
}

table seqUtil::readPrimers(const std::string &idFileName,
                           const std::string &fileDelim,
													 bool forceRead) {
	table inTab(idFileName, fileDelim);
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	//inTab.outPutContentOrganized(std::cout);
	inTab.removeEmpty(false);
  table ans;
  ans.hasHeader_ = true;
  ans.columnNames_ = {"geneName", "forwardPrimer", "reversePrimer"};
  bool readingGene = forceRead;
  bool readingBarcode = false;
  for (const auto &fIter : inTab.content_) {
    if ("gene" == stringToLowerReturn(fIter[0]) ||
    		"target" == stringToLowerReturn(fIter[0])) {
      readingGene = true;
      readingBarcode = false;
      continue;
    } else if (stringToLowerReturn(fIter[0]) == "id") {
      readingGene = false;
      readingBarcode = true;
      continue;
    } else {
      if (readingGene) {
        ans.content_.push_back(fIter);
      } else if (readingBarcode) {

      }
    }
  }
  ans.setColNamePositions();

  return ans;
}


VecStr seqUtil::getBarcodesInOrderTheyAppear(const std::string &fileName) {
	table inTab(fileName);
  VecStr ans;
  bool readingGene = false;
  bool readingBarcodes = true;
  for (const auto &fIter : inTab.content_) {
    if (fIter[0] == "gene") {
      readingGene = true;
      readingBarcodes = false;
      continue;
    }
    if (fIter[0] == "id") {
      readingGene = false;
      readingBarcodes = true;
      continue;
    }
    if (readingBarcodes && !readingGene) {
      ans.push_back(fIter[0]);
    }
  }

  return ans;
}

std::map<std::string, int> seqUtil::makeDNAKmerCompMap(int kLength) {
  std::string dna = "A C G T";
  VecStr dnaStrings = stringToVector<std::string>(dna);
  std::vector<VecStr> strVectors = permuteVector(dnaStrings, kLength);
  std::map<std::string, int> ans;
  for (const auto &iter : strVectors) {
    ans.insert(std::make_pair(vectorToString(iter, ""), 0));
  }
  return ans;
}
std::map<std::string, kmer> seqUtil::makeDNAKmerMap(int kLength) {
  std::map<std::string, int> stringMap = makeDNAKmerCompMap(kLength);
  std::map<std::string, kmer> kmerMapAns;
  for (const auto &strMap : stringMap) {
    kmerMapAns.emplace(strMap.first, kmer(strMap.first, -1));
  }
  return kmerMapAns;
}

double seqUtil::getAverageErrorRate(const std::vector<int> &qual) {
  std::vector<int>::const_iterator itr;
  double sum = 0;
  for (itr = qual.begin(); itr < qual.end(); itr++) {
    sum += pow(10.0, -(*itr) / 10.0);
  }
  return sum / qual.size();
}

VecStr seqUtil::findLongestShardedMotif(VecStr dnaStrings) {
  stringSorter strSorter = stringSorter();
  strSorter.sortStrByLength(dnaStrings);
  std::string shortString = dnaStrings[0];
  VecStr longestVec;
  for (size_t i = shortString.length(); i != 1; --i) {
    bool doesContain = true;
    bool found = false;
    size_t pos = 0;
    while (pos + i <= shortString.length()) {
      doesContain = true;
      for (VecStr::iterator sIter = dnaStrings.begin();
           sIter != dnaStrings.end(); ++sIter) {
        if (sIter->find(shortString.substr(pos, i)) == std::string::npos) {
          doesContain = false;
          break;
        }
      }
      if (doesContain) {
        found = true;
        longestVec.push_back(shortString.substr(pos, i));
      }
      pos++;
    }
    if (found) {
      break;
    }
  }
  return longestVec;
}
std::string seqUtil::findLongestSharedSubString(VecStr dnaStrings) {
  stringSorter strSorter = stringSorter();
  strSorter.sortStrByLength(dnaStrings);
  std::string shortString = dnaStrings[0];
  std::string longestString;
  for (size_t i = shortString.length(); i != 1; --i) {
    bool doesContain = true;
    bool found = false;
    size_t pos = 0;
    while (pos + i <= shortString.length()) {
      doesContain = true;
      for (VecStr::iterator sIter = dnaStrings.begin();
           sIter != dnaStrings.end(); ++sIter) {
        if (sIter->find(shortString.substr(pos, i)) == std::string::npos) {
          doesContain = false;
          break;
        }
      }
      if (doesContain) {
        found = true;
        longestString = shortString.substr(pos, i);
      }
      pos++;
    }
    if (found) {
      break;
    }
  }
  return longestString;
}

std::string seqUtil::getSeqFromFlow(const std::vector<double> &flows,
                                    const std::string &flowSeq) {
  /**
   * @brief Convert a flow gram to sequence.
   *
   * @param flows The flows to be translated
   * @param flowSeq The flowOrder, defaults to TACG
   * @return The translated sequence
   */
  std::string ans = "";
  int flowNumber = 0;
  uint32_t roundsWithOutSignal = 0;
  for (const auto &flowIter : flows) {
    int signal = flowIter + 0.5;
    char base = flowSeq[flowNumber % flowSeq.size()];
    if (signal > 0) {
      roundsWithOutSignal = 0;
      for (int i = 0; i < signal; ++i) {
        ans.push_back(base);
      }
    } else {
      roundsWithOutSignal++;
      if (roundsWithOutSignal == flowSeq.size()) {
        ans.push_back('N');
      }
    }
    flowNumber++;
  }
  return ans;
}

void seqUtil::processQualityWindowString(const std::string &qualityWindowString,
		uint32_t &qualityWindowLength, uint32_t &qualityWindowStep,
		uint8_t &qualityWindowThres) {
	VecStr toks = tokenizeString(qualityWindowString, ",");
	if (toks.size() != 3) {
		std::cout << "Error qualityWindow must be given with three values separatd "
				"by two commas like \'50,5,20\' not \'" << qualityWindowString << "\'"
				<< std::endl;
	} else {
		qualityWindowLength = atoi(toks[0].c_str());
		qualityWindowStep = atoi(toks[1].c_str());
		qualityWindowThres = atoi(toks[2].c_str());
	}
}

bool seqUtil::checkQualityWindow(int windowSize, int minimumAverageQaul,
                                 int stepSize,
                                 const std::vector<uint8_t> &quality) {
  bool pass = true;
  uint32_t currentPos = 0;
  while ((windowSize + currentPos) < quality.size() && pass) {
    uint32_t sum = 0;
    for (const auto qPos : iter::range(currentPos, currentPos + windowSize)) {
      sum += quality[qPos];
    }
    if ((static_cast<double>(sum) / windowSize) < minimumAverageQaul) {
      pass = false;
      break;
    } else {
      currentPos += stepSize;
    }
  }
  return pass;
}
size_t seqUtil::checkQualityWindowPos(int windowSize, int minimumAverageQaul,
                                      int stepSize,
                                      const std::vector<uint8_t> &quality) {
	bool pass = true;
  uint32_t currentPos = 0;
  while ((windowSize + currentPos) < quality.size() && pass) {
    uint32_t sum = 0;
    for (const auto qPos : iter::range(currentPos, currentPos + windowSize)) {
      sum += quality[qPos];
    }
    if ((static_cast<double>(sum) / windowSize) < minimumAverageQaul) {
      pass = false;
      break;
    } else {
      currentPos += stepSize;
    }
  }
  if (pass) {
    return quality.size() - 1;
  } else {
    return currentPos;
  }
}

bool seqUtil::doesSequenceContainDegenerativeBase(const std::string &seq) {
  if (seq.find("N") != std::string::npos ||
      seq.find("R") != std::string::npos ||
      seq.find("S") != std::string::npos ||
      seq.find("Y") != std::string::npos ||
      seq.find("W") != std::string::npos ||
      seq.find("K") != std::string::npos ||
      seq.find("M") != std::string::npos) {
    return true;
  } else {
    return false;
  }
}

std::string seqUtil::removeIntronsThenTranslate(std::string dna,
                                                const VecStr &introns) {
  for (const auto &sIter : introns) {
    dna = njh::replaceString(dna, sIter, "");
  }
  return convertToProtein(dna, 0, true);
}

double seqUtil::calculateTransitionTransversionRatio(const std::string &seq1,
                                                     const std::string &seq2) {
  int transitions = 0;
  int transversions = 0;
  if (seq1.size() != seq2.size()) {
    std::cout << "Sequences need to be same length" << std::endl;
    return -1;
  } else {
    for (size_t i = 0; i < seq1.size(); ++i) {

      if (seq1[i] != seq2[i]) {
        if (toupper(seq1[i]) == 'G' || toupper(seq1[i]) == 'A') {
          if (toupper(seq2[i]) == 'G' || toupper(seq2[i]) == 'A') {
            ++transitions;
          } else if (toupper(seq2[i]) == 'C' || toupper(seq2[i]) == 'T') {
            ++transversions;
          } else {
            std::cout << "Unrecognized base " << seq2[i] << std::endl;
          }
        } else if (toupper(seq1[i]) == 'C' || toupper(seq1[i]) == 'T') {
          if (toupper(seq2[i]) == 'G' || toupper(seq2[i]) == 'A') {
            ++transversions;
          } else if (toupper(seq2[i]) == 'C' || toupper(seq2[i]) == 'T') {
            ++transitions;
          } else {
            std::cout << "Unrecognized base " << seq1[i] << std::endl;
          }
        } else {
          std::cout << "Unrecognized base " << seq1[i] << std::endl;
        }
      }
    }
  }
  if (transversions == 0) {
    return (double)std::string::npos;
  } else {
    return (double)transitions / transversions;
  }
}

bool seqUtil::isHomopolyer(const std::string &seq) {
  if (seq.size() == 1) {
    return true;
  } else {
    for (size_t i = 1; i < seq.size(); ++i) {
      if (seq[i] != seq[i - 1]) {
        return false;
      }
    }
  }
  return true;
}

bool seqUtil::checkTwoEqualSeqs(const std::string &seq1,
                                const std::string &seq2,
                                int allowableMismatches) {
  int currentMismatches = 0;
  for (auto i : iter::range(seq1.size())) {
    if (seq1[i] != seq2[i]) {
      ++currentMismatches;
      if (currentMismatches > allowableMismatches) {
        return false;
      }
    }
  }
  return true;
}





uint64_t seqUtil::getNumberOfPossibleDNAStrandsFromProtein(
    const std::string &protein) {
  uint64_t ans = 1;
  for (const auto &c : protein) {
    ans *= aminoAcidInfo::infos::allInfo.at(c).numCodons_;
  }
  if (protein.back() != '*') {
    ans *= 3;
  }
  return ans;
}

void seqUtil::removeLowerCase(std::string &sequence,
                              std::vector<uint8_t> &quality) {
  for (uint32_t i = 0; i < sequence.size(); i++) {
    if (islower(sequence[i])) {
      sequence.erase(sequence.begin() + i);
      quality.erase(quality.begin() + (i));
      if (i == sequence.size()) {
        break;
      }
      i--;
    }
  }
}
std::pair<std::string, std::vector<uint8_t>> seqUtil::removeLowerCaseReturn(
    std::string sequence, std::vector<uint8_t> quality) {
  for (uint32_t i = 0; i < sequence.size(); ++i) {
    if (islower(sequence[i])) {
      sequence.erase(sequence.begin() + i);
      quality.erase(quality.begin() + (i));
      if (i == sequence.size()) {
        break;
      }
      --i;
    }
  }
  return {sequence, quality};
}



std::string seqUtil::removeGapsReturn(const std::string &seq) {
  return removeCharReturn(seq, '-');
}
void seqUtil::removeGaps(std::string &seq) {
  removeChar(seq, '-');
  return;
}

std::unordered_map<uint64_t, std::string> seqUtil::findMinimumHammingDistance(
    const std::string &seq, const std::string &subSeq, int kLength) {
  std::unordered_map<uint64_t, std::string> ans;
  int minimumHammingDistance = (int)subSeq.size() + 1;
  for (auto i : iter::range(seq.size() - kLength + 1)) {
    std::string currentSubString = seq.substr(i, kLength);
    int currentHammingDistance =
        computeHammingDistance(currentSubString, subSeq);
    // std::cout<<"hamDis: "<<currentHammingDistance<<" seq:
    // "<<currentSubString<<" subSeq: "<<subSeq<<std::endl;
    if (currentHammingDistance == minimumHammingDistance) {
      ans[i] = currentSubString;
    } else if (currentHammingDistance < minimumHammingDistance) {
      minimumHammingDistance = currentHammingDistance;
      ans.clear();
      ans[i] = currentSubString;
    }
  }
  return ans;
}
std::string seqUtil::createDegenerativeString(const VecStr &dnaString) {
  // size check, should be same size
  uint64_t firstSize = dnaString.front().size();
  for (const auto &dna : dnaString) {
    if (dna.size() != firstSize) {
      std::cout << "All strings should be same size" << std::endl;
      return "";
    }
  }

  std::vector<charCounter> counters(firstSize);
  for (auto i : iter::range(firstSize)) {
    for (const auto &dna : dnaString) {
      counters[i].increaseCountOfBase(dna[i]);
    }
  }
  std::string ans;
  for (const auto &count : counters) {
    ans.push_back(count.getDegenativeBase());
  }
  return ans;
}
std::string seqUtil::genMotifStrAccountDegenBase(const std::string & primer){
	std::unordered_map<char, std::vector<char>> degenBaseExapnd = {
			{ 'N', std::vector<char>{'A', 'C','G', 'T'}},
			{ 'K', std::vector<char>{'G', 'T'}},
			{ 'Y', std::vector<char>{'C', 'T'}},
			{ 'W', std::vector<char>{'A', 'T'}},
			{ 'S', std::vector<char>{'C', 'G'}},
			{ 'R', std::vector<char>{'A', 'G'}},
			{ 'M', std::vector<char>{'A', 'C'}},
			{ 'B', std::vector<char>{'C', 'G', 'T'}},
			{ 'D', std::vector<char>{'A', 'G', 'T'}},
			{ 'H', std::vector<char>{'A', 'C', 'T'}},
			{ 'V', std::vector<char>{'A', 'C', 'G'}}
	};
	std::string primerMotif = "";
	for(const auto base : primer){
		if(njh::in(base, degenBaseExapnd)){
			std::string degen = "[";
			for(const auto & d : degenBaseExapnd.at(base)){
				degen.push_back(d);
			}
			degen.push_back(']');
			primerMotif.append(degen);
		}else{
			primerMotif.push_back(base);
		}
	}
	return primerMotif;
}





std::vector<uint32_t> seqUtil::getQualPositions(const std::string &consensus,
                                                const std::string &compare) {
  std::vector<uint32_t> ans;
  uint32_t pos = 0;
  for (auto i : iter::range(consensus.size())) {
    if (consensus[i] == '-') {
      ++pos;
    } else if (compare[i] == '-') {
      ans.push_back(pos);
    } else {
      ans.push_back(pos);
      ++pos;
    }
  }
  return ans;
}
std::vector<uint32_t> seqUtil::rearrangeQuals(
    const std::vector<uint8_t> &qual, const std::vector<uint32_t> &positions) {
  std::vector<uint32_t> ans;
  ans.reserve(positions.size());
  for (const auto pos : iter::range(positions.size() - 1)) {
    if (positions[pos] != positions[pos + 1]) {
      ans.push_back(qual[positions[pos]]);
    } else {
      if (positions[pos] == qual.size()) {
        ans.push_back(qual.back());
      } else if (positions[pos] != 0) {
        ans.push_back((qual[positions[pos] - 1] + qual[positions[pos]]) / 2);
      } else if (positions[pos] == 0) {
        ans.push_back(qual[0]);
      } else {
        // i dont' think this should happen...
      }
    }
  }
  if (positions.back() != qual.size()) {
    ans.push_back(qual[positions.back()]);
  } else {
    ans.push_back(qual.back());
  }
  return ans;
}
uint32_t seqUtil::countMismatchesInAlignment(const std::string &ref,
                                             const std::string &compare,
                                             const char &ignore) {
  uint32_t count = 0;
  auto mis = std::make_pair(ref.begin(), compare.begin());
  while (mis.first != ref.end()) {
    mis = std::mismatch(mis.first, ref.end(), mis.second);
    if (mis.first != ref.end()) {
      if (*mis.first != ignore && *mis.second != ignore) {
        ++count;
      }
      ++mis.first;
      ++mis.second;
    }
  }
  return count;
}




void seqUtil::rstripRead(std::string & str,
		std::vector<uint8_t> & qual, char c){
	uint32_t pos = str.size();
	while (pos != 0 && str[pos - 1] == c){
		--pos;
	}
	if(pos != str.size()){
		qual.erase(qual.begin() + pos, qual.end());
		str.erase(str.begin() + pos, str.end());
	}
}




}  // namespace njhseq
