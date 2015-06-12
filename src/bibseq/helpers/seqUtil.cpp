#include "seqUtil.hpp"
#include "bibseq/simulation/mutator.hpp"
#include "bibseq/IO/fileUtils.hpp"

namespace bibseq {



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
        ss << "Error in reverseComplement.  Unknown base character = "
                  << thisBase << ".  Exiting." << std::endl;
        throw std::runtime_error{ss.str()};  // Stop right away.
        break;
    }
    outSeq[numChar - 1 - i] = outBase;
  }
  return outSeq;
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
  size_t numChar = seq.size();
  std::string outSeq("");
  outSeq.resize(numChar / 3);  // numChar is exactly divisible by 3.
  // In below, cB is currentBase.
  char cB[3], newBase;
  std::string cBstring = "   ";
  for (size_t i = start; i < numChar; i += 3) {
    cB[0] = static_cast<char>(toupper(seq[i]));
    cB[1] = static_cast<char>(toupper(seq[i + 1]));
    cB[2] = static_cast<char>(toupper(seq[i + 2]));
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

    if (forceStartM && (i == start)) {
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





table seqUtil::readBarcodes(const std::string &idFileName,
                            const std::string &fileDelim, int &barcodeSize,
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
                           const std::string &fileDelim, bool forceRead) {
	table inTab(idFileName, fileDelim);
	inTab.removeEmpty(false);
  table ans;
  ans.hasHeader_ = true;
  ans.columnNames_ = {"geneName", "forwardPrimer", "reversePrimer"};
  bool readingGene = forceRead;
  bool readingBarcode = false;
  for (const auto &fIter : inTab.content_) {
    if (stringToLowerReturn(fIter[0]) == "gene") {
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




void seqUtil::processQualityWindowString(const std::string &qualityWindowString,
                                         int &qualityWindowLength,
                                         int &qualityWindowStep,
                                         int &qualityWindowThres) {
  VecStr toks = tokenizeString(qualityWindowString, ",");
  if (toks.size() != 3) {
    std::cout << "Error qualityWindow must be given with three values separatd "
                 "by two commas like \'50,5,20\' not \'" << qualityWindowString
              << "\'" << std::endl;
  } else {
    qualityWindowLength = atoi(toks[0].c_str());
    qualityWindowStep = atoi(toks[1].c_str());
    qualityWindowThres = atoi(toks[2].c_str());
  }
}

bool seqUtil::checkQualityWindow(int windowSize, int minimumAverageQaul,
                                 int stepSize,
                                 const std::vector<uint32_t> &quality) {
  bool pass = true;
  uint32_t currentPos = 0;
  while ((windowSize + currentPos) < quality.size() && pass) {
    uint32_t sum = 0;
    for (const auto & qPos : iter::range(currentPos, currentPos + windowSize)) {
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
                                      const std::vector<uint32_t> &quality) {
	bool pass = true;
  uint32_t currentPos = 0;
  while ((windowSize + currentPos) < quality.size() && pass) {
    uint32_t sum = 0;
    for (const auto & qPos : iter::range(currentPos, currentPos + windowSize)) {
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


void seqUtil::removeLowerCase(std::string &sequence,
                              std::vector<uint32_t> &quality) {
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
std::pair<std::string, std::vector<uint32_t>> seqUtil::removeLowerCaseReturn(
    std::string sequence, std::vector<uint32_t> quality) {
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
  return {sequence, quality};
}

std::string seqUtil::removeGapsReturn(const std::string &seq) {
  return removeCharReturn(seq, '-');
}
void seqUtil::removeGaps(std::string &seq) {
  removeChar(seq, '-');
  return;
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
    const std::vector<uint32_t> &qual, const std::vector<uint32_t> &positions) {
  std::vector<uint32_t> ans;
  ans.reserve(positions.size());
  for (const auto &pos : iter::range(positions.size() - 1)) {
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

std::map<std::string, std::unordered_map<std::string, std::vector<double>>>
seqUtil::getCountsForModel(
    std::map<std::string, std::map<double, uint32_t>> counts,
    std::map<std::string, std::map<double, uint32_t>> mismatchCounts) {
  std::map<std::string, std::unordered_map<std::string, std::vector<double>>>
      ans;
  for (const auto &count : counts) {
    if (bib::in(count.first, {"mean", "median", "base", "min"})) {
      ans[count.first] = getCountsForSpecificModel(
          count.second, mismatchCounts.at(count.first));
    }
  }
  return ans;
}
std::unordered_map<std::string, std::vector<double>>
seqUtil::getCountsForSpecificModel(std::map<double, uint32_t> counts,
                                   std::map<double, uint32_t> mismatchCounts) {
  std::unordered_map<std::string, std::vector<double>> currentCounts;
  std::ofstream currentBigCountFile;
  currentBigCountFile << "error\tqual\tweight" << std::endl;
  for (const auto &subCount : counts) {
    if (bib::has(mismatchCounts, subCount.first)) {
      currentCounts["error"].emplace_back(0);
      currentCounts["qual"].emplace_back(subCount.first);
      currentCounts["weight"]
          .emplace_back(subCount.second - mismatchCounts.at(subCount.first));
      currentCounts["error"].emplace_back(1);
      currentCounts["qual"].emplace_back(subCount.first);
      currentCounts["weight"].emplace_back(mismatchCounts.at(subCount.first));
      // currentBigCountFile << 0 << "\t" << subCount.first << "\t" <<
      // subCount.second - mismatchCounts[count.first].at(subCount.first) <<
      // std::endl;
      // currentBigCountFile << 1 << "\t" << subCount.first << "\t" <<
      // mismatchCounts[count.first].at(subCount.first) << std::endl;
    } else {
      // no mismatches for this quality
      currentCounts["error"].emplace_back(0);
      currentCounts["qual"].emplace_back(subCount.first);
      currentCounts["weight"].emplace_back(subCount.second);
      currentCounts["error"].emplace_back(1);
      currentCounts["qual"].emplace_back(subCount.first);
      currentCounts["weight"].emplace_back(0);
    }
  }
  return currentCounts;
}

std::map<std::string, std::unordered_map<double, double>>
seqUtil::getTrueErrorRate(
    std::map<std::string, std::map<double, uint32_t>> counts,
    std::map<std::string, std::map<double, uint32_t>> mismatchCounts) {
  std::map<std::string, std::unordered_map<double, double>> ans;
  for (const auto &count : counts) {
    if (bib::in(count.first, {"mean", "median", "base", "min"})) {
      ans[count.first] = getTrueErrorRateSpecific(
          count.second, mismatchCounts.at(count.first));
    }
  }
  return ans;
}
std::unordered_map<double, double> seqUtil::getTrueErrorRateSpecific(
    std::map<double, uint32_t> counts,
    std::map<double, uint32_t> mismatchCounts) {
  std::unordered_map<double, double> currentCounts;
  for (const auto &subCount : counts) {
    if (bib::has(mismatchCounts, subCount.first)) {
      currentCounts[subCount.first] =
          (double)mismatchCounts.at(subCount.first) / subCount.second;
    } else {
      // no mismatches for this quality
      currentCounts[subCount.first] = 0;
    }
  }
  return currentCounts;
}


void seqUtil::rstripRead(std::string & str,
		std::vector<uint32_t> & qual, char c){
	uint32_t pos = str.size();
	while (pos != 0 && str[pos - 1] == c){
		--pos;
	}
	if(pos != str.size()){
		qual.erase(qual.begin() + pos, qual.end());
		str.erase(str.begin() + pos, str.end());
	}
}



\
}  // namespace bibseq
