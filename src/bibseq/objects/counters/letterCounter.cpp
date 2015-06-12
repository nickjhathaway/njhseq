#include "letterCounter.hpp"
#include <iostream>
#include <cmath>

namespace bibseq {

letterCounter::letterCounter(const std::string &seq,
                             const std::vector<uint32_t> &quals) {
	gcContent_ = 0;
  letters['A'] = 0;
  letters['T'] = 0;
  letters['G'] = 0;
  letters['C'] = 0;
  letters['-'] = 0;
  qualities['A'] = 0;
  qualities['T'] = 0;
  qualities['G'] = 0;
  qualities['C'] = 0;
  qualities['-'] = 0;
  for (auto i : seq) {
    ++letters[i];
    qualities[i] += quals[i];
    allQualities[i].push_back(quals[i]);
  }
}

letterCounter::letterCounter(const std::string &seq) {
	gcContent_ = 0;
  for (auto i : seq) {
    ++letters[i];
  }
}



void letterCounter::increaseCountOfBase(const char &base) {
  ++letters[base];
}

void letterCounter::increaseCountByString(const std::string &seq) {
  for (const auto &c : seq) {
    increaseCountOfBase(c);
  }
}
void letterCounter::increaseCountOfBase(const char &base, double cnt) {
	letters[base] += static_cast<uint32_t>(cnt);
}


void letterCounter::increaseCountByString(const std::string &seq, double cnt) {
  for (const auto &c : seq) {
    increaseCountOfBase(c, cnt);
  }
}

void letterCounter::getBest(char &letter) {
  int bestCount = 0;
  for (const auto & let : letters) {
    if (let.second > bestCount) {
      letter = let.first;
      bestCount = let.second;
    }
  }
}
char letterCounter::outputBestLetter() {
  int bestCount = 0;
  char bestBase = ' ';
  for (const auto & let : letters) {
    if (let.second > bestCount) {
      bestBase = let.first;
      bestCount = let.second;
    }
  }
  return bestBase;
}

void letterCounter::getBest(char &letter, uint32_t &quality) {
  int bestCount = 0;
  char bestBase = ' ';
  uint32_t bestQuality = 0;

  if (letters['A'] > bestCount) {
    bestBase = 'A';
    bestQuality = qualities['A'];
    bestCount = letters['A'];
  }
  if (letters['C'] > bestCount) {
    bestBase = 'C';
    bestQuality = qualities['C'];
    bestCount = letters['C'];
  }
  if (letters['G'] > bestCount) {
    bestBase = 'G';
    bestQuality = qualities['G'];
    bestCount = letters['G'];
  }
  if (letters['T'] > bestCount) {
    bestBase = 'T';
    bestQuality = qualities['T'];
    bestCount = letters['T'];
  }

  if (letters['-'] > bestCount) {
    bestBase = '-';
    bestQuality = qualities['-'];
    bestCount = letters['-'];
  }
  letter = bestBase;
  quality = bestQuality;
}

void letterCounter::getBest(char &letter, uint32_t &quality,
                            uint32_t size) {
  char bestBase = ' ';
  uint32_t bestQuality = 0;
  int totalCount = 0;
  totalCount =
      letters['A'] + letters['G'] + letters['C'] + letters['T'] + letters['-'];
  int bestCount = size - totalCount;
  if (letters['A'] > bestCount) {
    bestBase = 'A';
    bestQuality = qualities['A'];
    bestCount = letters['A'];
  }
  if (letters['G'] > bestCount) {
    bestBase = 'G';
    bestQuality = qualities['G'];
    bestCount = letters['G'];
  }
  if (letters['T'] > bestCount) {
    bestBase = 'T';
    bestQuality = qualities['T'];
    bestCount = letters['T'];
  }
  if (letters['C'] > bestCount) {
    bestBase = 'C';
    bestQuality = qualities['C'];
    bestCount = letters['C'];
  }
  if (letters['-'] > bestCount) {
    bestBase = '-';
    bestQuality = qualities['-'];
    bestCount = letters['-'];
  }
  letter = bestBase;
  quality = bestQuality;
}

void letterCounter::outPutInfo(std::ostream &out, bool ifQualities) const {
  out << "letter\tcount\tfraction" << std::endl;
  ;
  double totalCount = getTotalCount();
  if (!ifQualities) {
    for (const auto &let : letters) {
      out << let.first << "\t" << let.second << "\t" << let.second / totalCount
          << std::endl;
    }
  } else {
    for (const auto &let : letters) {
      out << let.first << "\t" << let.second << "\t" << qualities.at(let.first)
          << let.second / totalCount << std::endl;
    }
  }
}

void letterCounter::outPutACGTInfo(std::ostream &out) const {
  if (letters.find('A') == letters.end()) {
    out << 0 << " ";
  } else {
    out << letters.at('A') << " ";
  }
  if (letters.find('C') == letters.end()) {
    out << 0 << " ";
  } else {
    out << letters.at('C') << " ";
  }
  if (letters.find('G') == letters.end()) {
    out << 0 << " ";
  } else {
    out << letters.at('G') << " ";
  }
  if (letters.find('T') == letters.end()) {
    out << 0 << " " << std::endl;
  } else {
    out << letters.at('T') << " " << std::endl;
  }
}

void letterCounter::outPutACGTFractionInfo(std::ostream &out) {
  setFractions();
  if (fractions.find('A') == fractions.end()) {
    out << 0 << " ";
  } else {
    out << fractions['A'] << " ";
  }
  if (fractions.find('C') == fractions.end()) {
    out << 0 << " ";
  } else {
    out << fractions['C'] << " ";
  }
  if (fractions.find('G') == fractions.end()) {
    out << 0 << " ";
  } else {
    out << fractions['G'] << " ";
  }
  if (fractions.find('T') == fractions.end()) {
    out << 0 << " " << std::endl;
  } else {
    out << fractions['T'] << " " << std::endl;
  }
}

void letterCounter::calcGcContent() {
  uint32_t numberOfGC = letters['C'] + letters['G'];
  gcContent_ = (double)numberOfGC / getTotalCount();
}
double letterCounter::computeEntrophy() {
  setFractions();
  double sum = 0;
  for (const auto &frac : fractions) {
    if (0 != frac.second) {
      if (1 == frac.second) {
        return 0;
      }
      sum += frac.second * std::log2(frac.second);
    }
  }
  return (-1 * sum);
}
char letterCounter::getDegenativeBase() const {
  if (letters.at('A') == 0 && letters.at('C') == 0 && letters.at('G') == 0 &&
      letters.at('T') > 0) {
    return 'T';
  } else if (letters.at('A') == 0 && letters.at('C') == 0 &&
             letters.at('G') > 0 && letters.at('T') == 0) {
    return 'G';
  } else if (letters.at('A') == 0 && letters.at('C') > 0 &&
             letters.at('G') == 0 && letters.at('T') == 0) {
    return 'C';
  } else if (letters.at('A') > 0 && letters.at('C') == 0 &&
             letters.at('G') == 0 && letters.at('T') == 0) {
    return 'A';
  } else if (letters.at('A') == 0 && letters.at('C') == 0 &&
             letters.at('G') > 0 && letters.at('T') > 0) {
    return 'K';
  } else if (letters.at('A') == 0 && letters.at('C') > 0 &&
             letters.at('G') == 0 && letters.at('T') > 0) {
    return 'Y';
  } else if (letters.at('A') > 0 && letters.at('C') == 0 &&
             letters.at('G') == 0 && letters.at('T') > 0) {
    return 'W';
  } else if (letters.at('A') == 0 && letters.at('C') > 0 &&
             letters.at('G') > 0 && letters.at('T') == 0) {
    return 'S';
  } else if (letters.at('A') > 0 && letters.at('C') == 0 &&
             letters.at('G') > 0 && letters.at('T') == 0) {
    return 'R';
  } else if (letters.at('A') > 0 && letters.at('C') > 0 &&
             letters.at('G') == 0 && letters.at('T') == 0) {
    return 'M';
  } else if (letters.at('A') == 0 && letters.at('C') > 0 &&
             letters.at('G') > 0 && letters.at('T') > 0) {
    return 'B';
  } else if (letters.at('A') > 0 && letters.at('C') == 0 &&
             letters.at('G') > 0 && letters.at('T') > 0) {
    return 'D';
  } else if (letters.at('A') > 0 && letters.at('C') > 0 &&
             letters.at('G') == 0 && letters.at('T') > 0) {
    return 'H';
  } else if (letters.at('A') > 0 && letters.at('C') > 0 &&
             letters.at('G') > 0 && letters.at('T') == 0) {
    return 'V';
  } else {
    return 'N';
  }
}
int letterCounter::getGcDifference() { return letters['G'] - letters['C']; }
void letterCounter::printDescription(std::ostream &out, bool deep) const {
  out << "letterCounter{" << std::endl;
  out << "letters:" << std::endl;
  out << "\tstd::map<std::string, int>" << std::endl;
  for (const auto &let : letters) {
    out << "\t" << let.first << ":" << let.second << std::endl;
  }
  out << "fractions:" << std::endl;
  out << "\tstd::map<std::string, double>" << std::endl;
  for (const auto &frac : fractions) {
    out << "\t" << frac.first << ":" << frac.second << std::endl;
  }
  out << "qualities:" << std::endl;
  out << "\tstd::map<std::string, uint32_t>" << std::endl;
  for (const auto &qual : qualities) {
    out << "\t" << qual.first << ":" << qual.second << std::endl;
  }
  out << "allQualities:" << std::endl;
  out << "\tstd::map<std::string, std::vector<uint32_t>>" << std::endl;
  for (const auto &qual : allQualities) {
    out << "\t" << qual.first << ":" << qual.second << std::endl;
  }
  out << "gcContent:" << gcContent_ << std::endl << "}" << std::endl;
}

std::multimap<double, char, std::less<double>>
letterCounter::createLikelihoodMaps(bool setFractionFirst) {
  if (setFractionFirst) {
    setFractions();
  }
  std::multimap<double, char, std::less<double>> likelihoods;
  for (const auto &let : fractions) {
    likelihoods.emplace(let.second, let.first);
  }
  return likelihoods;
}

}  // namespace bib
