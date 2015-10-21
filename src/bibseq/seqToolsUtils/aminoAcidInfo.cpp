// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//

/*
 * aminoAcidInfo.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: nickhathaway
 */

#include "bibseq/seqToolsUtils/aminoAcidInfo.hpp"



namespace bibseq {
namespace aminoAcidInfo {

const std::unordered_map<char, aminoAcid> infos::allInfo {
	 {'A', aminoAcid(VecStr{"GCA","GCC","GCG","GCT"}, VecStr{"GCA","GCC","GCG","GCU"}, 4, 'A', "ala","Alanine", "nonpolar", 71.03711, 1.8)},
	 {'C', aminoAcid(VecStr{"TGC","TGT"}, VecStr{"UGC","UGU"}, 2, 'C',"cys","Cysteine", "polar", 103.00919, 2.5)},
	 {'D', aminoAcid(VecStr{"GAC","GAT"}, VecStr{"GAC","GAU"}, 2, 'D', "asp","Aspartic_acid", "acidic",115.02694, -3.5)},
	 {'E', aminoAcid(VecStr{"GAA","GAG"}, VecStr{"GAA","GAG"}, 2, 'E', "glu","Glutamic_acid", "acidic", 129.04259, -3.5)},
	 {'F', aminoAcid(VecStr{"TTC","TTT"}, VecStr{"UUC","UUU"}, 2, 'F', "phe","Phenylalanine", "nonpolar",147.06841, 2.8)},
	 {'G', aminoAcid(VecStr{"GGA","GGC","GGG","GGT"}, VecStr{"GGA","GGC","GGG","GGU"}, 4, 'G', "gly","Glycine", "nonpolar",57.02146, -0.4)},
	 {'H', aminoAcid(VecStr{"CAC","CAT"}, VecStr{"CAC","CAU"}, 2, 'H', "his","Histidine", "basic",137.05891, -3.2)},
	 {'I', aminoAcid(VecStr{"ATA","ATC","ATT"}, VecStr{"AUA","AUC","AUU"}, 3, 'I', "ile","Isoleucine", "nonpolar",113.08406, 4.5)},
	 {'K', aminoAcid(VecStr{"AAA","AAG"}, VecStr{"AAA","AAG"}, 2, 'K', "lys","Lysine", "basic",128.09496, -3.9)},
	 {'L', aminoAcid(VecStr{"CTA","CTC","CTG","CTT","TTA","TTG"}, VecStr{"CUA","CUC","CUG","CUU","UUA","UUG"}, 6, 'L', "leu","Leucine", "nonpolar",113.08406, 3.8)},
	 {'M', aminoAcid(VecStr{"ATG"}, VecStr{"AUG"}, 1, 'M',  "met","Methionine", "nonpolar",131.04049, 1.9)},
	 {'N', aminoAcid(VecStr{"AAC","AAT"}, VecStr{"AAC","AAU"}, 2, 'N', "asn","Asparagine", "polar",114.04293, -3.5)},
	 {'P', aminoAcid(VecStr{"CCA","CCC","CCG","CCT"}, VecStr{"CCA","CCC","CCG","CCU"}, 4, 'P', "pro","Proline", "nonpolar",97.05276, -1.6)},
	 {'Q', aminoAcid(VecStr{"CAA","CAG"}, VecStr{"CAA","CAG"}, 2, 'Q', "gln","Glutamine", "polar",128.05858, -3.5)},
	 {'R', aminoAcid(VecStr{"CGA","CGC","CGG","CGT","AGA","AGG"}, VecStr{"CGA","CGC","CGG","CGU","AGA","AGG"}, 6, 'R', "arg","Arginine", "basic",156.10111, -4.5)},
	 {'S', aminoAcid(VecStr{"TCA","TCC","TCG","TCT","AGC","AGT"}, VecStr{"UCA","UCC","UCG","UCU","AGC","AGU"}, 6, 'S', "ser","Serine", "polar",87.03203, -0.8)},
	 {'T', aminoAcid(VecStr{"ACA","ACC","ACG","ACT"}, VecStr{"ACA","ACC","ACG","ACU"}, 4, 'T', "thr","Threonine", "polar", 101.04768, -0.7)},
	 {'V', aminoAcid(VecStr{"GTA","GTC","GTG","GTT"}, VecStr{"GUA","GUC","GUG","GUU"}, 4, 'V', "val","Valine", "nonpolar", 99.06841, 4.2)},
	 {'W', aminoAcid(VecStr{"TGG"}, VecStr{"UGG"}, 1, 'W', "trp","Tryptophan", "nonpolar",186.07931, -0.9)},
	 {'Y', aminoAcid(VecStr{"TAC","TAT"}, VecStr{"UAC","UAU"}, 2, 'Y', "tyr","Tyrosine", "polar",163.06333, -1.3)},
	 {'*', aminoAcid(VecStr{"TAA","TAG","TGA"}, VecStr{"UAA","UAG","UGA"}, 3, '*', "sto","Stop", "stop", 0, 0)}
};
const std::unordered_map<std::string, std::string> infos::aaClassColorCode {
	{"nonpolar", "ffe641"},
	{"polar", "a4d9b2"},
	{"basic", "adafda"},
	{"acidic", "f6a4c9"},
	{"stop", "a0a0a0"},
};


const std::map<int, std::vector<char>> infos::weightIntToAminoAcid{
    {57, {'G'}},
    {71, {'A'}},
    {87, {'S'}},
    {97, {'P'}},
    {99, {'V'}},
    {101, {'T'}},
    {103, {'C'}},
    {113, {'I', 'L'}},
    {114, {'N'}},
    {115, {'D'}},
    {128, {'K', 'Q'}},
    {129, {'E'}},
    {131, {'M'}},
    {137, {'H'}},
    {147, {'F'}},
    {156, {'R'}},
    {163, {'Y'}},
    {186, {'W'}}};
const std::unordered_map<std::string, char> infos::rnaCodonToAminoACid = {{"GCA", 'A'},
                                                                {"GCC", 'A'},
                                                                {"GCG", 'A'},
                                                                {"GCU", 'A'},
                                                                {"CGA", 'R'},
                                                                {"CGC", 'R'},
                                                                {"CGG", 'R'},
                                                                {"CGU", 'R'},
                                                                {"AGA", 'R'},
                                                                {"AGG", 'R'},
                                                                {"AAC", 'N'},
                                                                {"AAU", 'N'},
                                                                {"GAC", 'D'},
                                                                {"GAU", 'D'},
                                                                {"UGC", 'C'},
                                                                {"UGU", 'C'},
                                                                {"CAA", 'Q'},
                                                                {"CAG", 'Q'},
                                                                {"GAA", 'E'},
                                                                {"GAG", 'E'},
                                                                {"GGA", 'G'},
                                                                {"GGC", 'G'},
                                                                {"GGG", 'G'},
                                                                {"GGU", 'G'},
                                                                {"CAC", 'H'},
                                                                {"CAU", 'H'},
                                                                {"AUA", 'I'},
                                                                {"AUC", 'I'},
                                                                {"AUU", 'I'},
                                                                {"CUA", 'L'},
                                                                {"CUC", 'L'},
                                                                {"CUG", 'L'},
                                                                {"CUU", 'L'},
                                                                {"UUA", 'L'},
                                                                {"UUG", 'L'},
                                                                {"AAA", 'K'},
                                                                {"AAG", 'K'},
                                                                {"AUG", 'M'},
                                                                {"UUC", 'F'},
                                                                {"UUU", 'F'},
                                                                {"CCA", 'P'},
                                                                {"CCC", 'P'},
                                                                {"CCG", 'P'},
                                                                {"CCU", 'P'},
                                                                {"UCA", 'S'},
                                                                {"UCC", 'S'},
                                                                {"UCG", 'S'},
                                                                {"UCU", 'S'},
                                                                {"AGC", 'S'},
                                                                {"AGU", 'S'},
                                                                {"ACA", 'T'},
                                                                {"ACC", 'T'},
                                                                {"ACG", 'T'},
                                                                {"ACU", 'T'},
                                                                {"UGG", 'W'},
                                                                {"UAC", 'Y'},
                                                                {"UAU", 'Y'},
                                                                {"GUA", 'V'},
                                                                {"GUC", 'V'},
                                                                {"GUG", 'V'},
                                                                {"GUU", 'V'},
                                                                {"UAA", '*'},
                                                                {"UAG", '*'},
                                                                {"UGA", '*'}};

const std::unordered_map<std::string, char> infos::dnaCodonToAminoAcid = {{"GCA", 'A'},
                                                                {"GCC", 'A'},
                                                                {"GCG", 'A'},
                                                                {"GCT", 'A'},
                                                                {"CGA", 'R'},
                                                                {"CGC", 'R'},
                                                                {"CGG", 'R'},
                                                                {"CGT", 'R'},
                                                                {"AGA", 'R'},
                                                                {"AGG", 'R'},
                                                                {"AAC", 'N'},
                                                                {"AAT", 'N'},
                                                                {"GAC", 'D'},
                                                                {"GAT", 'D'},
                                                                {"TGC", 'C'},
                                                                {"TGT", 'C'},
                                                                {"CAA", 'Q'},
                                                                {"CAG", 'Q'},
                                                                {"GAA", 'E'},
                                                                {"GAG", 'E'},
                                                                {"GGA", 'G'},
                                                                {"GGC", 'G'},
                                                                {"GGG", 'G'},
                                                                {"GGT", 'G'},
                                                                {"CAC", 'H'},
                                                                {"CAT", 'H'},
                                                                {"ATA", 'I'},
                                                                {"ATC", 'I'},
                                                                {"ATT", 'I'},
                                                                {"CTA", 'L'},
                                                                {"CTC", 'L'},
                                                                {"CTG", 'L'},
                                                                {"CTT", 'L'},
                                                                {"TTA", 'L'},
                                                                {"TTG", 'L'},
                                                                {"AAA", 'K'},
                                                                {"AAG", 'K'},
                                                                {"ATG", 'M'},
                                                                {"TTC", 'F'},
                                                                {"TTT", 'F'},
                                                                {"CCA", 'P'},
                                                                {"CCC", 'P'},
                                                                {"CCG", 'P'},
                                                                {"CCT", 'P'},
                                                                {"TCA", 'S'},
                                                                {"TCC", 'S'},
                                                                {"TCG", 'S'},
                                                                {"TCT", 'S'},
                                                                {"AGC", 'S'},
                                                                {"AGT", 'S'},
                                                                {"ACA", 'T'},
                                                                {"ACC", 'T'},
                                                                {"ACG", 'T'},
                                                                {"ACT", 'T'},
                                                                {"TGG", 'W'},
                                                                {"TAC", 'Y'},
                                                                {"TAT", 'Y'},
                                                                {"GTA", 'V'},
                                                                {"GTC", 'V'},
                                                                {"GTG", 'V'},
                                                                {"GTT", 'V'},
                                                                {"TAA", '*'},
                                                                {"TAG", '*'},
                                                                {"TGA", '*'}};

const std::map<int, VecStr> infos::wieghtToSimilarDoubles = {
    {114, {"GG"}},
    {128, {"GA"}},
    {156, {"GV"}},
    {186, {"GE"}},
    {128, {"AG"}},
    {186, {"AD"}},
    {186, {"SV"}},
    {156, {"VG"}},
    {186, {"VS"}},
    {186, {"DA"}},
    {186, {"EG"}},
    {114, {"GG"}},
    {128, {"GA", "AG"}},
    {156, {"GV", "VG"}},
    {186, {"GE", "AD", "SV", "VS", "DA", "EG"}}};

void codonUsageCounter::increaseCountByString(const std::string &seq, double cnt){
	//check to see if stop is a divisible by three to prevent trying to count a codon less than three
	uint64_t stop = seq.size();
	while (stop % 3 != 0){
		--stop;
	}
	for(const auto & pos : iter::range<uint64_t> (0, stop, 3)){
		counts_[seq.substr(pos, 3)] += cnt;
	}
}




} // aminoAcidInfo
} // bib
