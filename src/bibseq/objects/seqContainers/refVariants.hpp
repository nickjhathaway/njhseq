#pragma once
//
//  baseContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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
#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include "bibseq/alignment/aligner.h"
#include "bibseq/seqToolsUtils/seqToolsUtils.hpp"


namespace bibseq {
class variant {
public:
	variant(const seqInfo & seqBase): seqBase_(seqBase){}
	seqInfo seqBase_;
	//maps keys are position in reference
	std::unordered_map<uint32_t, mismatch> mismatches_;
	std::unordered_map<uint32_t, gap> insertions_;
	std::unordered_map<uint32_t, gap> deletions_;

	void outputInfo(std::ofstream & out, const seqInfo & ref)const{
		//mismatches
		for (const auto & m : mismatches_) {
			out << ref.name_ << "\t" << m.first << "\t" << seqBase_.name_ << "\t"
					<< m.second.refBase << "\t" << m.second.seqBase << "\t" <<seqInfo::getFastqString({m.second.seqQual},
							SangerQualOffset)
					<< std::endl;
		}
		//insertions
		for (const auto & i : insertions_) {
			out << ref.name_ << "\t" << i.first << "\t" << seqBase_.name_ << "\t"
					<< "-" << "\t" << i.second.gapedSequence_ << "\t"
					<< seqInfo::getFastqString(i.second.qualities_,
							SangerQualOffset) << std::endl;
		}
		//deletions
		for (const auto & d : deletions_) {
			out << ref.name_ << "\t" << d.first << "\t" << seqBase_.name_ << "\t"
					<< d.second.gapedSequence_ << "\t" << "_" << "\t"
					<< seqInfo::getFastqString(d.second.qualities_,
							SangerQualOffset) << std::endl;
		}
	}
	/*
	std::string cigarString()const{

		 * R'(\d+)(\w)'

	}


	 * M	BAM_CMATCH	0
I	BAM_CINS	1
D	BAM_CDEL	2
N	BAM_CREF_SKIP	3
S	BAM_CSOFT_CLIP	4
H	BAM_CHARD_CLIP	5
P	BAM_CPAD	6
=	BAM_CEQUAL	7
X	BAM_CDIFF	8

	std::vector<std::pair<uint32_t, uint32_t>> cigarRepresentation()const{

	}*/

};

class refVariants {
public:

	refVariants(const seqInfo & seqBase): seqBase_(seqBase){}
	seqInfo seqBase_;
	std::vector<variant> variants_;

	void addVariant(const seqInfo & var, aligner & alignerObj,
			bool weighHomopolymer) {
		variant varObj(var);
		alignerObj.alignCache(seqBase_, var, false);
		alignerObj.profilePrimerAlignment(seqBase_, var, weighHomopolymer);
		for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
			varObj.mismatches_.emplace(m.second.refBasePos, m.second);
		}
		for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
			if (g.second.ref_) {
				//insertion
				varObj.insertions_.emplace(getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, g.first),
						g.second);
			} else {
				//deletion
				varObj.deletions_.emplace(getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, g.first),
						g.second);
			}
		}
		variants_.emplace_back(varObj);
	}
	std::vector<uint32_t> getVariantSnpLoci()const{
		std::set<uint32_t> loci;
		for(const auto & v : variants_){
			for(const auto & m : v.mismatches_){
				loci.emplace(m.second.refBasePos);
			}
		}
		return std::vector<uint32_t> {loci.begin(), loci.end()};
	}

	std::map<uint32_t, std::vector<char>> getVariantSnpLociMap()const{
		std::map<uint32_t, std::set<char>> lociChars;
		for(const auto & v : variants_){
			for(const auto & m : v.mismatches_){
				lociChars[m.first].emplace(m.second.seqBase);
			}
		}
		std::map<uint32_t, std::vector<char>> ret;
		for(const auto & l : lociChars){
			ret[l.first] =std::vector<char> {l.second.begin(), l.second.end()};
		}
		return ret;
	}

	std::vector<uint32_t> getVariantSnpLoci(VecStr names, uint32_t expand = 0)const{
		std::set<uint32_t> loci;
		for(const auto & v : variants_){
			if(bib::in(v.seqBase_.name_, names)){
				for(const auto & m : v.mismatches_){
					loci.emplace(m.second.refBasePos);
					if(expand > 0){
						for(const auto & e : iter::range<uint32_t>(1,expand + 1)){
							if(e <=m.second.refBasePos){
								loci.emplace(m.second.refBasePos - e);
							}
							if(e + m.second.refBasePos < seqBase_.seq_.size()){
								loci.emplace(m.second.refBasePos + e);
							}
						}
					}
				}
			}
		}
		return std::vector<uint32_t> {loci.begin(), loci.end()};
	}

	std::map<uint32_t, std::vector<char>> getVariantSnpLociMap(VecStr names, uint32_t expand = 0)const{
		std::map<uint32_t, std::set<char>> lociChars;
		for(const auto & v : variants_){
			if(bib::in(v.seqBase_.name_, names)){
				for(const auto & m : v.mismatches_){
					lociChars[m.first].emplace(m.second.seqBase);
					if(expand > 0){
						for(const auto & e : iter::range<uint32_t>(1,expand + 1)){
							if(e <=m.second.refBasePos){
								lociChars[m.first - e].emplace(m.second.seqBase);
							}
							if(e + m.second.refBasePos < seqBase_.seq_.size()){
								lociChars[m.first + e].emplace(m.second.seqBase);
							}
						}
					}
				}
			}
		}
		std::map<uint32_t, std::vector<char>> ret;
		for(const auto & l : lociChars){
			ret[l.first] =std::vector<char> {l.second.begin(), l.second.end()};
		}
		return ret;
	}

	void outPut(std::ofstream & out)const {
		for(const auto & v : variants_){
			v.outputInfo(out, seqBase_);
		}
	}

};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "refVariants.cpp"
#endif
