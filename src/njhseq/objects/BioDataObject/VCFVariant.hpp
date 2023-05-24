#pragma once
//
// Created by Nicholas Hathaway on 5/21/23.
//

#include <utility>

#include "njhseq/common.h"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"


namespace njhseq {

class VCFVariant {
public:
	enum class VarType {
		SNP, INSERTION, DELETION
	};


	VCFVariant(GenomicRegion reg, const std::string& ref, const std::string& seq,
						 double freq) :
					region_(std::move(reg)), ref_(ref), variant_(seq), freq_(freq) {
		if (ref.size() == 1 && seq.size() == 1) {
			vtype_ = VCFVariant::VarType::SNP;
		} else if (ref.size() == 1 && seq.size() > 1) {
			vtype_ = VCFVariant::VarType::INSERTION;
		} else if (ref.size() > 1 && seq.size() == 1) {
			vtype_ = VCFVariant::VarType::DELETION;
		}
	}
	GenomicRegion region_;

	std::string ref_;
	std::string variant_;


	VarType vtype_;
	double freq_;
	static std::vector<VCFVariant> readVCFLine(const std::string line);
};


/*
 *
Pf3D7_14_v3	860870	.	G	A	40	PASS	DP=8466;NS=7618;AC=2017;AF=0.238247
Pf3D7_14_v3	861007	.	TAATAATAAT	T	40	PASS	DP=8466;NS=7618;AC=126;AF=0.0148831
Pf3D7_14_v3	861010	.	TAATAAT	T	40	PASS	DP=8466;NS=7618;AC=4159;AF=0.491259
Pf3D7_14_v3	861016	.	T	TAAT	40	PASS	DP=8466;NS=7618;AC=248;AF=0.0292936
 */




} //namespace njhseq

