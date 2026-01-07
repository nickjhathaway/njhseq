#pragma once

//
// Created by Nicholas Hathaway on 1/2/26.
//

#include "njhseq/utils.h"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"
#include "njhseq/objects/Meta/MetaDataInName.hpp"


namespace njhseq{


class MummerOutput {
public:
	//excepts mummer to be run like with the -F option, so expects 4 columns each time regardless of the number of seqs in reference or query
	//e.g. mummer -b -c -F ref.fasta query.fasta

	// mummer -b -c -F -l 200 -maxmatch, maxmatch does all matches regarless of how unique they are in each sequence
	class MummerTargetHits {
	public:
		class MummerRecord {
		public:
			MummerRecord(std::string query, bool reverse, std::string target, uint32_t targetStart, uint32_t queryStart, uint32_t size)
							:  query_(std::move(query)), reverse_(reverse), target_(std::move(target)),
								 targetStart_(targetStart), queryStart_(queryStart), size_(size) {

			}
			MummerRecord(std::string query, bool reverse, std::string line): query_(std::move(query)), reverse_(reverse){

				auto toks = tokenizeString(njh::trim(line), "whitespace");
				if(toks.size()!= 4){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "error, expected 4 items but got: " << toks.size() << "\n";
					throw std::runtime_error{ss.str()};
				}
				target_ = toks[0];
				targetStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				queryStart_ =  njh::StrToNumConverter::stoToNum<uint32_t>(toks[2]);
				size_ =        njh::StrToNumConverter::stoToNum<uint32_t>(toks[3]);
			}
			std::string query_;
			bool reverse_{false};
			std::string target_;
			uint32_t targetStart_{
							std::numeric_limits<uint32_t>::max()}; /**< 1 based positing since this is from the file specs */
			uint32_t queryStart_{
							std::numeric_limits<uint32_t>::max()};/**<1 based positing since this is from the file specs */
			uint32_t size_{std::numeric_limits<uint32_t>::max()};

			[[nodiscard]] Bed3RecordCore genBed3()const{
				Bed3RecordCore ret{target_, targetStart_ -1, targetStart_ -1 + size_};
				MetaDataInName meta;
				meta.addMeta("query", query_);
				if(reverse_){
					meta.addMeta("queryStart", queryStart_ - size_);
					meta.addMeta("queryEnd", queryStart_);
				}else{
					meta.addMeta("queryStart", queryStart_ - 1);
					meta.addMeta("queryEnd", queryStart_ - 1 + size_);
				}
				ret.extraFields_.emplace_back(meta.createMetaName());
				return ret;
			}
			[[nodiscard]] Bed6RecordCore genBed6() const {
				Bed6RecordCore ret{target_, targetStart_ -1, targetStart_ -1 + size_, query_, static_cast<double>(size_), reverse_ ? '-' : '+'};
				MetaDataInName meta;
				meta.addMeta("query", query_);
				if(reverse_){
					meta.addMeta("queryStart", queryStart_ - size_);
					meta.addMeta("queryEnd", queryStart_);
				}else{
					meta.addMeta("queryStart", queryStart_ - 1);
					meta.addMeta("queryEnd", queryStart_ - 1 + size_);
				}
				ret.extraFields_.emplace_back(meta.createMetaName());
				return ret;
			}
		};
		explicit MummerTargetHits(std::string target): target_(std::move(target)){

		};
		std::string target_;
		std::vector<MummerRecord> hits_;
	};

	std::vector<MummerTargetHits> targets_;
};

class NucmerShowCoordsRecord{
public:

	explicit NucmerShowCoordsRecord(const std::string & line){
		auto toks = tokenizeString(line, "\t");

		if(toks.size() < 13){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "should be at least 13 toks, " << toks.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
		refStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[0]);
		refEnd_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
		queryStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[2]);
		queryEnd_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[3]);
		refHitLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[4]);
		queryHitLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);
		perId_ = njh::StrToNumConverter::stoToNum<double>(toks[6]);
		refFullLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[7]);
		queryFullLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[8]);
		refLenCov_ = njh::StrToNumConverter::stoToNum<double>(toks[9]);
		queryLenCov_ = njh::StrToNumConverter::stoToNum<double>(toks[10]);
		refName_ = toks[11];
		queryName_ = toks[12];

	}
	//13 toks
	//show-coords -T -l  -c -H out.delta
	//[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[TAGS]
	// positions are 1 based
	//0 ref start
	//1 ref end
	//2 query start
	//3 query end
	//4 ref hit len
	//5 query hit len
	//6 identity
	//7 full ref length
	//8 full query length
	//9 ref coverage percentage
	//10 query coverage percentage
	//11 reference name
	//12 query name

	uint32_t refStart_{std::numeric_limits<uint32_t>::max()};
	uint32_t refEnd_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryStart_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryEnd_{std::numeric_limits<uint32_t>::max()};
	uint32_t refHitLen_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryHitLen_{std::numeric_limits<uint32_t>::max()};
	double perId_{std::numeric_limits<double>::max()};
	uint32_t refFullLen_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryFullLen_{std::numeric_limits<uint32_t>::max()};
	double refLenCov_{std::numeric_limits<double>::max()};
	double queryLenCov_{std::numeric_limits<double>::max()};
	std::string refName_;
	std::string queryName_;

	[[nodiscard]] bool reverseStrand() const{
		return queryEnd_ < queryStart_;
	}
	[[nodiscard]] Bed6RecordCore genBed6() const{
		MetaDataInName meta;
		meta.addMeta("perID", perId_);
		meta.addMeta("refLenCov", refLenCov_);
		meta.addMeta("queryLenCov", queryLenCov_);
		uint32_t actualStart = reverseStrand() ? queryEnd_ : queryStart_;
		--actualStart;
		uint32_t actualEnd = reverseStrand() ? queryStart_ : queryEnd_;
		meta.addMeta("actualStart", actualStart);
		meta.addMeta("actualEnd", actualEnd);
    meta.addMeta("queryName", queryName_);
		meta.addMeta("queryStart", queryStart_);
		meta.addMeta("queryEnd", queryEnd_);

//		meta.addMeta("queryStart", queryStart_);
//		meta.addMeta("queryEnd", queryEnd_);

		Bed6RecordCore ret(refName_,
											 refStart_ -1 ,
											 refEnd_ ,
//											 reverseStrand() ? refEnd_ -1 : refStart_ -1,
//											 reverseStrand() ? refStart_ : refEnd_,
											 njh::pasteAsStr(queryName_, "-",actualStart , "-", actualEnd),
											 uAbsdiff(refStart_, refEnd_) + 1,
											 reverseStrand() ? '-' : '+');
		ret.extraFields_.emplace_back(meta.createMetaName());
		return ret;
	}

};

}//namespace njhseq
