#pragma once

//
// Created by Nicholas Hathaway on 1/30/24.
//

#include "njhseq/objects/BioDataObject.h"


namespace njhseq {
class VCFOutput {

	/**
	 * Class was defined trying to conform to the standards laid out by https://samtools.github.io/hts-specs/VCFv4.2.pdf, accessed on 2023-12-17
	 **/

	/**
*Header line syntax
The header line names the 8 fixed, mandatory columns. These columns are as follows: 1. #CHROM
2. POS
3. ID
4. REF
5. ALT
6. QUAL
7. FILTER
8. INFO
If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number of sample IDs. Duplicate sample IDs are not allowed. The header line is tab-delimited.

Info field: - additional information: (String, no whitespace, semicolons, or equals-signs permitted; commas are permitted only as delimiters for lists of values)

	*/
public:
	std::string vcfFormatVersion_{"VCFv4.0"};
	struct InfoEntry {
		InfoEntry() = default;
		InfoEntry(		std::string id,
		std::string number,
		std::string type,
		std::string description): id_(std::move(id)),
			number_(std::move(number)),
			type_(std::move(type)),
			description_(std::move(description)) {

		}
		std::string id_;
		/**
		 * \brief description of the number of elements in this entry
		 *
		 * per VCR 4.2 standards
		 * • If the field has one value per alternate allele then this value should be ‘A’.
       • If the field has one value for each possible allele (including the reference), then this value should be ‘R’.
       • If the field has one value for each possible genotype (more relevant to the FORMAT tags) then this value should be ‘G’.
       • If the number of possible values varies, is unknown, or is unbounded, then this value should be ‘.’.
		 */
		std::string number_;
		std::string type_;/**< the data type, should be of Integer, Float, Flag, Character, and String */
		std::string description_;

		std::string source_;/**< optional */
		std::string version_;/**< optional */

		bool operator==(const InfoEntry& other) const {
			return this->id_ == other.id_ &&
			       this->number_ == other.number_ &&
			       this->type_ == other.type_ &&
			       this->description_ == other.description_ &&
			       this->source_ == other.source_ &&
			       this->version_ == other.version_;
		}

		[[nodiscard]] Json::Value toJson() const ;
	};

	struct FormatEntry {
		FormatEntry() = default;
		FormatEntry(		std::string id,
		std::string number,
		std::string type,
		std::string description): id_(std::move(id)),
			number_(std::move(number)),
			type_(std::move(type)),
			description_(std::move(description)) {

		}
		std::string id_;
		/**
		 * \brief description of the number of elements in this entry
		 *
		 * per VCR 4.2 standards
		 * • If the field has one value per alternate allele then this value should be ‘A’.
			 • If the field has one value for each possible allele (including the reference), then this value should be ‘R’.
			 • If the field has one value for each possible genotype (more relevant to the FORMAT tags) then this value should be ‘G’.
			 • If the number of possible values varies, is unknown, or is unbounded, then this value should be ‘.’.
		 */
		std::string number_;
		std::string type_;/**< the data type, should be of Integer, Float, Character, and String */
		std::string description_;
		bool operator==(const FormatEntry& other) const {
			return this->id_ == other.id_ &&
						 this->number_ == other.number_ &&
						 this->type_ == other.type_ &&
						 this->description_ == other.description_;
		}
		[[nodiscard]] Json::Value toJson() const ;
	};

	struct FilterEntry {
		FilterEntry() = default;

		FilterEntry(std::string id,
		            std::string description): id_(std::move(id)),
		                                      description_(std::move(description)) {
		}

		std::string id_;
		std::string description_;
		bool operator==(const FilterEntry& other) const {
			return this->id_ == other.id_ &&
						 this->description_ == other.description_ ;
		}
		[[nodiscard]] Json::Value toJson() const ;
	};


	struct ContigEntry {
		//contig=<ID=22,length=49691432,assembly=B36,md5=2041e6a0c914b48dd537922cca63acb8,species="Homo sapiens">
		ContigEntry() = default;
		ContigEntry(std::string id, uint32_t length): id_(std::move(id)), length_(length) {

		}
		std::string id_;/**<chromosome name */
		uint32_t length_{std::numeric_limits<uint32_t>::max()};/**<length of chromosome */

		std::string assembly_;/**<associated assembly name, can be left blank */
		std::string md5_;/**<the md5 sum to check the sequence, can be left blank */
		std::string species_;/**<sepcies name, can be left blank */

		std::unordered_map<std::string, std::string> otherKeysValues_;/**<anything that isn't ID, length, assembly, md5, or species */
		bool operator==(const ContigEntry& other) const {
			return this->id_ == other.id_ &&
						 this->length_ == other.length_ &&
						 this->assembly_ == other.assembly_ &&
						 this->md5_ == other.md5_ &&
						 this->species_ == other.species_ &&
						 this->otherKeysValues_ == other.otherKeysValues_;
		}
		[[nodiscard]] Json::Value toJson() const ;
	};

	std::map<std::string,FormatEntry> formatEntries_;
	std::map<std::string,InfoEntry> infoEntries_;
	std::vector<FilterEntry> filterEntries_;
	std::map<std::string,ContigEntry> contigEntries_;

	class VCFRecord {
	public:
		VCFRecord() = default;
		std::string chrom_;/**<chromosome name */
		uint32_t pos_{std::numeric_limits<uint32_t>::max()};/**<position, 1-based to conform with VCF standards */
		std::string id_{"."};/**<id for this variant */
		std::string ref_;/**<reference sequence for this variant */
		VecStr alts_;/**<all the alternative alleles */
		uint32_t qual_{40};/**<the quality for this variant */
		std::string filter_{"PASS"};/**<whether or not this variant passes QC checks */

		MetaDataInName info_;/**<the info entry for this variant */
		std::map<std::string, MetaDataInName> sampleFormatInfos_;/**<the sample info for this variant */

		[[nodiscard]] uint32_t getNumberOfAlleles() const ;

		[[nodiscard]] Json::Value toJson() const ;

		/**
		 * \brief Generate a genomic region for this record, will only have the genomic region info, none of the sample info if loaded
		 * \return a genomic region that covers this variant, will be zero-based positioning
		 */
		GenomicRegion genRegion() const;
	};

	std::vector<VCFRecord> records_;/**< the variant records*/

	std::multimap<std::string, MetaDataInName> otherHeaderMetaFields_;/**< meta fields that aren't currently handled, meta fields being lines with ##[META]=<>*/
	std::multimap<std::string, std::string> otherHeaderValuePairs_;/**< meta fields that are just ##key=value*/

	VecStr headerNonSampleFields_;/**<the header that begins with #CHROM */
	VecStr samples_;/**< the samples in the order they appear in the VCF file */

	void sortRecords();

	void writeOutHeaderFieldsOtherThanFormat(std::ostream & out) const;

	/**
	 * \brief write out the header and the CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO fields
	 * \param out output stream
	 * \param selectRegions if supplied, will only output regions that intersect with these regions
	 */
	void writeOutFixedOnly(std::ostream & out, const std::vector<GenomicRegion> & selectRegions = {}) const;




	/**
	 * \brief write out the header and the CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO fields and then also FORMAT column and the sample columns after that
	 * \param out the output stream to write to
	 * \param selectRegions if supplied, will only output regions that intersect with these regions
	 */
	void writeOutFixedAndSampleMeta(std::ostream & out, const std::vector<GenomicRegion> & selectRegions = {}) const;


	/**
	 * \brief read in the header of a vcf file and initialize
	 * \param fnp the file name path to the input VCF file
	 * \return return a VCFOutput with the header info filed in, no records read in yet
	 */
	static VCFOutput readInHeader(const bfs::path & fnp);

	/**
	 * \brief add in records from vcf file, will skip all lines beginning with "#" and does not check header compatibility, this should be handled elsewere
	 * \param in the input stream of a vcf file
	 */
	void addInRecordsFromFile(std::istream & in);

	/**
 * \brief add in records from vcf file, will skip all lines beginning with "#" and does not check header compatibility, this should be handled elsewere, this function reads in just the first 8 columns, ignores the sample meta data
 * \param in the input stream of a vcf file
 */
	void addInRecordsFixedDataFromFile(std::istream & in);
	[[nodiscard]] VCFRecord processRecordLineForFixedData(const std::string & line) const;
	[[nodiscard]] VCFRecord processRecordLineForFixedDataAndSampleMetaData(const std::string & line) const;

	[[nodiscard]] size_t expectedColumnNumber() const;

	void addInBlnaksForAnyMissingSamples(const std::set<std::string>& samples);


	[[nodiscard]] Json::Value headerToJson() const;



	static VCFOutput comnbineVCFs(const std::vector<bfs::path> & vcfsFnps,
		const std::set<std::string> & sampleNamesSet,
		bool doNotRescueVariantCallsAccrossTargets);

};


} // namespace njhseq
