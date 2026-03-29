// Auto-generated header from JSON Schema
#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <optional>
#include <map>
#include <set>
#include <unordered_set>
#include <utility>
#include <stdexcept>
#include <regex>
#include <limits>
#include <type_traits>
#include <nlohmann/json.hpp>

#ifndef PMO_NA_TOKENS
#define PMO_NA_TOKENS "N/A","NA","Not Applicable",""
#endif
inline constexpr const char* PMO_NA_TOKENS_ARR[] = { PMO_NA_TOKENS };
inline constexpr std::size_t PMO_NA_TOKENS_COUNT = sizeof(PMO_NA_TOKENS_ARR) / sizeof(const char*);

namespace njhseq::pmo {
class BioMethod;
class BioinformaticsMethodInfo;
class BioinformaticsRunInfo;
class MicrohaplotypeForTarget;
class DetectedMicrohaplotypesForTarget;
class DetectedMicrohaplotypesForSample;
class DetectedMicrohaplotypes;
class GenomeInfo;
class GenomicLocation;
class ParasiteDensity;
class PlateInfo;
class LibrarySampleInfo;
class MarkerOfInterest;
class MaskingInfo;
class ReactionInfo;
class PanelInfo;
class PmoGenerationMethod;
class PmoHeader;
class ProjectInfo;
class PrimerInfo;
class TargetInfo;
class ProteinVariant;
class Pseudocigar;
class RepresentativeMicrohaplotype;
class RepresentativeMicrohaplotypesForTarget;
class RepresentativeMicrohaplotypes;
class StageReadCounts;
class ReadCountsByStageForTarget;
class ReadCountsByStageForLibrarySample;
class ReadCountsByStage;
class SequencingInfo;
class TravelInfo;
class SpecimenInfo;
class PortableMicrohaplotypeObject;

const std::unordered_set<std::string>& PMO_NA_STRINGS();

class BioMethod {
public:
  BioMethod() = default;
  std::optional<std::vector<std::string>> additional_argument_;
  std::string program_;
  std::optional<std::string> program_description_;
  std::optional<std::string> program_url_;
  std::string program_version_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static BioMethod from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class BioinformaticsMethodInfo {
public:
  BioinformaticsMethodInfo() = default;
  std::vector<BioMethod> methods_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static BioinformaticsMethodInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class BioinformaticsRunInfo {
public:
  BioinformaticsRunInfo() = default;
  uint32_t bioinformatics_methods_id_{std::numeric_limits<uint32_t>::max()};
  std::string bioinformatics_run_name_;
  std::optional<std::string> run_date_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static BioinformaticsRunInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class MicrohaplotypeForTarget {
public:
  MicrohaplotypeForTarget() = default;
  uint32_t mhap_id_{std::numeric_limits<uint32_t>::max()};
  uint32_t reads_{std::numeric_limits<uint32_t>::max()};
  std::optional<uint32_t> umis_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static MicrohaplotypeForTarget from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class DetectedMicrohaplotypesForTarget {
public:
  DetectedMicrohaplotypesForTarget() = default;
  std::vector<MicrohaplotypeForTarget> mhaps_;
  uint32_t mhaps_target_id_{std::numeric_limits<uint32_t>::max()};
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static DetectedMicrohaplotypesForTarget from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class DetectedMicrohaplotypesForSample {
public:
  DetectedMicrohaplotypesForSample() = default;
  uint32_t library_sample_id_{std::numeric_limits<uint32_t>::max()};
  std::vector<DetectedMicrohaplotypesForTarget> target_results_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static DetectedMicrohaplotypesForSample from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class DetectedMicrohaplotypes {
public:
  DetectedMicrohaplotypes() = default;
  std::optional<uint32_t> bioinformatics_run_id_;
  std::vector<DetectedMicrohaplotypesForSample> library_samples_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static DetectedMicrohaplotypes from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class GenomeInfo {
public:
  GenomeInfo() = default;
  std::optional<std::vector<std::string>> chromosomes_;
  std::string genome_version_;
  std::optional<std::string> gff_url_;
  std::string name_;
  std::vector<uint32_t> taxon_id_;
  std::string url_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static GenomeInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class GenomicLocation {
public:
  GenomicLocation() = default;
  std::optional<std::string> alt_seq_;
  std::string chrom_;
  uint32_t end_{std::numeric_limits<uint32_t>::max()};
  uint32_t genome_id_{std::numeric_limits<uint32_t>::max()};
  std::optional<std::string> ref_seq_;
  uint32_t start_{std::numeric_limits<uint32_t>::max()};
  std::optional<std::string> strand_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static GenomicLocation from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ParasiteDensity {
public:
  ParasiteDensity() = default;
  std::optional<std::string> date_measured_;
  std::optional<std::string> density_method_comments_;
  float parasite_density_{std::numeric_limits<float>::max()};
  std::string parasite_density_method_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ParasiteDensity from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class PlateInfo {
public:
  PlateInfo() = default;
  uint32_t plate_col_{std::numeric_limits<uint32_t>::max()};
  std::string plate_name_;
  std::string plate_row_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static PlateInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class LibrarySampleInfo {
public:
  LibrarySampleInfo() = default;
  std::optional<std::vector<std::string>> alternate_identifiers_;
  std::optional<std::string> experiment_accession_;
  std::optional<std::string> fastqs_loc_;
  std::optional<PlateInfo> library_prep_plate_info_;
  std::string library_sample_name_;
  uint32_t panel_id_{std::numeric_limits<uint32_t>::max()};
  std::optional<std::vector<ParasiteDensity>> qpcr_parasite_density_info_;
  std::optional<std::string> run_accession_;
  std::optional<uint32_t> sequencing_info_id_;
  uint32_t specimen_id_{std::numeric_limits<uint32_t>::max()};
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static LibrarySampleInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class MarkerOfInterest {
public:
  MarkerOfInterest() = default;
  std::optional<std::vector<std::string>> associations_;
  GenomicLocation marker_location_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static MarkerOfInterest from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class MaskingInfo {
public:
  MaskingInfo() = default;
  std::optional<std::string> masking_generation_description_;
  uint32_t replacement_size_{std::numeric_limits<uint32_t>::max()};
  uint32_t seq_segment_size_{std::numeric_limits<uint32_t>::max()};
  uint32_t seq_start_{std::numeric_limits<uint32_t>::max()};
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static MaskingInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ReactionInfo {
public:
  ReactionInfo() = default;
  std::vector<uint32_t> panel_targets_;
  std::string reaction_name_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ReactionInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class PanelInfo {
public:
  PanelInfo() = default;
  std::string panel_name_;
  std::vector<ReactionInfo> reactions_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static PanelInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class PmoGenerationMethod {
public:
  PmoGenerationMethod() = default;
  std::string program_name_;
  std::string program_version_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static PmoGenerationMethod from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class PmoHeader {
public:
  PmoHeader() = default;
  std::optional<std::string> creation_date_;
  std::optional<PmoGenerationMethod> generation_method_;
  std::string pmo_version_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static PmoHeader from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ProjectInfo {
public:
  ProjectInfo() = default;
  std::optional<std::string> BioProject_accession_;
  std::optional<std::string> project_collector_chief_scientist_;
  std::optional<std::vector<std::string>> project_contributors_;
  std::string project_description_;
  std::string project_name_;
  std::optional<std::string> project_type_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ProjectInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class PrimerInfo {
public:
  PrimerInfo() = default;
  std::optional<GenomicLocation> location_;
  std::string seq_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static PrimerInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class TargetInfo {
public:
  TargetInfo() = default;
  PrimerInfo forward_primer_;
  std::optional<std::string> gene_name_;
  std::optional<GenomicLocation> insert_location_;
  std::optional<std::vector<MarkerOfInterest>> markers_of_interest_;
  PrimerInfo reverse_primer_;
  std::optional<std::vector<std::string>> target_attributes_;
  std::string target_name_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static TargetInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ProteinVariant {
public:
  ProteinVariant() = default;
  std::optional<std::string> alternative_gene_name_;
  std::optional<GenomicLocation> codon_genomic_location_;
  std::optional<std::string> gene_name_;
  GenomicLocation protein_location_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ProteinVariant from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class Pseudocigar {
public:
  Pseudocigar() = default;
  std::optional<std::string> pseudocigar_generation_description_;
  std::string pseudocigar_seq_;
  GenomicLocation ref_loc_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static Pseudocigar from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class RepresentativeMicrohaplotype {
public:
  RepresentativeMicrohaplotype() = default;
  std::optional<std::vector<std::string>> alt_annotations_;
  std::optional<std::vector<ProteinVariant>> associated_protein_variants_;
  std::optional<std::vector<GenomicLocation>> associated_seq_variants_;
  std::optional<std::vector<MaskingInfo>> masking_;
  std::optional<std::string> microhaplotype_name_;
  std::optional<Pseudocigar> pseudocigar_;
  std::optional<std::string> quality_;
  std::string seq_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static RepresentativeMicrohaplotype from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class RepresentativeMicrohaplotypesForTarget {
public:
  RepresentativeMicrohaplotypesForTarget() = default;
  std::optional<GenomicLocation> mhap_location_;
  std::vector<RepresentativeMicrohaplotype> microhaplotypes_;
  uint32_t target_id_{std::numeric_limits<uint32_t>::max()};
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static RepresentativeMicrohaplotypesForTarget from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class RepresentativeMicrohaplotypes {
public:
  RepresentativeMicrohaplotypes() = default;
  std::vector<RepresentativeMicrohaplotypesForTarget> targets_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static RepresentativeMicrohaplotypes from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class StageReadCounts {
public:
  StageReadCounts() = default;
  uint32_t reads_{std::numeric_limits<uint32_t>::max()};
  std::string stage_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static StageReadCounts from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ReadCountsByStageForTarget {
public:
  ReadCountsByStageForTarget() = default;
  std::vector<StageReadCounts> stages_;
  uint32_t target_id_{std::numeric_limits<uint32_t>::max()};
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ReadCountsByStageForTarget from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ReadCountsByStageForLibrarySample {
public:
  ReadCountsByStageForLibrarySample() = default;
  uint32_t library_sample_id_{std::numeric_limits<uint32_t>::max()};
  std::optional<std::vector<ReadCountsByStageForTarget>> read_counts_for_targets_;
  uint32_t total_raw_count_{std::numeric_limits<uint32_t>::max()};
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ReadCountsByStageForLibrarySample from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class ReadCountsByStage {
public:
  ReadCountsByStage() = default;
  std::optional<uint32_t> bioinformatics_run_id_;
  std::vector<ReadCountsByStageForLibrarySample> read_counts_by_library_sample_by_stage_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static ReadCountsByStage from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class SequencingInfo {
public:
  SequencingInfo() = default;
  std::optional<std::string> library_kit_;
  std::string library_layout_;
  std::optional<std::string> library_screen_;
  std::string library_selection_;
  std::string library_source_;
  std::string library_strategy_;
  std::optional<std::string> nucl_acid_amp_;
  std::optional<std::string> nucl_acid_amp_date_;
  std::optional<std::string> nucl_acid_ext_;
  std::optional<std::string> nucl_acid_ext_date_;
  std::optional<std::string> pcr_cond_;
  std::optional<std::string> seq_center_;
  std::optional<std::string> seq_date_;
  std::string seq_instrument_model_;
  std::string seq_platform_;
  std::string sequencing_info_name_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static SequencingInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class TravelInfo {
public:
  TravelInfo() = default;
  std::optional<float> bed_net_usage_;
  std::optional<std::string> geo_admin1_;
  std::optional<std::string> geo_admin2_;
  std::optional<std::string> geo_admin3_;
  std::optional<std::string> lat_lon_;
  std::string travel_country_;
  std::string travel_end_date_;
  std::string travel_start_date_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static TravelInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class SpecimenInfo {
public:
  SpecimenInfo() = default;
  std::optional<std::vector<std::string>> alternate_identifiers_;
  std::optional<bool> blood_meal_;
  std::optional<std::string> collection_country_;
  std::optional<std::string> collection_date_;
  std::optional<std::vector<std::string>> drug_usage_;
  std::optional<std::string> env_broad_scale_;
  std::optional<std::string> env_local_scale_;
  std::optional<std::string> env_medium_;
  std::optional<std::string> geo_admin1_;
  std::optional<std::string> geo_admin2_;
  std::optional<std::string> geo_admin3_;
  std::optional<bool> gravid_;
  std::optional<uint32_t> gravidity_;
  std::optional<bool> has_travel_out_six_month_;
  std::optional<float> host_age_;
  std::optional<std::string> host_sex_;
  std::optional<std::string> host_subject_name_;
  std::optional<uint32_t> host_taxon_id_;
  std::optional<std::string> lat_lon_;
  std::optional<std::vector<ParasiteDensity>> parasite_density_info_;
  std::optional<uint32_t> project_id_;
  std::optional<std::string> specimen_accession_;
  std::optional<std::string> specimen_collect_device_;
  std::optional<std::vector<std::string>> specimen_comments_;
  std::string specimen_name_;
  std::optional<std::string> specimen_store_loc_;
  std::optional<std::vector<uint32_t>> specimen_taxon_id_;
  std::optional<std::string> specimen_type_;
  std::optional<PlateInfo> storage_plate_info_;
  std::optional<std::vector<TravelInfo>> travel_out_six_month_;
  std::optional<std::vector<std::string>> treatment_status_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static SpecimenInfo from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

class PortableMicrohaplotypeObject {
public:
  PortableMicrohaplotypeObject() = default;
  std::optional<std::vector<BioinformaticsMethodInfo>> bioinformatics_methods_info_;
  std::optional<std::vector<BioinformaticsRunInfo>> bioinformatics_run_info_;
  std::vector<DetectedMicrohaplotypes> detected_microhaplotypes_;
  std::vector<LibrarySampleInfo> library_sample_info_;
  std::vector<PanelInfo> panel_info_;
  PmoHeader pmo_header_;
  std::optional<std::vector<ProjectInfo>> project_info_;
  std::optional<std::vector<ReadCountsByStage>> read_counts_by_stage_;
  RepresentativeMicrohaplotypes representative_microhaplotypes_;
  std::optional<std::vector<SequencingInfo>> sequencing_info_;
  std::vector<SpecimenInfo> specimen_info_;
  std::vector<TargetInfo> target_info_;
  std::optional<std::vector<GenomeInfo>> targeted_genomes_;
  std::map<std::string, nlohmann::json> extras_;

  [[nodiscard]] static PortableMicrohaplotypeObject from_json(const nlohmann::json& j);
  [[nodiscard]] nlohmann::json to_json() const;
  void validate() const;
};

}
