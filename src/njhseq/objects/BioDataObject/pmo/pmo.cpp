
#include "pmo.hpp"

#include <utility>
#include <regex>
#include <set>


namespace njhseq::pmo {
const std::unordered_set<std::string>& PMO_NA_STRINGS() {
    static const std::unordered_set<std::string> s = [] {
        std::unordered_set<std::string> tmp;
        for (auto i : PMO_NA_TOKENS_ARR) tmp.emplace(i);
        return tmp;
    }();
    return s;
}

BioMethod BioMethod::from_json(const nlohmann::json& j) {
    BioMethod ret;
    std::set<std::string> _known;
    _known.insert("additional_argument");
    if (auto it = j.find("additional_argument"); it != j.end() && !it->is_null()) ret.additional_argument_ = it->get<std::vector<std::string>>();
    _known.insert("program");
    j.at("program").get_to(ret.program_);
    _known.insert("program_description");
    if (auto it = j.find("program_description"); it != j.end() && !it->is_null()) ret.program_description_ = it->get<std::string>();
    _known.insert("program_url");
    if (auto it = j.find("program_url"); it != j.end() && !it->is_null()) ret.program_url_ = it->get<std::string>();
    _known.insert("program_version");
    j.at("program_version").get_to(ret.program_version_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json BioMethod::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (additional_argument_.has_value()) j["additional_argument"] = additional_argument_.value();
    j["program"] = program_;
    if (program_description_.has_value()) j["program_description"] = program_description_.value();
    if (program_url_.has_value()) j["program_url"] = program_url_.value();
    j["program_version"] = program_version_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void BioMethod::validate() const {
    if (!std::regex_match(program_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: BioMethod.program pattern");
    if (!std::regex_match(program_version_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: BioMethod.program_version pattern");
}

BioinformaticsMethodInfo BioinformaticsMethodInfo::from_json(const nlohmann::json& j) {
    BioinformaticsMethodInfo ret;
    std::set<std::string> _known;
    _known.insert("methods");
    { std::vector<BioMethod> _vec; for (const auto& _el : j.at("methods")) _vec.emplace_back(BioMethod::from_json(_el)); ret.methods_ = std::move(_vec); }
    _known.insert("pipeline");
    if (auto it = j.find("pipeline"); it != j.end() && !it->is_null()) ret.pipeline_ = BioMethod::from_json(*it);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json BioinformaticsMethodInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : methods_) _arr.emplace_back(_el.to_json()); j["methods"] = std::move(_arr); }
    if (pipeline_.has_value()) j["pipeline"] = pipeline_.value().to_json();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void BioinformaticsMethodInfo::validate() const {
    for (const auto& _el : methods_) _el.validate();
    if (pipeline_.has_value()) pipeline_.value().validate();
}

BioinformaticsRunInfo BioinformaticsRunInfo::from_json(const nlohmann::json& j) {
    BioinformaticsRunInfo ret;
    std::set<std::string> _known;
    _known.insert("bioinformatics_methods_id");
    if (auto it = j.find("bioinformatics_methods_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.bioinformatics_methods_id_ = std::numeric_limits<uint32_t>::max(); else ret.bioinformatics_methods_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("bioinformatics_methods_id").get_to(ret.bioinformatics_methods_id_);
    _known.insert("bioinformatics_run_name");
    j.at("bioinformatics_run_name").get_to(ret.bioinformatics_run_name_);
    _known.insert("run_date");
    if (auto it = j.find("run_date"); it != j.end() && !it->is_null()) ret.run_date_ = it->get<std::string>();
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json BioinformaticsRunInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (bioinformatics_methods_id_ == std::numeric_limits<uint32_t>::max()) j["bioinformatics_methods_id"] = "NA"; else j["bioinformatics_methods_id"] = bioinformatics_methods_id_;
    j["bioinformatics_run_name"] = bioinformatics_run_name_;
    if (run_date_.has_value()) j["run_date"] = run_date_.value();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void BioinformaticsRunInfo::validate() const {
    if (bioinformatics_methods_id_ < 0) throw std::runtime_error("Validation failed: BioinformaticsRunInfo.bioinformatics_methods_id minimum");
    if (!std::regex_match(bioinformatics_run_name_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: BioinformaticsRunInfo.bioinformatics_run_name pattern");
}

MicrohaplotypeForTarget MicrohaplotypeForTarget::from_json(const nlohmann::json& j) {
    MicrohaplotypeForTarget ret;
    std::set<std::string> _known;
    _known.insert("mhap_id");
    if (auto it = j.find("mhap_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.mhap_id_ = std::numeric_limits<uint32_t>::max(); else ret.mhap_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("mhap_id").get_to(ret.mhap_id_);
    _known.insert("reads");
    if (auto it = j.find("reads"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.reads_ = std::numeric_limits<uint32_t>::max(); else ret.reads_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("reads").get_to(ret.reads_);
    _known.insert("umis");
    if (auto it = j.find("umis"); it != j.end() && !it->is_null()) {
        if (it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.umis_ = std::numeric_limits<uint32_t>::max(); else ret.umis_ = static_cast<uint32_t>(std::stod(_s)); }
        else ret.umis_ = it->get<uint32_t>();
    }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json MicrohaplotypeForTarget::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (mhap_id_ == std::numeric_limits<uint32_t>::max()) j["mhap_id"] = "NA"; else j["mhap_id"] = mhap_id_;
    if (reads_ == std::numeric_limits<uint32_t>::max()) j["reads"] = "NA"; else j["reads"] = reads_;
    if (umis_.has_value()) { if (umis_.value() == std::numeric_limits<uint32_t>::max()) j["umis"] = "NA"; else j["umis"] = umis_.value(); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void MicrohaplotypeForTarget::validate() const {
    if (mhap_id_ < 0) throw std::runtime_error("Validation failed: MicrohaplotypeForTarget.mhap_id minimum");
    if (reads_ < 0) throw std::runtime_error("Validation failed: MicrohaplotypeForTarget.reads minimum");
}

DetectedMicrohaplotypesForTarget DetectedMicrohaplotypesForTarget::from_json(const nlohmann::json& j) {
    DetectedMicrohaplotypesForTarget ret;
    std::set<std::string> _known;
    _known.insert("mhaps");
    { std::vector<MicrohaplotypeForTarget> _vec; for (const auto& _el : j.at("mhaps")) _vec.emplace_back(MicrohaplotypeForTarget::from_json(_el)); ret.mhaps_ = std::move(_vec); }
    _known.insert("mhaps_target_id");
    if (auto it = j.find("mhaps_target_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.mhaps_target_id_ = std::numeric_limits<uint32_t>::max(); else ret.mhaps_target_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("mhaps_target_id").get_to(ret.mhaps_target_id_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json DetectedMicrohaplotypesForTarget::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : mhaps_) _arr.emplace_back(_el.to_json()); j["mhaps"] = std::move(_arr); }
    if (mhaps_target_id_ == std::numeric_limits<uint32_t>::max()) j["mhaps_target_id"] = "NA"; else j["mhaps_target_id"] = mhaps_target_id_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void DetectedMicrohaplotypesForTarget::validate() const {
    if (mhaps_target_id_ < 0) throw std::runtime_error("Validation failed: DetectedMicrohaplotypesForTarget.mhaps_target_id minimum");
    for (const auto& _el : mhaps_) _el.validate();
}

DetectedMicrohaplotypesForSample DetectedMicrohaplotypesForSample::from_json(const nlohmann::json& j) {
    DetectedMicrohaplotypesForSample ret;
    std::set<std::string> _known;
    _known.insert("library_sample_id");
    if (auto it = j.find("library_sample_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.library_sample_id_ = std::numeric_limits<uint32_t>::max(); else ret.library_sample_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("library_sample_id").get_to(ret.library_sample_id_);
    _known.insert("target_results");
    { std::vector<DetectedMicrohaplotypesForTarget> _vec; for (const auto& _el : j.at("target_results")) _vec.emplace_back(DetectedMicrohaplotypesForTarget::from_json(_el)); ret.target_results_ = std::move(_vec); }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json DetectedMicrohaplotypesForSample::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (library_sample_id_ == std::numeric_limits<uint32_t>::max()) j["library_sample_id"] = "NA"; else j["library_sample_id"] = library_sample_id_;
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : target_results_) _arr.emplace_back(_el.to_json()); j["target_results"] = std::move(_arr); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void DetectedMicrohaplotypesForSample::validate() const {
    if (library_sample_id_ < 0) throw std::runtime_error("Validation failed: DetectedMicrohaplotypesForSample.library_sample_id minimum");
    for (const auto& _el : target_results_) _el.validate();
}

DetectedMicrohaplotypes DetectedMicrohaplotypes::from_json(const nlohmann::json& j) {
    DetectedMicrohaplotypes ret;
    std::set<std::string> _known;
    _known.insert("bioinformatics_run_id");
    if (auto it = j.find("bioinformatics_run_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.bioinformatics_run_id_ = std::numeric_limits<uint32_t>::max(); else ret.bioinformatics_run_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("bioinformatics_run_id").get_to(ret.bioinformatics_run_id_);
    _known.insert("library_samples");
    { std::vector<DetectedMicrohaplotypesForSample> _vec; for (const auto& _el : j.at("library_samples")) _vec.emplace_back(DetectedMicrohaplotypesForSample::from_json(_el)); ret.library_samples_ = std::move(_vec); }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json DetectedMicrohaplotypes::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (bioinformatics_run_id_ == std::numeric_limits<uint32_t>::max()) j["bioinformatics_run_id"] = "NA"; else j["bioinformatics_run_id"] = bioinformatics_run_id_;
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : library_samples_) _arr.emplace_back(_el.to_json()); j["library_samples"] = std::move(_arr); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void DetectedMicrohaplotypes::validate() const {
    if (bioinformatics_run_id_ < 0) throw std::runtime_error("Validation failed: DetectedMicrohaplotypes.bioinformatics_run_id minimum");
    for (const auto& _el : library_samples_) _el.validate();
}

GenomeInfo GenomeInfo::from_json(const nlohmann::json& j) {
    GenomeInfo ret;
    std::set<std::string> _known;
    _known.insert("chromosomes");
    if (auto it = j.find("chromosomes"); it != j.end() && !it->is_null()) ret.chromosomes_ = it->get<std::vector<std::string>>();
    _known.insert("genome_version");
    j.at("genome_version").get_to(ret.genome_version_);
    _known.insert("gff_url");
    if (auto it = j.find("gff_url"); it != j.end() && !it->is_null()) ret.gff_url_ = it->get<std::string>();
    _known.insert("name");
    j.at("name").get_to(ret.name_);
    _known.insert("taxon_id");
    { std::vector<uint32_t> _vec; for (const auto& _el : j.at("taxon_id")) { if (_el.is_string()) { auto _s = _el.get<std::string>(); if (PMO_NA_STRINGS().count(_s)) _vec.emplace_back(std::numeric_limits<uint32_t>::max()); else _vec.emplace_back(static_cast<uint32_t>(std::stod(_s))); } else _vec.emplace_back(_el.get<uint32_t>()); } ret.taxon_id_ = std::move(_vec); }
    _known.insert("url");
    j.at("url").get_to(ret.url_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json GenomeInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (chromosomes_.has_value()) j["chromosomes"] = chromosomes_.value();
    j["genome_version"] = genome_version_;
    if (gff_url_.has_value()) j["gff_url"] = gff_url_.value();
    j["name"] = name_;
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _v : taxon_id_) { if (_v == std::numeric_limits<uint32_t>::max()) _arr.emplace_back("NA"); else _arr.emplace_back(_v); } j["taxon_id"] = std::move(_arr); }
    j["url"] = url_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void GenomeInfo::validate() const {
    if (!std::regex_match(genome_version_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: GenomeInfo.genome_version pattern");
    if (!std::regex_match(name_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: GenomeInfo.name pattern");
    if (!std::regex_match(url_, std::regex(R"(^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$)"))) throw std::runtime_error("Validation failed: GenomeInfo.url pattern");
}

GenomicLocation GenomicLocation::from_json(const nlohmann::json& j) {
    GenomicLocation ret;
    std::set<std::string> _known;
    _known.insert("alt_seq");
    if (auto it = j.find("alt_seq"); it != j.end() && !it->is_null()) ret.alt_seq_ = it->get<std::string>();
    _known.insert("chrom");
    j.at("chrom").get_to(ret.chrom_);
    _known.insert("end");
    if (auto it = j.find("end"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.end_ = std::numeric_limits<uint32_t>::max(); else ret.end_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("end").get_to(ret.end_);
    _known.insert("genome_id");
    if (auto it = j.find("genome_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.genome_id_ = std::numeric_limits<uint32_t>::max(); else ret.genome_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("genome_id").get_to(ret.genome_id_);
    _known.insert("ref_seq");
    if (auto it = j.find("ref_seq"); it != j.end() && !it->is_null()) ret.ref_seq_ = it->get<std::string>();
    _known.insert("start");
    if (auto it = j.find("start"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.start_ = std::numeric_limits<uint32_t>::max(); else ret.start_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("start").get_to(ret.start_);
    _known.insert("strand");
    if (auto it = j.find("strand"); it != j.end() && !it->is_null()) ret.strand_ = it->get<std::string>();
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json GenomicLocation::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (alt_seq_.has_value()) j["alt_seq"] = alt_seq_.value();
    j["chrom"] = chrom_;
    if (end_ == std::numeric_limits<uint32_t>::max()) j["end"] = "NA"; else j["end"] = end_;
    if (genome_id_ == std::numeric_limits<uint32_t>::max()) j["genome_id"] = "NA"; else j["genome_id"] = genome_id_;
    if (ref_seq_.has_value()) j["ref_seq"] = ref_seq_.value();
    if (start_ == std::numeric_limits<uint32_t>::max()) j["start"] = "NA"; else j["start"] = start_;
    if (strand_.has_value()) j["strand"] = strand_.value();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void GenomicLocation::validate() const {
    if (!std::regex_match(chrom_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: GenomicLocation.chrom pattern");
    if (end_ < 0) throw std::runtime_error("Validation failed: GenomicLocation.end minimum");
    if (genome_id_ < 0) throw std::runtime_error("Validation failed: GenomicLocation.genome_id minimum");
    if (start_ < 0) throw std::runtime_error("Validation failed: GenomicLocation.start minimum");
}

PlateInfo PlateInfo::from_json(const nlohmann::json& j) {
    PlateInfo ret;
    std::set<std::string> _known;
    _known.insert("plate_col");
    if (auto it = j.find("plate_col"); it != j.end() && !it->is_null()) {
        if (it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.plate_col_ = std::numeric_limits<uint32_t>::max(); else ret.plate_col_ = static_cast<uint32_t>(std::stod(_s)); }
        else ret.plate_col_ = it->get<uint32_t>();
    }
    _known.insert("plate_name");
    if (auto it = j.find("plate_name"); it != j.end() && !it->is_null()) ret.plate_name_ = it->get<std::string>();
    _known.insert("plate_row");
    if (auto it = j.find("plate_row"); it != j.end() && !it->is_null()) ret.plate_row_ = it->get<std::string>();
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json PlateInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (plate_col_.has_value()) { if (plate_col_.value() == std::numeric_limits<uint32_t>::max()) j["plate_col"] = "NA"; else j["plate_col"] = plate_col_.value(); }
    if (plate_name_.has_value()) j["plate_name"] = plate_name_.value();
    if (plate_row_.has_value()) j["plate_row"] = plate_row_.value();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void PlateInfo::validate() const {
}

ParasiteDensity ParasiteDensity::from_json(const nlohmann::json& j) {
    ParasiteDensity ret;
    std::set<std::string> _known;
    _known.insert("date_measured");
    if (auto it = j.find("date_measured"); it != j.end() && !it->is_null()) ret.date_measured_ = it->get<std::string>();
    _known.insert("density_method_comments");
    if (auto it = j.find("density_method_comments"); it != j.end() && !it->is_null()) ret.density_method_comments_ = it->get<std::string>();
    _known.insert("parasite_density");
    if (auto it = j.find("parasite_density"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.parasite_density_ = std::numeric_limits<float>::max(); else ret.parasite_density_ = static_cast<float>(std::stod(_s)); }
    else j.at("parasite_density").get_to(ret.parasite_density_);
    _known.insert("parasite_density_method");
    j.at("parasite_density_method").get_to(ret.parasite_density_method_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ParasiteDensity::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (date_measured_.has_value()) j["date_measured"] = date_measured_.value();
    if (density_method_comments_.has_value()) j["density_method_comments"] = density_method_comments_.value();
    if (parasite_density_ == std::numeric_limits<float>::max()) j["parasite_density"] = "NA"; else j["parasite_density"] = parasite_density_;
    j["parasite_density_method"] = parasite_density_method_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ParasiteDensity::validate() const {
    if (parasite_density_ < 0) throw std::runtime_error("Validation failed: ParasiteDensity.parasite_density minimum");
    if (!std::regex_match(parasite_density_method_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: ParasiteDensity.parasite_density_method pattern");
}

LibrarySampleInfo LibrarySampleInfo::from_json(const nlohmann::json& j) {
    LibrarySampleInfo ret;
    std::set<std::string> _known;
    _known.insert("alternate_identifiers");
    if (auto it = j.find("alternate_identifiers"); it != j.end() && !it->is_null()) ret.alternate_identifiers_ = it->get<std::vector<std::string>>();
    _known.insert("experiment_accession");
    if (auto it = j.find("experiment_accession"); it != j.end() && !it->is_null()) ret.experiment_accession_ = it->get<std::string>();
    _known.insert("fastqs_loc");
    if (auto it = j.find("fastqs_loc"); it != j.end() && !it->is_null()) ret.fastqs_loc_ = it->get<std::string>();
    _known.insert("library_prep_plate_info");
    if (auto it = j.find("library_prep_plate_info"); it != j.end() && !it->is_null()) ret.library_prep_plate_info_ = PlateInfo::from_json(*it);
    _known.insert("library_sample_name");
    j.at("library_sample_name").get_to(ret.library_sample_name_);
    _known.insert("panel_id");
    if (auto it = j.find("panel_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.panel_id_ = std::numeric_limits<uint32_t>::max(); else ret.panel_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("panel_id").get_to(ret.panel_id_);
    _known.insert("qpcr_parasite_density_info");
    if (auto it = j.find("qpcr_parasite_density_info"); it != j.end() && !it->is_null()) { std::vector<ParasiteDensity> _vec; for (const auto& _el : *it) _vec.emplace_back(ParasiteDensity::from_json(_el)); ret.qpcr_parasite_density_info_ = std::move(_vec); }
    _known.insert("run_accession");
    if (auto it = j.find("run_accession"); it != j.end() && !it->is_null()) ret.run_accession_ = it->get<std::string>();
    _known.insert("sequencing_info_id");
    if (auto it = j.find("sequencing_info_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.sequencing_info_id_ = std::numeric_limits<uint32_t>::max(); else ret.sequencing_info_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("sequencing_info_id").get_to(ret.sequencing_info_id_);
    _known.insert("specimen_id");
    if (auto it = j.find("specimen_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.specimen_id_ = std::numeric_limits<uint32_t>::max(); else ret.specimen_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("specimen_id").get_to(ret.specimen_id_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json LibrarySampleInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (alternate_identifiers_.has_value()) j["alternate_identifiers"] = alternate_identifiers_.value();
    if (experiment_accession_.has_value()) j["experiment_accession"] = experiment_accession_.value();
    if (fastqs_loc_.has_value()) j["fastqs_loc"] = fastqs_loc_.value();
    if (library_prep_plate_info_.has_value()) j["library_prep_plate_info"] = library_prep_plate_info_.value().to_json();
    j["library_sample_name"] = library_sample_name_;
    if (panel_id_ == std::numeric_limits<uint32_t>::max()) j["panel_id"] = "NA"; else j["panel_id"] = panel_id_;
    if (qpcr_parasite_density_info_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : qpcr_parasite_density_info_.value()) _arr.emplace_back(_el.to_json()); j["qpcr_parasite_density_info"] = std::move(_arr); }
    if (run_accession_.has_value()) j["run_accession"] = run_accession_.value();
    if (sequencing_info_id_ == std::numeric_limits<uint32_t>::max()) j["sequencing_info_id"] = "NA"; else j["sequencing_info_id"] = sequencing_info_id_;
    if (specimen_id_ == std::numeric_limits<uint32_t>::max()) j["specimen_id"] = "NA"; else j["specimen_id"] = specimen_id_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void LibrarySampleInfo::validate() const {
    if (!std::regex_match(library_sample_name_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: LibrarySampleInfo.library_sample_name pattern");
    if (panel_id_ < 0) throw std::runtime_error("Validation failed: LibrarySampleInfo.panel_id minimum");
    if (sequencing_info_id_ < 0) throw std::runtime_error("Validation failed: LibrarySampleInfo.sequencing_info_id minimum");
    if (specimen_id_ < 0) throw std::runtime_error("Validation failed: LibrarySampleInfo.specimen_id minimum");
    if (library_prep_plate_info_.has_value()) library_prep_plate_info_.value().validate();
    if (qpcr_parasite_density_info_.has_value()) for (const auto& _el : qpcr_parasite_density_info_.value()) _el.validate();
}

MarkerOfInterest MarkerOfInterest::from_json(const nlohmann::json& j) {
    MarkerOfInterest ret;
    std::set<std::string> _known;
    _known.insert("associations");
    if (auto it = j.find("associations"); it != j.end() && !it->is_null()) ret.associations_ = it->get<std::vector<std::string>>();
    _known.insert("marker_location");
    ret.marker_location_ = GenomicLocation::from_json(j.at("marker_location"));
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json MarkerOfInterest::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (associations_.has_value()) j["associations"] = associations_.value();
    j["marker_location"] = marker_location_.to_json();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void MarkerOfInterest::validate() const {
    marker_location_.validate();
}

MaskingInfo MaskingInfo::from_json(const nlohmann::json& j) {
    MaskingInfo ret;
    std::set<std::string> _known;
    _known.insert("masking_generation_description");
    if (auto it = j.find("masking_generation_description"); it != j.end() && !it->is_null()) ret.masking_generation_description_ = it->get<std::string>();
    _known.insert("replacement_size");
    if (auto it = j.find("replacement_size"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.replacement_size_ = std::numeric_limits<uint32_t>::max(); else ret.replacement_size_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("replacement_size").get_to(ret.replacement_size_);
    _known.insert("seq_segment_size");
    if (auto it = j.find("seq_segment_size"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.seq_segment_size_ = std::numeric_limits<uint32_t>::max(); else ret.seq_segment_size_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("seq_segment_size").get_to(ret.seq_segment_size_);
    _known.insert("seq_start");
    if (auto it = j.find("seq_start"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.seq_start_ = std::numeric_limits<uint32_t>::max(); else ret.seq_start_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("seq_start").get_to(ret.seq_start_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json MaskingInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (masking_generation_description_.has_value()) j["masking_generation_description"] = masking_generation_description_.value();
    if (replacement_size_ == std::numeric_limits<uint32_t>::max()) j["replacement_size"] = "NA"; else j["replacement_size"] = replacement_size_;
    if (seq_segment_size_ == std::numeric_limits<uint32_t>::max()) j["seq_segment_size"] = "NA"; else j["seq_segment_size"] = seq_segment_size_;
    if (seq_start_ == std::numeric_limits<uint32_t>::max()) j["seq_start"] = "NA"; else j["seq_start"] = seq_start_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void MaskingInfo::validate() const {
    if (replacement_size_ < 0) throw std::runtime_error("Validation failed: MaskingInfo.replacement_size minimum");
    if (seq_segment_size_ < 0) throw std::runtime_error("Validation failed: MaskingInfo.seq_segment_size minimum");
    if (seq_start_ < 0) throw std::runtime_error("Validation failed: MaskingInfo.seq_start minimum");
}

ReactionInfo ReactionInfo::from_json(const nlohmann::json& j) {
    ReactionInfo ret;
    std::set<std::string> _known;
    _known.insert("panel_targets");
    { std::vector<uint32_t> _vec; for (const auto& _el : j.at("panel_targets")) { if (_el.is_string()) { auto _s = _el.get<std::string>(); if (PMO_NA_STRINGS().count(_s)) _vec.emplace_back(std::numeric_limits<uint32_t>::max()); else _vec.emplace_back(static_cast<uint32_t>(std::stod(_s))); } else _vec.emplace_back(_el.get<uint32_t>()); } ret.panel_targets_ = std::move(_vec); }
    _known.insert("reaction_name");
    j.at("reaction_name").get_to(ret.reaction_name_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ReactionInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _v : panel_targets_) { if (_v == std::numeric_limits<uint32_t>::max()) _arr.emplace_back("NA"); else _arr.emplace_back(_v); } j["panel_targets"] = std::move(_arr); }
    j["reaction_name"] = reaction_name_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ReactionInfo::validate() const {
    if (!std::regex_match(reaction_name_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: ReactionInfo.reaction_name pattern");
}

PanelInfo PanelInfo::from_json(const nlohmann::json& j) {
    PanelInfo ret;
    std::set<std::string> _known;
    _known.insert("panel_name");
    j.at("panel_name").get_to(ret.panel_name_);
    _known.insert("reactions");
    { std::vector<ReactionInfo> _vec; for (const auto& _el : j.at("reactions")) _vec.emplace_back(ReactionInfo::from_json(_el)); ret.reactions_ = std::move(_vec); }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json PanelInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    j["panel_name"] = panel_name_;
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : reactions_) _arr.emplace_back(_el.to_json()); j["reactions"] = std::move(_arr); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void PanelInfo::validate() const {
    if (!std::regex_match(panel_name_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: PanelInfo.panel_name pattern");
    for (const auto& _el : reactions_) _el.validate();
}

PmoGenerationMethod PmoGenerationMethod::from_json(const nlohmann::json& j) {
    PmoGenerationMethod ret;
    std::set<std::string> _known;
    _known.insert("program_name");
    j.at("program_name").get_to(ret.program_name_);
    _known.insert("program_version");
    j.at("program_version").get_to(ret.program_version_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json PmoGenerationMethod::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    j["program_name"] = program_name_;
    j["program_version"] = program_version_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void PmoGenerationMethod::validate() const {
    if (!std::regex_match(program_name_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: PmoGenerationMethod.program_name pattern");
    if (!std::regex_match(program_version_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: PmoGenerationMethod.program_version pattern");
}

PmoHeader PmoHeader::from_json(const nlohmann::json& j) {
    PmoHeader ret;
    std::set<std::string> _known;
    _known.insert("creation_date");
    if (auto it = j.find("creation_date"); it != j.end() && !it->is_null()) ret.creation_date_ = it->get<std::string>();
    _known.insert("generation_method");
    if (auto it = j.find("generation_method"); it != j.end() && !it->is_null()) ret.generation_method_ = PmoGenerationMethod::from_json(*it);
    _known.insert("pmo_version");
    j.at("pmo_version").get_to(ret.pmo_version_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json PmoHeader::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (creation_date_.has_value()) j["creation_date"] = creation_date_.value();
    if (generation_method_.has_value()) j["generation_method"] = generation_method_.value().to_json();
    j["pmo_version"] = pmo_version_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void PmoHeader::validate() const {
    if (!std::regex_match(pmo_version_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: PmoHeader.pmo_version pattern");
    if (generation_method_.has_value()) generation_method_.value().validate();
}

ProteinVariant ProteinVariant::from_json(const nlohmann::json& j) {
    ProteinVariant ret;
    std::set<std::string> _known;
    _known.insert("alternative_gene_name");
    if (auto it = j.find("alternative_gene_name"); it != j.end() && !it->is_null()) ret.alternative_gene_name_ = it->get<std::string>();
    _known.insert("codon_genomic_location");
    if (auto it = j.find("codon_genomic_location"); it != j.end() && !it->is_null()) ret.codon_genomic_location_ = GenomicLocation::from_json(*it);
    _known.insert("gene_name");
    if (auto it = j.find("gene_name"); it != j.end() && !it->is_null()) ret.gene_name_ = it->get<std::string>();
    _known.insert("protein_location");
    ret.protein_location_ = GenomicLocation::from_json(j.at("protein_location"));
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ProteinVariant::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (alternative_gene_name_.has_value()) j["alternative_gene_name"] = alternative_gene_name_.value();
    if (codon_genomic_location_.has_value()) j["codon_genomic_location"] = codon_genomic_location_.value().to_json();
    if (gene_name_.has_value()) j["gene_name"] = gene_name_.value();
    j["protein_location"] = protein_location_.to_json();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ProteinVariant::validate() const {
    if (codon_genomic_location_.has_value()) codon_genomic_location_.value().validate();
    protein_location_.validate();
}

Pseudocigar Pseudocigar::from_json(const nlohmann::json& j) {
    Pseudocigar ret;
    std::set<std::string> _known;
    _known.insert("pseudocigar_generation_description");
    if (auto it = j.find("pseudocigar_generation_description"); it != j.end() && !it->is_null()) ret.pseudocigar_generation_description_ = it->get<std::string>();
    _known.insert("pseudocigar_seq");
    j.at("pseudocigar_seq").get_to(ret.pseudocigar_seq_);
    _known.insert("ref_loc");
    ret.ref_loc_ = GenomicLocation::from_json(j.at("ref_loc"));
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json Pseudocigar::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (pseudocigar_generation_description_.has_value()) j["pseudocigar_generation_description"] = pseudocigar_generation_description_.value();
    j["pseudocigar_seq"] = pseudocigar_seq_;
    j["ref_loc"] = ref_loc_.to_json();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void Pseudocigar::validate() const {
    if (!std::regex_match(pseudocigar_seq_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: Pseudocigar.pseudocigar_seq pattern");
    ref_loc_.validate();
}

RepresentativeMicrohaplotype RepresentativeMicrohaplotype::from_json(const nlohmann::json& j) {
    RepresentativeMicrohaplotype ret;
    std::set<std::string> _known;
    _known.insert("alt_annotations");
    if (auto it = j.find("alt_annotations"); it != j.end() && !it->is_null()) ret.alt_annotations_ = it->get<std::vector<std::string>>();
    _known.insert("associated_protein_variants");
    if (auto it = j.find("associated_protein_variants"); it != j.end() && !it->is_null()) { std::vector<ProteinVariant> _vec; for (const auto& _el : *it) _vec.emplace_back(ProteinVariant::from_json(_el)); ret.associated_protein_variants_ = std::move(_vec); }
    _known.insert("associated_seq_variants");
    if (auto it = j.find("associated_seq_variants"); it != j.end() && !it->is_null()) { std::vector<GenomicLocation> _vec; for (const auto& _el : *it) _vec.emplace_back(GenomicLocation::from_json(_el)); ret.associated_seq_variants_ = std::move(_vec); }
    _known.insert("masking");
    if (auto it = j.find("masking"); it != j.end() && !it->is_null()) { std::vector<MaskingInfo> _vec; for (const auto& _el : *it) _vec.emplace_back(MaskingInfo::from_json(_el)); ret.masking_ = std::move(_vec); }
    _known.insert("microhaplotype_name");
    if (auto it = j.find("microhaplotype_name"); it != j.end() && !it->is_null()) ret.microhaplotype_name_ = it->get<std::string>();
    _known.insert("pseudocigar");
    if (auto it = j.find("pseudocigar"); it != j.end() && !it->is_null()) ret.pseudocigar_ = Pseudocigar::from_json(*it);
    _known.insert("quality");
    if (auto it = j.find("quality"); it != j.end() && !it->is_null()) ret.quality_ = it->get<std::string>();
    _known.insert("seq");
    j.at("seq").get_to(ret.seq_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json RepresentativeMicrohaplotype::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (alt_annotations_.has_value()) j["alt_annotations"] = alt_annotations_.value();
    if (associated_protein_variants_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : associated_protein_variants_.value()) _arr.emplace_back(_el.to_json()); j["associated_protein_variants"] = std::move(_arr); }
    if (associated_seq_variants_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : associated_seq_variants_.value()) _arr.emplace_back(_el.to_json()); j["associated_seq_variants"] = std::move(_arr); }
    if (masking_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : masking_.value()) _arr.emplace_back(_el.to_json()); j["masking"] = std::move(_arr); }
    if (microhaplotype_name_.has_value()) j["microhaplotype_name"] = microhaplotype_name_.value();
    if (pseudocigar_.has_value()) j["pseudocigar"] = pseudocigar_.value().to_json();
    if (quality_.has_value()) j["quality"] = quality_.value();
    j["seq"] = seq_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void RepresentativeMicrohaplotype::validate() const {
    if (!std::regex_match(seq_, std::regex("^[A-z]+$"))) throw std::runtime_error("Validation failed: RepresentativeMicrohaplotype.seq pattern");
    if (associated_protein_variants_.has_value()) for (const auto& _el : associated_protein_variants_.value()) _el.validate();
    if (associated_seq_variants_.has_value()) for (const auto& _el : associated_seq_variants_.value()) _el.validate();
    if (masking_.has_value()) for (const auto& _el : masking_.value()) _el.validate();
    if (pseudocigar_.has_value()) pseudocigar_.value().validate();
}

RepresentativeMicrohaplotypesForTarget RepresentativeMicrohaplotypesForTarget::from_json(const nlohmann::json& j) {
    RepresentativeMicrohaplotypesForTarget ret;
    std::set<std::string> _known;
    _known.insert("mhap_location");
    if (auto it = j.find("mhap_location"); it != j.end() && !it->is_null()) ret.mhap_location_ = GenomicLocation::from_json(*it);
    _known.insert("microhaplotypes");
    { std::vector<RepresentativeMicrohaplotype> _vec; for (const auto& _el : j.at("microhaplotypes")) _vec.emplace_back(RepresentativeMicrohaplotype::from_json(_el)); ret.microhaplotypes_ = std::move(_vec); }
    _known.insert("target_id");
    if (auto it = j.find("target_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.target_id_ = std::numeric_limits<uint32_t>::max(); else ret.target_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("target_id").get_to(ret.target_id_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json RepresentativeMicrohaplotypesForTarget::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (mhap_location_.has_value()) j["mhap_location"] = mhap_location_.value().to_json();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : microhaplotypes_) _arr.emplace_back(_el.to_json()); j["microhaplotypes"] = std::move(_arr); }
    if (target_id_ == std::numeric_limits<uint32_t>::max()) j["target_id"] = "NA"; else j["target_id"] = target_id_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void RepresentativeMicrohaplotypesForTarget::validate() const {
    if (target_id_ < 0) throw std::runtime_error("Validation failed: RepresentativeMicrohaplotypesForTarget.target_id minimum");
    if (mhap_location_.has_value()) mhap_location_.value().validate();
    for (const auto& _el : microhaplotypes_) _el.validate();
}

RepresentativeMicrohaplotypes RepresentativeMicrohaplotypes::from_json(const nlohmann::json& j) {
    RepresentativeMicrohaplotypes ret;
    std::set<std::string> _known;
    _known.insert("targets");
    { std::vector<RepresentativeMicrohaplotypesForTarget> _vec; for (const auto& _el : j.at("targets")) _vec.emplace_back(RepresentativeMicrohaplotypesForTarget::from_json(_el)); ret.targets_ = std::move(_vec); }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json RepresentativeMicrohaplotypes::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : targets_) _arr.emplace_back(_el.to_json()); j["targets"] = std::move(_arr); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void RepresentativeMicrohaplotypes::validate() const {
    for (const auto& _el : targets_) _el.validate();
}

ProjectInfo ProjectInfo::from_json(const nlohmann::json& j) {
    ProjectInfo ret;
    std::set<std::string> _known;
    _known.insert("BioProject_accession");
    if (auto it = j.find("BioProject_accession"); it != j.end() && !it->is_null()) ret.BioProject_accession_ = it->get<std::string>();
    _known.insert("project_collector_chief_scientist");
    if (auto it = j.find("project_collector_chief_scientist"); it != j.end() && !it->is_null()) ret.project_collector_chief_scientist_ = it->get<std::string>();
    _known.insert("project_contributors");
    if (auto it = j.find("project_contributors"); it != j.end() && !it->is_null()) ret.project_contributors_ = it->get<std::vector<std::string>>();
    _known.insert("project_description");
    j.at("project_description").get_to(ret.project_description_);
    _known.insert("project_name");
    j.at("project_name").get_to(ret.project_name_);
    _known.insert("project_type");
    if (auto it = j.find("project_type"); it != j.end() && !it->is_null()) ret.project_type_ = it->get<std::string>();
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ProjectInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (BioProject_accession_.has_value()) j["BioProject_accession"] = BioProject_accession_.value();
    if (project_collector_chief_scientist_.has_value()) j["project_collector_chief_scientist"] = project_collector_chief_scientist_.value();
    if (project_contributors_.has_value()) j["project_contributors"] = project_contributors_.value();
    j["project_description"] = project_description_;
    j["project_name"] = project_name_;
    if (project_type_.has_value()) j["project_type"] = project_type_.value();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ProjectInfo::validate() const {
    if (!std::regex_match(project_name_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: ProjectInfo.project_name pattern");
}

PrimerInfo PrimerInfo::from_json(const nlohmann::json& j) {
    PrimerInfo ret;
    std::set<std::string> _known;
    _known.insert("location");
    if (auto it = j.find("location"); it != j.end() && !it->is_null()) ret.location_ = GenomicLocation::from_json(*it);
    _known.insert("seq");
    j.at("seq").get_to(ret.seq_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json PrimerInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (location_.has_value()) j["location"] = location_.value().to_json();
    j["seq"] = seq_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void PrimerInfo::validate() const {
    if (!std::regex_match(seq_, std::regex("^[A-z]+$"))) throw std::runtime_error("Validation failed: PrimerInfo.seq pattern");
    if (location_.has_value()) location_.value().validate();
}

TargetInfo TargetInfo::from_json(const nlohmann::json& j) {
    TargetInfo ret;
    std::set<std::string> _known;
    _known.insert("forward_primer");
    ret.forward_primer_ = PrimerInfo::from_json(j.at("forward_primer"));
    _known.insert("gene_name");
    if (auto it = j.find("gene_name"); it != j.end() && !it->is_null()) ret.gene_name_ = it->get<std::string>();
    _known.insert("insert_location");
    if (auto it = j.find("insert_location"); it != j.end() && !it->is_null()) ret.insert_location_ = GenomicLocation::from_json(*it);
    _known.insert("markers_of_interest");
    if (auto it = j.find("markers_of_interest"); it != j.end() && !it->is_null()) { std::vector<MarkerOfInterest> _vec; for (const auto& _el : *it) _vec.emplace_back(MarkerOfInterest::from_json(_el)); ret.markers_of_interest_ = std::move(_vec); }
    _known.insert("reverse_primer");
    ret.reverse_primer_ = PrimerInfo::from_json(j.at("reverse_primer"));
    _known.insert("target_attributes");
    if (auto it = j.find("target_attributes"); it != j.end() && !it->is_null()) ret.target_attributes_ = it->get<std::vector<std::string>>();
    _known.insert("target_name");
    j.at("target_name").get_to(ret.target_name_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json TargetInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    j["forward_primer"] = forward_primer_.to_json();
    if (gene_name_.has_value()) j["gene_name"] = gene_name_.value();
    if (insert_location_.has_value()) j["insert_location"] = insert_location_.value().to_json();
    if (markers_of_interest_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : markers_of_interest_.value()) _arr.emplace_back(_el.to_json()); j["markers_of_interest"] = std::move(_arr); }
    j["reverse_primer"] = reverse_primer_.to_json();
    if (target_attributes_.has_value()) j["target_attributes"] = target_attributes_.value();
    j["target_name"] = target_name_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void TargetInfo::validate() const {
    if (!std::regex_match(target_name_, std::regex("^[A-z-._0-9]+$"))) throw std::runtime_error("Validation failed: TargetInfo.target_name pattern");
    forward_primer_.validate();
    if (insert_location_.has_value()) insert_location_.value().validate();
    if (markers_of_interest_.has_value()) for (const auto& _el : markers_of_interest_.value()) _el.validate();
    reverse_primer_.validate();
}

TravelInfo TravelInfo::from_json(const nlohmann::json& j) {
    TravelInfo ret;
    std::set<std::string> _known;
    _known.insert("bed_net_usage");
    if (auto it = j.find("bed_net_usage"); it != j.end() && !it->is_null()) {
        if (it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.bed_net_usage_ = std::numeric_limits<float>::max(); else ret.bed_net_usage_ = static_cast<float>(std::stod(_s)); }
        else ret.bed_net_usage_ = it->get<float>();
    }
    _known.insert("geo_admin1");
    if (auto it = j.find("geo_admin1"); it != j.end() && !it->is_null()) ret.geo_admin1_ = it->get<std::string>();
    _known.insert("geo_admin2");
    if (auto it = j.find("geo_admin2"); it != j.end() && !it->is_null()) ret.geo_admin2_ = it->get<std::string>();
    _known.insert("geo_admin3");
    if (auto it = j.find("geo_admin3"); it != j.end() && !it->is_null()) ret.geo_admin3_ = it->get<std::string>();
    _known.insert("lat_lon");
    if (auto it = j.find("lat_lon"); it != j.end() && !it->is_null()) ret.lat_lon_ = it->get<std::string>();
    _known.insert("travel_country");
    j.at("travel_country").get_to(ret.travel_country_);
    _known.insert("travel_end_date");
    j.at("travel_end_date").get_to(ret.travel_end_date_);
    _known.insert("travel_start_date");
    j.at("travel_start_date").get_to(ret.travel_start_date_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json TravelInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (bed_net_usage_.has_value()) { if (bed_net_usage_.value() == std::numeric_limits<float>::max()) j["bed_net_usage"] = "NA"; else j["bed_net_usage"] = bed_net_usage_.value(); }
    if (geo_admin1_.has_value()) j["geo_admin1"] = geo_admin1_.value();
    if (geo_admin2_.has_value()) j["geo_admin2"] = geo_admin2_.value();
    if (geo_admin3_.has_value()) j["geo_admin3"] = geo_admin3_.value();
    if (lat_lon_.has_value()) j["lat_lon"] = lat_lon_.value();
    j["travel_country"] = travel_country_;
    j["travel_end_date"] = travel_end_date_;
    j["travel_start_date"] = travel_start_date_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void TravelInfo::validate() const {
    if (!std::regex_match(travel_country_, std::regex("^[A-Za-z0-9 ,._:'–-]+$"))) throw std::runtime_error("Validation failed: TravelInfo.travel_country pattern");
    if (!std::regex_match(travel_end_date_, std::regex("\\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?"))) throw std::runtime_error("Validation failed: TravelInfo.travel_end_date pattern");
    if (!std::regex_match(travel_start_date_, std::regex("\\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?"))) throw std::runtime_error("Validation failed: TravelInfo.travel_start_date pattern");
}

SpecimenInfo SpecimenInfo::from_json(const nlohmann::json& j) {
    SpecimenInfo ret;
    std::set<std::string> _known;
    _known.insert("alternate_identifiers");
    if (auto it = j.find("alternate_identifiers"); it != j.end() && !it->is_null()) ret.alternate_identifiers_ = it->get<std::vector<std::string>>();
    _known.insert("blood_meal");
    if (auto it = j.find("blood_meal"); it != j.end() && !it->is_null()) ret.blood_meal_ = it->get<bool>();
    _known.insert("collection_country");
    j.at("collection_country").get_to(ret.collection_country_);
    _known.insert("collection_date");
    j.at("collection_date").get_to(ret.collection_date_);
    _known.insert("drug_usage");
    if (auto it = j.find("drug_usage"); it != j.end() && !it->is_null()) ret.drug_usage_ = it->get<std::vector<std::string>>();
    _known.insert("env_broad_scale");
    if (auto it = j.find("env_broad_scale"); it != j.end() && !it->is_null()) ret.env_broad_scale_ = it->get<std::string>();
    _known.insert("env_local_scale");
    if (auto it = j.find("env_local_scale"); it != j.end() && !it->is_null()) ret.env_local_scale_ = it->get<std::string>();
    _known.insert("env_medium");
    if (auto it = j.find("env_medium"); it != j.end() && !it->is_null()) ret.env_medium_ = it->get<std::string>();
    _known.insert("geo_admin1");
    if (auto it = j.find("geo_admin1"); it != j.end() && !it->is_null()) ret.geo_admin1_ = it->get<std::string>();
    _known.insert("geo_admin2");
    if (auto it = j.find("geo_admin2"); it != j.end() && !it->is_null()) ret.geo_admin2_ = it->get<std::string>();
    _known.insert("geo_admin3");
    if (auto it = j.find("geo_admin3"); it != j.end() && !it->is_null()) ret.geo_admin3_ = it->get<std::string>();
    _known.insert("gravid");
    if (auto it = j.find("gravid"); it != j.end() && !it->is_null()) ret.gravid_ = it->get<bool>();
    _known.insert("gravidity");
    if (auto it = j.find("gravidity"); it != j.end() && !it->is_null()) {
        if (it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.gravidity_ = std::numeric_limits<uint32_t>::max(); else ret.gravidity_ = static_cast<uint32_t>(std::stod(_s)); }
        else ret.gravidity_ = it->get<uint32_t>();
    }
    _known.insert("has_travel_out_six_month");
    if (auto it = j.find("has_travel_out_six_month"); it != j.end() && !it->is_null()) ret.has_travel_out_six_month_ = it->get<bool>();
    _known.insert("host_age");
    if (auto it = j.find("host_age"); it != j.end() && !it->is_null()) {
        if (it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.host_age_ = std::numeric_limits<float>::max(); else ret.host_age_ = static_cast<float>(std::stod(_s)); }
        else ret.host_age_ = it->get<float>();
    }
    _known.insert("host_sex");
    if (auto it = j.find("host_sex"); it != j.end() && !it->is_null()) ret.host_sex_ = it->get<std::string>();
    _known.insert("host_subject_id");
    if (auto it = j.find("host_subject_id"); it != j.end() && !it->is_null()) {
        if (it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.host_subject_id_ = std::numeric_limits<uint32_t>::max(); else ret.host_subject_id_ = static_cast<uint32_t>(std::stod(_s)); }
        else ret.host_subject_id_ = it->get<uint32_t>();
    }
    _known.insert("host_taxon_id");
    if (auto it = j.find("host_taxon_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.host_taxon_id_ = std::numeric_limits<uint32_t>::max(); else ret.host_taxon_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("host_taxon_id").get_to(ret.host_taxon_id_);
    _known.insert("lat_lon");
    if (auto it = j.find("lat_lon"); it != j.end() && !it->is_null()) ret.lat_lon_ = it->get<std::string>();
    _known.insert("parasite_density_info");
    if (auto it = j.find("parasite_density_info"); it != j.end() && !it->is_null()) { std::vector<ParasiteDensity> _vec; for (const auto& _el : *it) _vec.emplace_back(ParasiteDensity::from_json(_el)); ret.parasite_density_info_ = std::move(_vec); }
    _known.insert("project_id");
    if (auto it = j.find("project_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.project_id_ = std::numeric_limits<uint32_t>::max(); else ret.project_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("project_id").get_to(ret.project_id_);
    _known.insert("specimen_accession");
    if (auto it = j.find("specimen_accession"); it != j.end() && !it->is_null()) ret.specimen_accession_ = it->get<std::string>();
    _known.insert("specimen_collect_device");
    if (auto it = j.find("specimen_collect_device"); it != j.end() && !it->is_null()) ret.specimen_collect_device_ = it->get<std::string>();
    _known.insert("specimen_comments");
    if (auto it = j.find("specimen_comments"); it != j.end() && !it->is_null()) ret.specimen_comments_ = it->get<std::vector<std::string>>();
    _known.insert("specimen_name");
    j.at("specimen_name").get_to(ret.specimen_name_);
    _known.insert("specimen_store_loc");
    if (auto it = j.find("specimen_store_loc"); it != j.end() && !it->is_null()) ret.specimen_store_loc_ = it->get<std::string>();
    _known.insert("specimen_taxon_id");
    { std::vector<uint32_t> _vec; for (const auto& _el : j.at("specimen_taxon_id")) { if (_el.is_string()) { auto _s = _el.get<std::string>(); if (PMO_NA_STRINGS().count(_s)) _vec.emplace_back(std::numeric_limits<uint32_t>::max()); else _vec.emplace_back(static_cast<uint32_t>(std::stod(_s))); } else _vec.emplace_back(_el.get<uint32_t>()); } ret.specimen_taxon_id_ = std::move(_vec); }
    _known.insert("specimen_type");
    if (auto it = j.find("specimen_type"); it != j.end() && !it->is_null()) ret.specimen_type_ = it->get<std::string>();
    _known.insert("storage_plate_info");
    if (auto it = j.find("storage_plate_info"); it != j.end() && !it->is_null()) ret.storage_plate_info_ = PlateInfo::from_json(*it);
    _known.insert("travel_out_six_month");
    if (auto it = j.find("travel_out_six_month"); it != j.end() && !it->is_null()) { std::vector<TravelInfo> _vec; for (const auto& _el : *it) _vec.emplace_back(TravelInfo::from_json(_el)); ret.travel_out_six_month_ = std::move(_vec); }
    _known.insert("treatment_status");
    if (auto it = j.find("treatment_status"); it != j.end() && !it->is_null()) ret.treatment_status_ = it->get<std::vector<std::string>>();
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json SpecimenInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (alternate_identifiers_.has_value()) j["alternate_identifiers"] = alternate_identifiers_.value();
    if (blood_meal_.has_value()) j["blood_meal"] = blood_meal_.value();
    j["collection_country"] = collection_country_;
    j["collection_date"] = collection_date_;
    if (drug_usage_.has_value()) j["drug_usage"] = drug_usage_.value();
    if (env_broad_scale_.has_value()) j["env_broad_scale"] = env_broad_scale_.value();
    if (env_local_scale_.has_value()) j["env_local_scale"] = env_local_scale_.value();
    if (env_medium_.has_value()) j["env_medium"] = env_medium_.value();
    if (geo_admin1_.has_value()) j["geo_admin1"] = geo_admin1_.value();
    if (geo_admin2_.has_value()) j["geo_admin2"] = geo_admin2_.value();
    if (geo_admin3_.has_value()) j["geo_admin3"] = geo_admin3_.value();
    if (gravid_.has_value()) j["gravid"] = gravid_.value();
    if (gravidity_.has_value()) { if (gravidity_.value() == std::numeric_limits<uint32_t>::max()) j["gravidity"] = "NA"; else j["gravidity"] = gravidity_.value(); }
    if (has_travel_out_six_month_.has_value()) j["has_travel_out_six_month"] = has_travel_out_six_month_.value();
    if (host_age_.has_value()) { if (host_age_.value() == std::numeric_limits<float>::max()) j["host_age"] = "NA"; else j["host_age"] = host_age_.value(); }
    if (host_sex_.has_value()) j["host_sex"] = host_sex_.value();
    if (host_subject_id_.has_value()) { if (host_subject_id_.value() == std::numeric_limits<uint32_t>::max()) j["host_subject_id"] = "NA"; else j["host_subject_id"] = host_subject_id_.value(); }
    if (host_taxon_id_ == std::numeric_limits<uint32_t>::max()) j["host_taxon_id"] = "NA"; else j["host_taxon_id"] = host_taxon_id_;
    if (lat_lon_.has_value()) j["lat_lon"] = lat_lon_.value();
    if (parasite_density_info_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : parasite_density_info_.value()) _arr.emplace_back(_el.to_json()); j["parasite_density_info"] = std::move(_arr); }
    if (project_id_ == std::numeric_limits<uint32_t>::max()) j["project_id"] = "NA"; else j["project_id"] = project_id_;
    if (specimen_accession_.has_value()) j["specimen_accession"] = specimen_accession_.value();
    if (specimen_collect_device_.has_value()) j["specimen_collect_device"] = specimen_collect_device_.value();
    if (specimen_comments_.has_value()) j["specimen_comments"] = specimen_comments_.value();
    j["specimen_name"] = specimen_name_;
    if (specimen_store_loc_.has_value()) j["specimen_store_loc"] = specimen_store_loc_.value();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _v : specimen_taxon_id_) { if (_v == std::numeric_limits<uint32_t>::max()) _arr.emplace_back("NA"); else _arr.emplace_back(_v); } j["specimen_taxon_id"] = std::move(_arr); }
    if (specimen_type_.has_value()) j["specimen_type"] = specimen_type_.value();
    if (storage_plate_info_.has_value()) j["storage_plate_info"] = storage_plate_info_.value().to_json();
    if (travel_out_six_month_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : travel_out_six_month_.value()) _arr.emplace_back(_el.to_json()); j["travel_out_six_month"] = std::move(_arr); }
    if (treatment_status_.has_value()) j["treatment_status"] = treatment_status_.value();
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void SpecimenInfo::validate() const {
    if (!std::regex_match(collection_country_, std::regex("^[A-Za-z0-9 ,._:'–-]+$"))) throw std::runtime_error("Validation failed: SpecimenInfo.collection_country pattern");
    if (!std::regex_match(collection_date_, std::regex("(?:\\d{4}(?:-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?)?|NA)"))) throw std::runtime_error("Validation failed: SpecimenInfo.collection_date pattern");
    if (host_taxon_id_ < 0) throw std::runtime_error("Validation failed: SpecimenInfo.host_taxon_id minimum");
    if (project_id_ < 0) throw std::runtime_error("Validation failed: SpecimenInfo.project_id minimum");
    if (!std::regex_match(specimen_name_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SpecimenInfo.specimen_name pattern");
    if (parasite_density_info_.has_value()) for (const auto& _el : parasite_density_info_.value()) _el.validate();
    if (storage_plate_info_.has_value()) storage_plate_info_.value().validate();
    if (travel_out_six_month_.has_value()) for (const auto& _el : travel_out_six_month_.value()) _el.validate();
}

SequencingInfo SequencingInfo::from_json(const nlohmann::json& j) {
    SequencingInfo ret;
    std::set<std::string> _known;
    _known.insert("library_kit");
    if (auto it = j.find("library_kit"); it != j.end() && !it->is_null()) ret.library_kit_ = it->get<std::string>();
    _known.insert("library_layout");
    j.at("library_layout").get_to(ret.library_layout_);
    _known.insert("library_screen");
    if (auto it = j.find("library_screen"); it != j.end() && !it->is_null()) ret.library_screen_ = it->get<std::string>();
    _known.insert("library_selection");
    j.at("library_selection").get_to(ret.library_selection_);
    _known.insert("library_source");
    j.at("library_source").get_to(ret.library_source_);
    _known.insert("library_strategy");
    j.at("library_strategy").get_to(ret.library_strategy_);
    _known.insert("nucl_acid_amp");
    if (auto it = j.find("nucl_acid_amp"); it != j.end() && !it->is_null()) ret.nucl_acid_amp_ = it->get<std::string>();
    _known.insert("nucl_acid_amp_date");
    if (auto it = j.find("nucl_acid_amp_date"); it != j.end() && !it->is_null()) ret.nucl_acid_amp_date_ = it->get<std::string>();
    _known.insert("nucl_acid_ext");
    if (auto it = j.find("nucl_acid_ext"); it != j.end() && !it->is_null()) ret.nucl_acid_ext_ = it->get<std::string>();
    _known.insert("nucl_acid_ext_date");
    if (auto it = j.find("nucl_acid_ext_date"); it != j.end() && !it->is_null()) ret.nucl_acid_ext_date_ = it->get<std::string>();
    _known.insert("pcr_cond");
    if (auto it = j.find("pcr_cond"); it != j.end() && !it->is_null()) ret.pcr_cond_ = it->get<std::string>();
    _known.insert("seq_center");
    if (auto it = j.find("seq_center"); it != j.end() && !it->is_null()) ret.seq_center_ = it->get<std::string>();
    _known.insert("seq_date");
    if (auto it = j.find("seq_date"); it != j.end() && !it->is_null()) ret.seq_date_ = it->get<std::string>();
    _known.insert("seq_instrument_model");
    j.at("seq_instrument_model").get_to(ret.seq_instrument_model_);
    _known.insert("seq_platform");
    j.at("seq_platform").get_to(ret.seq_platform_);
    _known.insert("sequencing_info_name");
    j.at("sequencing_info_name").get_to(ret.sequencing_info_name_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json SequencingInfo::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (library_kit_.has_value()) j["library_kit"] = library_kit_.value();
    j["library_layout"] = library_layout_;
    if (library_screen_.has_value()) j["library_screen"] = library_screen_.value();
    j["library_selection"] = library_selection_;
    j["library_source"] = library_source_;
    j["library_strategy"] = library_strategy_;
    if (nucl_acid_amp_.has_value()) j["nucl_acid_amp"] = nucl_acid_amp_.value();
    if (nucl_acid_amp_date_.has_value()) j["nucl_acid_amp_date"] = nucl_acid_amp_date_.value();
    if (nucl_acid_ext_.has_value()) j["nucl_acid_ext"] = nucl_acid_ext_.value();
    if (nucl_acid_ext_date_.has_value()) j["nucl_acid_ext_date"] = nucl_acid_ext_date_.value();
    if (pcr_cond_.has_value()) j["pcr_cond"] = pcr_cond_.value();
    if (seq_center_.has_value()) j["seq_center"] = seq_center_.value();
    if (seq_date_.has_value()) j["seq_date"] = seq_date_.value();
    j["seq_instrument_model"] = seq_instrument_model_;
    j["seq_platform"] = seq_platform_;
    j["sequencing_info_name"] = sequencing_info_name_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void SequencingInfo::validate() const {
    if (!std::regex_match(library_layout_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.library_layout pattern");
    if (!std::regex_match(library_selection_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.library_selection pattern");
    if (!std::regex_match(library_source_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.library_source pattern");
    if (!std::regex_match(library_strategy_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.library_strategy pattern");
    if (!std::regex_match(seq_instrument_model_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.seq_instrument_model pattern");
    if (!std::regex_match(seq_platform_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.seq_platform pattern");
    if (!std::regex_match(sequencing_info_name_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: SequencingInfo.sequencing_info_name pattern");
}

StageReadCounts StageReadCounts::from_json(const nlohmann::json& j) {
    StageReadCounts ret;
    std::set<std::string> _known;
    _known.insert("read_count");
    if (auto it = j.find("read_count"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.read_count_ = std::numeric_limits<uint32_t>::max(); else ret.read_count_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("read_count").get_to(ret.read_count_);
    _known.insert("stage");
    j.at("stage").get_to(ret.stage_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json StageReadCounts::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (read_count_ == std::numeric_limits<uint32_t>::max()) j["read_count"] = "NA"; else j["read_count"] = read_count_;
    j["stage"] = stage_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void StageReadCounts::validate() const {
    if (read_count_ < 0) throw std::runtime_error("Validation failed: StageReadCounts.read_count minimum");
    if (!std::regex_match(stage_, std::regex("^[A-z-._0-9 ]+$"))) throw std::runtime_error("Validation failed: StageReadCounts.stage pattern");
}

ReadCountsByStageForTarget ReadCountsByStageForTarget::from_json(const nlohmann::json& j) {
    ReadCountsByStageForTarget ret;
    std::set<std::string> _known;
    _known.insert("stages");
    { std::vector<StageReadCounts> _vec; for (const auto& _el : j.at("stages")) _vec.emplace_back(StageReadCounts::from_json(_el)); ret.stages_ = std::move(_vec); }
    _known.insert("target_id");
    if (auto it = j.find("target_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.target_id_ = std::numeric_limits<uint32_t>::max(); else ret.target_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("target_id").get_to(ret.target_id_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ReadCountsByStageForTarget::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : stages_) _arr.emplace_back(_el.to_json()); j["stages"] = std::move(_arr); }
    if (target_id_ == std::numeric_limits<uint32_t>::max()) j["target_id"] = "NA"; else j["target_id"] = target_id_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ReadCountsByStageForTarget::validate() const {
    if (target_id_ < 0) throw std::runtime_error("Validation failed: ReadCountsByStageForTarget.target_id minimum");
    for (const auto& _el : stages_) _el.validate();
}

ReadCountsByStageForLibrarySample ReadCountsByStageForLibrarySample::from_json(const nlohmann::json& j) {
    ReadCountsByStageForLibrarySample ret;
    std::set<std::string> _known;
    _known.insert("library_sample_id");
    if (auto it = j.find("library_sample_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.library_sample_id_ = std::numeric_limits<uint32_t>::max(); else ret.library_sample_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("library_sample_id").get_to(ret.library_sample_id_);
    _known.insert("read_counts_for_targets");
    if (auto it = j.find("read_counts_for_targets"); it != j.end() && !it->is_null()) { std::vector<ReadCountsByStageForTarget> _vec; for (const auto& _el : *it) _vec.emplace_back(ReadCountsByStageForTarget::from_json(_el)); ret.read_counts_for_targets_ = std::move(_vec); }
    _known.insert("total_raw_count");
    if (auto it = j.find("total_raw_count"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.total_raw_count_ = std::numeric_limits<uint32_t>::max(); else ret.total_raw_count_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("total_raw_count").get_to(ret.total_raw_count_);
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ReadCountsByStageForLibrarySample::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (library_sample_id_ == std::numeric_limits<uint32_t>::max()) j["library_sample_id"] = "NA"; else j["library_sample_id"] = library_sample_id_;
    if (read_counts_for_targets_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : read_counts_for_targets_.value()) _arr.emplace_back(_el.to_json()); j["read_counts_for_targets"] = std::move(_arr); }
    if (total_raw_count_ == std::numeric_limits<uint32_t>::max()) j["total_raw_count"] = "NA"; else j["total_raw_count"] = total_raw_count_;
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ReadCountsByStageForLibrarySample::validate() const {
    if (library_sample_id_ < 0) throw std::runtime_error("Validation failed: ReadCountsByStageForLibrarySample.library_sample_id minimum");
    if (total_raw_count_ < 0) throw std::runtime_error("Validation failed: ReadCountsByStageForLibrarySample.total_raw_count minimum");
    if (read_counts_for_targets_.has_value()) for (const auto& _el : read_counts_for_targets_.value()) _el.validate();
}

ReadCountsByStage ReadCountsByStage::from_json(const nlohmann::json& j) {
    ReadCountsByStage ret;
    std::set<std::string> _known;
    _known.insert("bioinformatics_run_id");
    if (auto it = j.find("bioinformatics_run_id"); it != j.end() && it->is_string()) { auto _s = it->get<std::string>(); if (PMO_NA_STRINGS().count(_s)) ret.bioinformatics_run_id_ = std::numeric_limits<uint32_t>::max(); else ret.bioinformatics_run_id_ = static_cast<uint32_t>(std::stod(_s)); }
    else j.at("bioinformatics_run_id").get_to(ret.bioinformatics_run_id_);
    _known.insert("read_counts_by_library_sample_by_stage");
    { std::vector<ReadCountsByStageForLibrarySample> _vec; for (const auto& _el : j.at("read_counts_by_library_sample_by_stage")) _vec.emplace_back(ReadCountsByStageForLibrarySample::from_json(_el)); ret.read_counts_by_library_sample_by_stage_ = std::move(_vec); }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json ReadCountsByStage::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    if (bioinformatics_run_id_ == std::numeric_limits<uint32_t>::max()) j["bioinformatics_run_id"] = "NA"; else j["bioinformatics_run_id"] = bioinformatics_run_id_;
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : read_counts_by_library_sample_by_stage_) _arr.emplace_back(_el.to_json()); j["read_counts_by_library_sample_by_stage"] = std::move(_arr); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void ReadCountsByStage::validate() const {
    if (bioinformatics_run_id_ < 0) throw std::runtime_error("Validation failed: ReadCountsByStage.bioinformatics_run_id minimum");
    for (const auto& _el : read_counts_by_library_sample_by_stage_) _el.validate();
}

PortableMicrohaplotypeObject PortableMicrohaplotypeObject::from_json(const nlohmann::json& j) {
    PortableMicrohaplotypeObject ret;
    std::set<std::string> _known;
    _known.insert("bioinformatics_methods_info");
    { std::vector<BioinformaticsMethodInfo> _vec; for (const auto& _el : j.at("bioinformatics_methods_info")) _vec.emplace_back(BioinformaticsMethodInfo::from_json(_el)); ret.bioinformatics_methods_info_ = std::move(_vec); }
    _known.insert("bioinformatics_run_info");
    { std::vector<BioinformaticsRunInfo> _vec; for (const auto& _el : j.at("bioinformatics_run_info")) _vec.emplace_back(BioinformaticsRunInfo::from_json(_el)); ret.bioinformatics_run_info_ = std::move(_vec); }
    _known.insert("detected_microhaplotypes");
    { std::vector<DetectedMicrohaplotypes> _vec; for (const auto& _el : j.at("detected_microhaplotypes")) _vec.emplace_back(DetectedMicrohaplotypes::from_json(_el)); ret.detected_microhaplotypes_ = std::move(_vec); }
    _known.insert("library_sample_info");
    { std::vector<LibrarySampleInfo> _vec; for (const auto& _el : j.at("library_sample_info")) _vec.emplace_back(LibrarySampleInfo::from_json(_el)); ret.library_sample_info_ = std::move(_vec); }
    _known.insert("panel_info");
    { std::vector<PanelInfo> _vec; for (const auto& _el : j.at("panel_info")) _vec.emplace_back(PanelInfo::from_json(_el)); ret.panel_info_ = std::move(_vec); }
    _known.insert("pmo_header");
    ret.pmo_header_ = PmoHeader::from_json(j.at("pmo_header"));
    _known.insert("project_info");
    { std::vector<ProjectInfo> _vec; for (const auto& _el : j.at("project_info")) _vec.emplace_back(ProjectInfo::from_json(_el)); ret.project_info_ = std::move(_vec); }
    _known.insert("read_counts_by_stage");
    if (auto it = j.find("read_counts_by_stage"); it != j.end() && !it->is_null()) { std::vector<ReadCountsByStage> _vec; for (const auto& _el : *it) _vec.emplace_back(ReadCountsByStage::from_json(_el)); ret.read_counts_by_stage_ = std::move(_vec); }
    _known.insert("representative_microhaplotypes");
    ret.representative_microhaplotypes_ = RepresentativeMicrohaplotypes::from_json(j.at("representative_microhaplotypes"));
    _known.insert("sequencing_info");
    { std::vector<SequencingInfo> _vec; for (const auto& _el : j.at("sequencing_info")) _vec.emplace_back(SequencingInfo::from_json(_el)); ret.sequencing_info_ = std::move(_vec); }
    _known.insert("specimen_info");
    { std::vector<SpecimenInfo> _vec; for (const auto& _el : j.at("specimen_info")) _vec.emplace_back(SpecimenInfo::from_json(_el)); ret.specimen_info_ = std::move(_vec); }
    _known.insert("target_info");
    { std::vector<TargetInfo> _vec; for (const auto& _el : j.at("target_info")) _vec.emplace_back(TargetInfo::from_json(_el)); ret.target_info_ = std::move(_vec); }
    _known.insert("targeted_genomes");
    { std::vector<GenomeInfo> _vec; for (const auto& _el : j.at("targeted_genomes")) _vec.emplace_back(GenomeInfo::from_json(_el)); ret.targeted_genomes_ = std::move(_vec); }
    for (auto it = j.begin(); it != j.end(); ++it) { if (_known.find(it.key()) == _known.end()) ret.extras_[it.key()] = it.value(); }
    return ret;
}

nlohmann::json PortableMicrohaplotypeObject::to_json() const {
    nlohmann::json j = nlohmann::json::object();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : bioinformatics_methods_info_) _arr.emplace_back(_el.to_json()); j["bioinformatics_methods_info"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : bioinformatics_run_info_) _arr.emplace_back(_el.to_json()); j["bioinformatics_run_info"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : detected_microhaplotypes_) _arr.emplace_back(_el.to_json()); j["detected_microhaplotypes"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : library_sample_info_) _arr.emplace_back(_el.to_json()); j["library_sample_info"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : panel_info_) _arr.emplace_back(_el.to_json()); j["panel_info"] = std::move(_arr); }
    j["pmo_header"] = pmo_header_.to_json();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : project_info_) _arr.emplace_back(_el.to_json()); j["project_info"] = std::move(_arr); }
    if (read_counts_by_stage_.has_value()) { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : read_counts_by_stage_.value()) _arr.emplace_back(_el.to_json()); j["read_counts_by_stage"] = std::move(_arr); }
    j["representative_microhaplotypes"] = representative_microhaplotypes_.to_json();
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : sequencing_info_) _arr.emplace_back(_el.to_json()); j["sequencing_info"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : specimen_info_) _arr.emplace_back(_el.to_json()); j["specimen_info"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : target_info_) _arr.emplace_back(_el.to_json()); j["target_info"] = std::move(_arr); }
    { nlohmann::json _arr = nlohmann::json::array(); for (const auto& _el : targeted_genomes_) _arr.emplace_back(_el.to_json()); j["targeted_genomes"] = std::move(_arr); }
    for (const auto& kv : extras_) j[kv.first] = kv.second;
    return j;
}

void PortableMicrohaplotypeObject::validate() const {
    for (const auto& _el : bioinformatics_methods_info_) _el.validate();
    for (const auto& _el : bioinformatics_run_info_) _el.validate();
    for (const auto& _el : detected_microhaplotypes_) _el.validate();
    for (const auto& _el : library_sample_info_) _el.validate();
    for (const auto& _el : panel_info_) _el.validate();
    pmo_header_.validate();
    for (const auto& _el : project_info_) _el.validate();
    if (read_counts_by_stage_.has_value()) for (const auto& _el : read_counts_by_stage_.value()) _el.validate();
    representative_microhaplotypes_.validate();
    for (const auto& _el : sequencing_info_) _el.validate();
    for (const auto& _el : specimen_info_) _el.validate();
    for (const auto& _el : target_info_) _el.validate();
    for (const auto& _el : targeted_genomes_) _el.validate();
}

} // namespace njhseq::pmo
