//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "RInsideUtils.hpp"

//////tools for dealing with multipleSampleCollapse
namespace bibseq {
/*
  void installLib (const std::string & lib,
  RInside & rSession){
  std::string txt = "suppressMessages(library(" + lib + "))";
  rSession.parseEvalQ(txt);              // load library, no return value
  }
  std::unordered_map<std::string, std::vector<double>> getModelQualError(RInside
  & R, const std::unordered_map<std::string, std::vector<double>>& data,
  uint32_t qualStart,
  uint32_t qualStop,
  double qualStep,
  const std::string & model){
  return getModelQualError(R, data.at("qual"), data.at("error"),
  data.at("weight"),
  qualStart, qualStop, qualStep, model);
  }
  std::unordered_map<std::string, std::vector<double>> getModelQualError(RInside
  & R, const std::vector<double> & qualVec,
  const std::vector<double> & errorVec,
  const std::vector<double> & weightVec,
  uint32_t qualStart,
  uint32_t qualStop,
  double qualStep,
  const std::string & model){
  Rcpp::DataFrame inData = Rcpp::DataFrame::create(Rcpp::Named("qual") =
  qualVec,
  Rcpp::Named("error") = errorVec,
  Rcpp::Named("weight") = weightVec);
  R["inData"] = inData;
  R.parseEvalQ("myModel <- glm(error ~ qual, family = binomial(link = \"" +
  model + "\"), data = inData, weights= weight);");
  std::string seqStrCmd = to_string(qualStart) + "," + to_string(qualStop) +
  ",by=" + to_string(qualStep);
  R.parseEvalQ("qTable = data.frame(qual = seq("+ seqStrCmd + "));");
  R.parseEvalQ("qTable[, c(\"p\", \"se\")] <- predict(myModel, qTable, type =
  \"response\", se.fit = TRUE)[-3];");
  //R.parseEvalQ("print(qTable);");
  R.parseEvalQ("qTable_mat = as.matrix(qTable);");
  //R.parseEvalQ("cat(qTable_mat, sep =" ");");
  std::string txt = "qTable_mat;";
  Rcpp::NumericMatrix outData = R.parseEval(txt);
  std::unordered_map<std::string, std::vector<double>> out;

  std::vector<double> qualOut;
  for(auto i : iter::range(outData.nrow())){
  qualOut.emplace_back(outData(i,0));
  }
  out["qual"] = qualOut;

  std::vector<double> predOut;
  for(auto i : iter::range(outData.nrow())){
  predOut.emplace_back(outData(i,1));
  }
  out["pred"] = predOut;
  std::vector<double> seOut;
  for(auto i : iter::range(outData.nrow())){
  seOut.emplace_back(outData(i,2));
  }
  out["se"] = seOut;
  return out;
  }*/
void ownRInside::installLib(const std::string& lib) {

  std::string txt = "suppressMessages(library(" + lib + "))";
  auto& r = get();

  r.parseEvalQ(txt);  // load library, no return value
}

std::unordered_map<std::string, std::vector<double>> getModelQualError(
    const std::unordered_map<std::string, std::vector<double>>& data,
    uint32_t qualStart, uint32_t qualStop, double qualStep,
    const std::string& model) {
  return getModelQualError(data.at("qual"), data.at("error"), data.at("weight"),
                           qualStart, qualStop, qualStep, model);
}

std::unordered_map<std::string, std::vector<double>> getModelQualError(
    const std::vector<double>& qualVec, const std::vector<double>& errorVec,
    const std::vector<double>& weightVec, uint32_t qualStart, uint32_t qualStop,
    double qualStep, const std::string& model) {
  ownRInside ownSession;

  Rcpp::DataFrame inData = Rcpp::DataFrame::create(
      Rcpp::Named("qual") = qualVec, Rcpp::Named("error") = errorVec,
      Rcpp::Named("weight") = weightVec);
  auto& r = ownSession.get();
  r["inData"] = inData;
  r.parseEvalQ("myModel <- glm(error ~ qual, family = binomial(link = \"" +
               model + "\"), data = inData, weights= weight);");
  std::string seqStrCmd = to_string(qualStart) + "," + to_string(qualStop) +
                          ",by=" + to_string(qualStep);
  r.parseEvalQ("qTable = data.frame(qual = seq(" + seqStrCmd + "));");
  r.parseEvalQ(
      "qTable[, c(\"p\", \"se\")] <- predict(myModel, qTable, type = "
      "\"response\", se.fit = TRUE)[-3];");
  // R.parseEvalQ("print(qTable);");
  r.parseEvalQ("qTable_mat = as.matrix(qTable);");
  // R.parseEvalQ("cat(qTable_mat, sep =" ");");
  std::string txt = "qTable_mat;";
  Rcpp::NumericMatrix outData = r.parseEval(txt);
  std::unordered_map<std::string, std::vector<double>> out;

  std::vector<double> qualOut;
  for (auto i : iter::range(outData.nrow())) {
    qualOut.emplace_back(outData(i, 0));
  }
  out["qual"] = qualOut;

  std::vector<double> predOut;
  for (auto i : iter::range(outData.nrow())) {
    predOut.emplace_back(outData(i, 1));
  }
  out["pred"] = predOut;
  std::vector<double> seOut;
  for (auto i : iter::range(outData.nrow())) {
    seOut.emplace_back(outData(i, 2));
  }
  out["se"] = seOut;
  return out;
}

std::unordered_map<double, double> createLikelihoodMap(
    const std::unordered_map<std::string, std::vector<double>>& predOut) {
  std::vector<double> qualVec = predOut.at("qual");
  std::vector<double> predVec = predOut.at("pred");
  std::unordered_map<double, double> ans;
  for (const auto& i : iter::range(qualVec.size())) {
    ans[roundDecPlaces(qualVec[i], 2)] = predVec[i];
  }
  return ans;
}

std::unique_ptr<RInside> ownRInside::rSession_ = std::unique_ptr<RInside>();

ownRInside::ownRInside() {
  if (!rSession_) {
    rSession_.reset(new RInside());
  }
}

ownRInside::~ownRInside() {
  if (rSession_) {
    rSession_->parseEvalQ("rm(list = ls())");

    // according this, we can't terminate the actual R instance via "quit()"
    // http://r.789695.n4.nabble.com/Terminating-and-restarting-an-embedded-R-instance-possible-td4641823.html
    // so this is not allowed:
    // rSession_->parseEvalQ("quit()");

    // can't do this because Rf_endEmbeddedR(0) doesn't actually terminate R
    // rSession_.reset();
  }
}
void ownRInside::multipleParseEvalQ(const VecStr& cmds) {
  auto& r = get();
  for (const auto& cmd : cmds) {
    r.parseEvalQ(cmd);
  }
}
void ownRInside::openPdfDevice(const std::string & pdfFilename, double pageWidth, double pageHeight){
  std::stringstream txt;
  std::string widthStr = to_string(pageWidth);
  std::string heightStr = to_string(pageHeight);
	if(endsWith(pdfFilename, ".pdf")){
		txt << "pdf(\"" << pdfFilename << "\", width=" << widthStr << ", height=" << heightStr << ");";
	}else{
		txt << "pdf(\"" << pdfFilename << ".pdf" << "\", width=" << widthStr << ", height=" << heightStr << ");";
	}
  auto& r = get();
  r.parseEvalQ(txt.str());
}
void ownRInside::closeCurrentDevice(){
  std::string txt = "dev.off();";
  auto& r = get();

  r.parseEvalQ(txt);
}

}  // namespace bib
