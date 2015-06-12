#include "runningParameters.hpp"
#include "bibseq/objects/dataContainers/table.hpp"

#include <bibcpp/bashUtils.h>

namespace bibseq {

namespace cto = bib::bashCT;
runningParameters::runningParameters(const std::vector<double>& parameter, uint32_t iterNumber,
		uint32_t clusterCount, bool onPerId) {
	if (onPerId) {
		if(parameter.size() != 3){
			std::stringstream ss;
			ss << cto::boldRed("Error in constructing runningParameters, size should be 3, not " )<< parameter.size() << std::endl;
			std::runtime_error{ss.str()};
		}
	  stopCheck_ = parameter[0];
	  if (stopCheck_ < 1.00 && stopCheck_ > 0.00) {
	    stopCheck_ = (int)(stopCheck_ * clusterCount);
	  } else if (stopCheck_ < 0.00) {
	    stopCheck_ = clusterCount;
	  }
	  smallCheckStop_ = parameter[1];
	  iterNumber_ = iterNumber;
	  // error parameters
	  errors_.distances_.percentIdentity_ = parameter[2];
	} else {
		if(parameter.size() != 8){
			std::stringstream ss;
			ss << cto::boldRed("Error in constructing runningParameters, size should be 8, not " )<< parameter.size() << std::endl;
			std::runtime_error{ss.str()};
		}
	  stopCheck_ = parameter[0];
	  if (stopCheck_ < 1.00 && stopCheck_ > 0.00) {
	    stopCheck_ = (int)(stopCheck_ * clusterCount);
	  } else if (stopCheck_ < 0.00) {
	    stopCheck_ = clusterCount;
	  }
	  smallCheckStop_ = parameter[1];
	  iterNumber_ = iterNumber;
	  // error parameters
	  errors_.oneBaseIndel_ = parameter[2];
	  errors_.twoBaseIndel_ = parameter[3];
	  errors_.largeBaseIndel_ = parameter[4];
	  errors_.hqMismatches_ = parameter[5];
	  errors_.lqMismatches_ = parameter[6];
	  errors_.lowKmerMismatches_ = parameter[7];
	}


}



std::string centerText(const std::string& text, uint32_t maxWidth) {
	uint32_t spaceWidth = maxWidth - text.size();
	uint32_t front = round(spaceWidth/2.0);
	uint32_t back = spaceWidth - front;
  return std::string(front, ' ') + text + std::string(back, ' ');
}

void runningParameters::printIterInfo(std::ostream & out, bool colorFormat, bool onPerId){
	if(onPerId){
		std::cout << cto::bold
				<< cto::addBGColor(145) << " Iteration "
				<< cto::addBGColor(188) << " StopAfter "
				<< cto::addBGColor(145) << " SizeCutOff "
				<< cto::addBGColor(188) << " PerIdCutOff "<< cto::reset << std::endl;
		std::cout
		<< cto::addBGColor(145) << " " << centerText(to_string(iterNumber_), 9) << " "
		<< cto::addBGColor(188) << " " << centerText(to_string(stopCheck_), 9) << " "
		<< cto::addBGColor(145) << " " << centerText(to_string(smallCheckStop_), 10) << " "
		<< cto::addBGColor(188) << " " << centerText(to_string(errors_.distances_.percentIdentity_ * 100) + "%", 11) << " "
		<< cto::reset << std::endl;
	}else{
		std::cout << cto::bold
				<< cto::addBGColor(145) << " Iteration "
				<< cto::addBGColor(188) << " StopAfter "
				<< cto::addBGColor(145) << " SizeCutOff "
				<< cto::addBGColor(188) << " 1BaseIndels "
				<< cto::addBGColor(145) << " 2BaseIndels "
				<< cto::addBGColor(188) << " >2BaseIndels "
				<< cto::addBGColor(145) << " HQMismatches "
				<< cto::addBGColor(188) << " LQMismatches "
				<< cto::addBGColor(145) << " LKMismatches "  << cto::reset << std::endl;
		std::cout
		<< cto::addBGColor(145) << " " << centerText(to_string(iterNumber_), 9) << " "
		<< cto::addBGColor(188) << " " << centerText(to_string(stopCheck_), 9) << " "
		<< cto::addBGColor(145) << " " << centerText(to_string(smallCheckStop_), 10) << " "
		<< cto::addBGColor(188) << " " << centerText(to_string(errors_.oneBaseIndel_), 11) << " "
		<< cto::addBGColor(145) << " " << centerText(to_string(errors_.twoBaseIndel_), 11) << " "
		<< cto::addBGColor(188) << " " << centerText(to_string(errors_.largeBaseIndel_), 12) << " "
		<< cto::addBGColor(145) << " " << centerText(to_string(errors_.hqMismatches_), 12) << " "
		<< cto::addBGColor(188) << " " << centerText(to_string(errors_.lqMismatches_), 12) << " "
		<< cto::addBGColor(145) << " " << centerText(to_string(errors_.lowKmerMismatches_), 12) << " "
		<< cto::reset << std::endl;
	}

}


std::map<int32_t, std::vector<double>> runningParameters::processParameters(const std::string & parametersFilename){
	std::map<int32_t, std::vector<double>> ret;
  table inTab(parametersFilename, ":");
  std::map<int, std::vector<double>>::iterator iteratorMapIter;
  int32_t iters = 1;
  for (const auto& row : inTab.content_) {
    if (row.empty() || row.front().front() == 's' || row.front().front() == '#') {
      continue;
    }
    if(row.size() != 8){
    	std::stringstream ss;
    	ss << "Error in parsing parameters file: " << parametersFilename << " row should be width of 8 not " << row.size() << std::endl;

    	throw std::runtime_error{cto::boldRed(ss.str())};
    }
    std::vector<double> tempVect;
    for (const auto & colPos : iter::range(row.size() )) {
      if (colPos == 0) {
        if (row[colPos].find("%") != std::string::npos) {
          double percent = std::stod(
          		row[colPos].substr(0, row[colPos].size() - 1));
          if (percent >= 100) {
            tempVect.emplace_back(0.99);
          } else if (percent > 0) {
            tempVect.emplace_back(percent / 100.00);
          } else {
            tempVect.emplace_back(0.10);
          }
        } else if (stringToLowerReturn(row[colPos]) == "all") {
          tempVect.emplace_back(-1.00);
        } else {
          tempVect.emplace_back(std::stod(row[colPos]));
        }
      } else {
      	tempVect.emplace_back(std::stod(row[colPos]));
      }
    }
    ret.insert(std::make_pair(iters, tempVect));
    iters++;
  }
  return ret;
}

std::map<int32_t, std::vector<double>> runningParameters::processParametersPerId(const std::string & parametersFilename){
	std::map<int32_t, std::vector<double>> ret;
  table inTab(parametersFilename, ":");
  std::map<int, std::vector<double>>::iterator iteratorMapIter;
  int32_t iters = 1;
  for (const auto& row : inTab.content_) {
    if (row.empty()  || row.front().front() == 's' || row.front().front() == '#') {
      continue;
    }
    if(row.size() != 3){
    	std::stringstream ss;
    	ss << "Error in parsing parameters file: " << parametersFilename << " row should be width of 3 not " << row.size() << std::endl;

    	throw std::runtime_error{cto::boldRed(ss.str())};
    }
    std::vector<double> tempVect;
  	//printVector(row);
    for (const auto & colPos : iter::range(row.size() )) {

      if (colPos == 0) {
        if (row[colPos].find("%") != std::string::npos) {
          double percent = std::stod(
          		row[colPos].substr(0, row[colPos].size() - 1));
          if (percent >= 100) {
            tempVect.emplace_back(0.99);
          } else if (percent > 0) {
            tempVect.emplace_back(percent / 100.00);
          } else {
            tempVect.emplace_back(0.10);
          }
        } else if (stringToLowerReturn(row[colPos]) == "all") {
          tempVect.emplace_back(-1.00);
        } else {
          tempVect.emplace_back(std::stod(row[colPos]));
        }
      } else {
      	tempVect.emplace_back(std::stod(row[colPos]));
      }
    }
    ret.insert(std::make_pair(iters, tempVect));
    iters++;
  }
  return ret;
}

}  // namespace bib
