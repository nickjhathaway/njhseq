#pragma once
/*
 * PopNamesInfo.hpp
 *
 *  Created on: Jul 31, 2016
 *      Author: nick
 */

#include "bibseq/utils.h"

namespace bibseq {

/**@brief simple class to hold the population name and the samples that go in it
 *
 */
class PopNamesInfo {

public:
	/**@brief construct with pop name and samples
	 *
	 * @param populationName population name
	 * @param samples the samples
	 */
	PopNamesInfo(std::string populationName, std::set<std::string> samples);
	/**@brief construct with population name and samples
	 *
	 * @param populationName the population name
	 * @param samples samples in vector that will be converted into set to get rid of duplicate names
	 */
	PopNamesInfo(std::string populationName, VecStr samples);

	std::string populationName_; /**< population name*/
	std::set<std::string> samples_; /**< samples in set so no duplicates*/

	/**@brief check to see if sample is samples_
	 *
	 * @param sample the sample to check for
	 * @return true if sample is in samples_
	 */
	bool hasSample(const std::string & sample) const;

	/**@brief output json object representation of class
	 *
	 * @return Json::Value with info of class
	 */
	Json::Value toJson() const;
private:
	/**@brief check formating of population name, it cannot contain ".", will throw if fails format requirments
	 *
	 */
	void checkPopNameThrow() const;



};


}  // namespace bibseq






