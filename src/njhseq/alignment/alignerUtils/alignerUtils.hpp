#pragma once
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
/*
 * alignerUtils.hpp
 *
 *  Created on: Nov 30, 2015
 *      Author: nick
 */


#include "njhseq/common.h"

namespace njhseq {

/**@brief Get The alignment position for the actual real seq position, accounting for gaps
 *
 * This will throw if seq size is 0 or if realSeqPos is greater than seq size or greater than a position real seq size
 *
 * @param seq The alignment with gpas
 * @param realSeqPos the actual seq position
 * @return The alignment position that corresponds to the actual position
 */
size_t getAlnPosForRealPos(const std::string & seq, size_t realSeqPos);

/**@brief Determine the actual position for the given alignment position accounting for gaps
 *
 * This will throw if seq size is 0 or if seqAlnPos is greater than or equal to seq size
 *
 * @param seq The alignment with gaps
 * @param seqAlnPos The alignment position
 * @return The actual position, if seqAlnPos is a front gap will return size_t max (std::string::npos), if seqAlnPos is in an end gap will return the actual seq size
 */
size_t getRealPosForAlnPos(const std::string & seq, size_t seqAlnPos);

/**@brief Get the actual position for the query sequence for the reference position
 *
 * @param ref The alignment string for ref
 * @param seq The alignment string for query seq
 * @param refPos The position in the reference ref (unaligned position)
 * @return the unaligned position of the query seq for the corresponding seq position
 */
size_t getRealPosForSeqForRefPos(const std::string & ref, const std::string & seq, size_t refPos);
/**@brief Get the actual position for the ref sequence for the query sequence reference position
 *
 * @param ref The alignment string for ref
 * @param seq The alignment string for query seq
 * @param seqPos The position in the query seq (unaligned position)
 * @return the unaligned position of the reference seq for the corresponding query seq position
 */
size_t getRealPosForRefForSeqPos(const std::string & ref, const std::string & seq, size_t seqPos);

}  // namespace njhseq

