#pragma once
//
//  readers.h
//
//  Created by Nicholas Hathaway on 7/20/13.
//
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
#include "njhseq/IO/SeqIO.h"
#include "njhseq/IO/cachedReader.hpp"
#include "njhseq/IO/IOUtils.hpp"
#include "njhseq/IO/fileUtils.hpp"
#include "njhseq/IO/FileWithTime.hpp"
#include "njhseq/IO/IOOptions.h"
#include "njhseq/IO/OutputStream.hpp"
#include "njhseq/IO/InputStream.hpp"

#include "njhseq/IO/MultiOutputStream.hpp"
#include "njhseq/IO/MultiOutputStreamCache.hpp"
