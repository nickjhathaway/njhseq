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
//
//  main.cpp
//
//  Created by Nicholas Hathaway on 8/11/13.
//


#include "njhseq.h"

int main(int argc, char* argv[]) {
	try {
	  std::cout << "Hello njhseq" << std::endl;
	  return 0;
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}

}
