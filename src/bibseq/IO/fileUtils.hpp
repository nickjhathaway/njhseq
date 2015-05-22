#pragma once
//
//  fileUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/15/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include <map>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include "bibseq/common/typedefs.hpp"

namespace bibseq {



int getdir(const std::string &dir,
           std::map<std::string, std::pair<std::string, bool>> &files);

std::map<std::string, std::pair<std::string, bool>> listFilesInDir(
    const std::string &directoryName, bool recursive);

// check to see if a file exists
bool fexists(const std::string &filename);

void openTextFile(std::ofstream &file, std::string filename,
                  std::string fileExtention, bool overWrite,
                  bool exitOnFailure);
template<typename OPTIONS>
void openTextFile(std::ofstream &file, std::string filename,
                  std::string fileExtention, OPTIONS outOptions){
	openTextFile(file, filename, fileExtention,
			outOptions.overWriteFile_, outOptions.exitOnFailureToWrite_);
}

// runLog stuff
void startRunLog(std::ofstream &runLog, const MapStrStr &inputCommands);

std::map<std::string, std::pair<std::string, bool>> getFiles(
    const std::string &directoryName, const std::string &contains,
    const std::string &filesOrDirectories, bool specific, bool recursive);

std::map<std::string, std::pair<std::string, bool>> getFiles(
    const std::string &directoryName, const VecStr &contains,
    const std::string &filesOrDirectories, bool specific, bool recursive);

VecStr collectFilenames(VecStr contains, VecStr doesNotContain,
                        bool recursive = true);
VecStr collectFilenamesWithExtention(const std::string &directory,
                                     const std::string &extention,
                                     bool recursive = true);
std::string cleanOutPerLine(const std::string &in, uint32_t width,
                            uint32_t indentAmount);
std::string cleanOut(const std::string &in, uint32_t width,
                     uint32_t indentAmount);




}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "fileUtils.cpp"
#endif