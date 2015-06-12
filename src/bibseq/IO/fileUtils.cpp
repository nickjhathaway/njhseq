#include "fileUtils.hpp"
#include "bibseq/utils/utils.hpp"
#include "bibseq/utils/stringUtils.hpp"
#include <bibcpp/files/fileUtilities.hpp>

namespace bibseq {


VecStr getNewestDirs(const std::string & dirName, const std::string & con){
  VecStr ret;
  auto dirs = bib::files::listAllFiles(dirName, false, {con});

  std::map<std::string, std::map<std::string,bib::files::bfs::path>> filesByTime;
  for(const auto & d : dirs){
  	if(d.second){
  		auto toks = tokenizeString(d.first.string(), "/");
  		auto fileToks = tokenizeString(toks.back(), con);
  		filesByTime[fileToks.front()][fileToks.back()] = d.first;
  	}
  }
  for(const auto & f : filesByTime){
  	ret.emplace_back(f.second.rbegin()->second.string());
  }
  return ret;
}




int getdir(const std::string &dir,
           std::map<std::string, std::pair<std::string, bool>> &files) {
  DIR *dp;
  struct dirent *dirp;
  if ((dp = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }
  int status = 0;
  while ((dirp = readdir(dp)) != NULL) {

    struct stat st_buf;
    std::stringstream tempName;
    tempName << dir << "/" << std::string(dirp->d_name).c_str();
    // status=stat(std::string(dirp->d_name).c_str(),&st_buf);
    status = stat(tempName.str().c_str(), &st_buf);
    if (std::string(dirp->d_name) == "." || std::string(dirp->d_name) == "..") {
    } else {
      if (S_ISDIR(st_buf.st_mode)) {
        files.insert(
            std::make_pair(tempName.str(), std::make_pair("directory", false)));
      } else {
        files.insert(
            std::make_pair(tempName.str(), std::make_pair("file", false)));
      }
    }
  }
  closedir(dp);
  return status;
}

std::map<std::string, std::pair<std::string, bool>> listFilesInDir(
    const std::string &directoryName, bool recursive) {
  std::map<std::string, std::pair<std::string, bool>> files;
  std::map<std::string, std::pair<std::string, bool>>::iterator fileIter;
  getdir(directoryName, files);
  if (recursive) {
    bool searching = true;
    while (searching) {
      for (fileIter = files.begin(); fileIter != files.end(); ++fileIter) {
        if (fileIter->second.first == "directory" && !fileIter->second.second) {
          getdir(fileIter->first, files);
          fileIter->second.second = true;
        }
      }
      for (fileIter = files.begin(); fileIter != files.end(); ++fileIter) {
        searching = false;
        if (fileIter->second.first == "directory" && !fileIter->second.second) {
          searching = true;
          break;
        }
      }
    }
  }
  return files;
}
// check to see if a file exists
bool fexists(const std::string &filename) {
  /*boost::filesystem::path p(filename);
  if(boost::filesystem::exists(p)){
          std::cout << "exists" << std::endl;
  }*/
  std::ifstream ifile(filename.c_str());
  if (ifile) {
    return true;
  } else {
    return false;
  }
}

void openTextFile(std::ofstream &file, std::string filename,
                  std::string fileExtention, bool overWrite,
                  bool exitOnFailure) {

  if (filename.find(fileExtention) == std::string::npos) {
    filename.append(fileExtention);
  } else if (filename.substr(filename.size() - fileExtention.size(),
                             fileExtention.size()) != fileExtention) {
    filename.append(fileExtention);
  }
  bib::files::openTextFile(file, filename, overWrite, false, exitOnFailure);
}
// runLog stuff
void startRunLog(std::ofstream &runLog, const MapStrStr &inputCommands) {
  runLog << "Ran " << getCurrentDate() << std::endl;
  runLog <<  inputCommands.find("-commandline")->second << std::endl;
}


std::map<std::string, std::pair<std::string, bool>> getFiles(
    const std::string &directoryName, const std::string &contains,
    const std::string &filesOrDirectories, bool specific, bool recursive) {
  std::map<std::string, std::pair<std::string, bool>> allFiles =
      listFilesInDir(directoryName, recursive);
  std::map<std::string, std::pair<std::string, bool>>::iterator fileIter;

  if (filesOrDirectories == "file") {
    std::map<std::string, std::pair<std::string, bool>> allFileFiles;
    for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
      if (fileIter->second.first == filesOrDirectories) {
        allFileFiles.insert(*fileIter);
      }
    }
    allFiles = allFileFiles;
  } else if (filesOrDirectories == "directory") {
    std::map<std::string, std::pair<std::string, bool>> allDirectoryFiles;
    for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
      if (fileIter->second.first == filesOrDirectories) {
        allDirectoryFiles.insert(*fileIter);
      }
    }
    allFiles = allDirectoryFiles;
  }

  if (specific) {
    std::map<std::string, std::pair<std::string, bool>> specificFiles;
    for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
      if (fileIter->first.find(contains) != std::string::npos) {
        specificFiles.insert(*fileIter);
      }
    }
    return specificFiles;
  } else {
    return allFiles;
  }
}
std::map<std::string, std::pair<std::string, bool>> getFiles(
    const std::string &directoryName, const VecStr &contains,
    const std::string &filesOrDirectories, bool specific, bool recursive) {
  std::map<std::string, std::pair<std::string, bool>> allFiles =
      listFilesInDir(directoryName, recursive);
  std::map<std::string, std::pair<std::string, bool>>::iterator fileIter;
  if (filesOrDirectories == "file") {
    std::map<std::string, std::pair<std::string, bool>> allFileFiles;
    for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
      if (fileIter->second.first == filesOrDirectories) {
        allFileFiles.insert(*fileIter);
      }
    }
    allFiles = allFileFiles;
  } else if (filesOrDirectories == "directory") {
    std::map<std::string, std::pair<std::string, bool>> allDirectoryFiles;
    for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
      if (fileIter->second.first == filesOrDirectories) {
        allDirectoryFiles.insert(*fileIter);
      }
    }
    allFiles = allDirectoryFiles;
  }

  if (specific) {
    std::map<std::string, std::pair<std::string, bool>> specificFiles;
    for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
      bool containsAll = true;
      for (VecStr::const_iterator conIter = contains.begin();
           conIter != contains.end(); ++conIter) {
        if (fileIter->first.find(*conIter) == std::string::npos) {
          containsAll = false;
          break;
        }
      }
      if (containsAll) {
        specificFiles.insert(*fileIter);
      }
    }
    return specificFiles;
  } else {
    return allFiles;
  }
}

VecStr collectFilenames(VecStr contains, VecStr doesNotContain,
                        bool recursive) {
  auto files = listFilesInDir(".", recursive);
  VecStr specificFiles;
  for (const auto &f : files) {
    bool include = true;
    for (const auto &con : contains) {
      if (f.first.find(con) == std::string::npos) {
        include = false;
        break;
      }
    }
    for (const auto &nott : contains) {
      if (f.first.find(nott) != std::string::npos) {
        include = false;
        break;
      }
    }
    if (include) {
      specificFiles.push_back(f.first);
    }
  }
  return specificFiles;
}
/*
VecStr collectFilenamesWithExtention(const std::string &directory,
                                     const std::string &extention,
                                     bool recursive) {
  auto files = listFilesInDir(directory, recursive);
  VecStr specificFiles;
  for (const auto &f : files) {
    if (getExtention(f.first) == extention) {
      specificFiles.push_back(f.first);
    }
  }
  /
  for (const auto &f : files) {
    bool include = true;
    for (const auto &con : contains) {
      if (f.first.find(con) == std::string::npos) {
        include = false;
        break;
      }
    }
    for (const auto &nott : contains) {
      if (f.first.find(nott) != std::string::npos) {
        include = false;
        break;
      }
    }
    if (include) {
      specificFiles.push_back(f.first);
    }
  }/
  return specificFiles;
}*/
std::string cleanOutPerLine(const std::string &in, uint32_t width,
                            uint32_t indentAmount) {
  std::stringstream tempStream(in);
  std::stringstream out;
  uint32_t characterLengthSum = 0;
  std::string str;
  tempStream >> str;
  out << str << " ";
  characterLengthSum += str.length() + 1;
  while (!tempStream.eof()) {
    tempStream >> str;
    if ((characterLengthSum + str.length() + 1) > width) {
      out << std::endl << std::string(indentAmount, ' ') << str << " ";
      characterLengthSum = indentAmount;
    } else {
      out << str << " ";
    }
    characterLengthSum += str.length() + 1;
  }
  return out.str();
}

std::string cleanOut(const std::string &in, uint32_t width,
                     uint32_t indentAmount) {
  std::istringstream inStream(in);
  std::stringstream out;
  for (std::string line; std::getline(inStream, line);) {
    out << cleanOutPerLine(line, width, indentAmount) << std::endl;
  }
  return out.str();
}


}  // namespace bib
