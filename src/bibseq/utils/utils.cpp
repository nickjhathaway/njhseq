#include "utils.hpp"
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//////durations and date
namespace bibseq {
std::map<uint64_t, uint32_t> printStringLengths(const VecStr &strings,
                                                std::ostream &out) {
  std::map<uint64_t, uint32_t> lengths;
  for (const auto &rs : strings) {
    ++lengths[rs.length()];
  }
  for (const auto &len : lengths) {
    out << len.first << ":" << len.second << std::endl;
  }
  return lengths;
}

std::string getCurrentDate() {
  time_t t = time(0);  // get time now
  struct tm *now = localtime(&t);
  std::stringstream timeStream;
  timeStream << (now->tm_year + 1900) << '-' << leftPadNumStr((now->tm_mon + 1))
             << '-' << leftPadNumStr(now->tm_mday) << '_'
             << leftPadNumStr(now->tm_hour) << '.'
             << leftPadNumStr(now->tm_min);
  return timeStream.str();
}

std::string convertBoolToString(bool convert) {
  if (convert) {
    return "true";
  } else {
    return "false";
  }
}


// with no header
void printTableOrganized(const std::vector<VecStr> &content,
                         std::ostream &out) {
  std::map<int, size_t> sizeMap;
  for (const auto &contentIter : content) {
    int count = 0;
    for (const auto &lineIter : contentIter) {
      if (sizeMap.find(count) == sizeMap.end()) {
        sizeMap.insert(std::make_pair(count, lineIter.size()));
      } else {
        if (sizeMap[count] < lineIter.size()) {
          sizeMap[count] = lineIter.size();
        }
      }
      ++count;
    }
  }
  for (const auto &contentIter : content) {
    int count = 0;
    for (const auto &lineIter : contentIter) {
      out << std::setw((int)sizeMap[count]) << std::left << lineIter;
      out << "\t";
      ++count;
    }
    out << std::endl;
  }
}
// with header
void printTableOrganized(const std::vector<VecStr> &content,
                         const VecStr &header, std::ostream &out) {
  std::map<int, size_t> sizeMap;
  {
    int count = 0;
    for (const auto &lineIter : header) {
      if (sizeMap.find(count) == sizeMap.end()) {
        sizeMap.insert(std::make_pair(count, lineIter.size()));
      } else {
        if (sizeMap[count] < lineIter.size()) {
          sizeMap[count] = lineIter.size();
        }
      }
      ++count;
    }
  }
  for (const auto &contentIter : content) {
    int count = 0;
    for (const auto &lineIter : contentIter) {
      if (sizeMap.find(count) == sizeMap.end()) {
        sizeMap.insert(std::make_pair(count, lineIter.size()));
      } else {
        if (sizeMap[count] < lineIter.size()) {
          sizeMap[count] = lineIter.size();
        }
      }
      ++count;
    }
  }
  {
    int count = 0;
    for (const auto &lineIter : header) {
      out << std::setw((int)sizeMap[count]) << std::left << lineIter;
      out << "\t";
      ++count;
    }
    out << std::endl;
  }
  for (const auto &contentIter : content) {
    int count = 0;
    for (const auto &lineIter : contentIter) {
      out << std::setw((int)sizeMap[count]) << std::left << lineIter;
      out << "\t";
      ++count;
    }
    out << std::endl;
  }
}


// from
// http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
size_t getPeakRSS() {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.PeakWorkingSetSize;

#elif(defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) ||       \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
  /* AIX and Solaris ------------------------------------------ */
  struct psinfo psinfo;
  int fd = -1;
  if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
    return (size_t)0L; /* Can't open? */
  if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo)) {
    close(fd);
    return (size_t)0L; /* Can't read? */
  }
  close(fd);
  return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || \
    (defined(__APPLE__) && defined(__MACH__))
  /* BSD, Linux, and OSX -------------------------------------- */
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
  return (size_t)rusage.ru_maxrss;
#else
  return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
  /* Unknown OS ----------------------------------------------- */
  return (size_t)0L; /* Unsupported. */
#endif
}
size_t getCurrentRSS() {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
  /* OSX ------------------------------------------------------ */
  struct mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t) & info,
                &infoCount) != KERN_SUCCESS)
    return (size_t)0L; /* Can't access? */
  return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
  /* Linux ---------------------------------------------------- */
  long rss = 0L;
  FILE *fp = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL)
    return (size_t)0L; /* Can't open? */
  if (fscanf(fp, "%*s%ld", &rss) != 1) {
    fclose(fp);
    return (size_t)0L; /* Can't read? */
  }
  fclose(fp);
  return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
  /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
  return (size_t)0L; /* Unsupported. */
#endif
}
}  // namespace bib
