#pragma once
//
//  bitSwaps.hpp
//  erpAnalysis
//
//  Created by Nick Hathaway on 11/9/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include <stdint.h>

//! Byte swap unsigned short
namespace bibseq {

inline uint16_t swap_uint16(uint16_t val) {
	return (val << 8) | (val >> 8);
}

//! Byte swap short
inline int16_t swap_int16(int16_t val) {
  return (val << 8) | ((val >> 8) & 0xFF);
}

//! Byte swap uint32_t
inline uint32_t swap_uint32(uint32_t val) {
  val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF);
  return (val << 16) | (val >> 16);
}

//! Byte swap int
inline int32_t swap_int32(int32_t val) {
  val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF);
  return (val << 16) | ((val >> 16) & 0xFFFF);
}

inline int64_t swap_int64(int64_t val) {
  val = ((val << 8) & 0xFF00FF00FF00FF00ULL) |
        ((val >> 8) & 0x00FF00FF00FF00FFULL);
  val = ((val << 16) & 0xFFFF0000FFFF0000ULL) |
        ((val >> 16) & 0x0000FFFF0000FFFFULL);
  return (val << 32) | ((val >> 32) & 0xFFFFFFFFULL);
}

inline uint64_t swap_uint64(uint64_t val) {
  val = ((val << 8) & 0xFF00FF00FF00FF00ULL) |
        ((val >> 8) & 0x00FF00FF00FF00FFULL);
  val = ((val << 16) & 0xFFFF0000FFFF0000ULL) |
        ((val >> 16) & 0x0000FFFF0000FFFFULL);
  return (val << 32) | (val >> 32);
}
}  // namespace bib