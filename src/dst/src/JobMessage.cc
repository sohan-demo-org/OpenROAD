// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "dst/JobMessage.h"

#include <cstddef>
#include <exception>
#include <sstream>
#include <string>
#include <string_view>

#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/hash2/sha2.hpp"
#include "boost/serialization/access.hpp"
#include "boost/serialization/unique_ptr.hpp"
#include "dst/BalancerJobDescription.h"

namespace dst {

namespace {

// Wire format of an authenticated message:
//   "DSTAUTH1 " <64 hex chars HMAC-SHA256(key, payload)> "\n" <payload>
constexpr std::string_view kAuthTag = "DSTAUTH1 ";
constexpr std::size_t kMacHexLength = 64;
constexpr std::size_t kAuthHeaderLength = kAuthTag.size() + kMacHexLength + 1;

std::string computeMac(const std::string& key, std::string_view payload)
{
  boost::hash2::hmac_sha2_256 hmac(key.data(), key.size());
  hmac.update(payload.data(), payload.size());
  const auto digest = hmac.result();
  static constexpr char kHex[] = "0123456789abcdef";
  std::string hex;
  hex.reserve(kMacHexLength);
  for (const unsigned char byte : digest) {
    hex.push_back(kHex[byte >> 4]);
    hex.push_back(kHex[byte & 0xF]);
  }
  return hex;
}

bool constantTimeEquals(std::string_view a, std::string_view b)
{
  if (a.size() != b.size()) {
    return false;
  }
  unsigned char diff = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    diff |= static_cast<unsigned char>(a[i] ^ b[i]);
  }
  return diff == 0;
}

// Splits an authenticated message into its payload after verifying the MAC.
bool verifyAndStrip(const std::string& key,
                    const std::string& str,
                    std::string_view& payload)
{
  if (str.size() < kAuthHeaderLength
      || std::string_view(str).substr(0, kAuthTag.size()) != kAuthTag
      || str[kAuthHeaderLength - 1] != '\n') {
    return false;
  }
  const std::string_view received_mac
      = std::string_view(str).substr(kAuthTag.size(), kMacHexLength);
  payload = std::string_view(str).substr(kAuthHeaderLength);
  return constantTimeEquals(received_mac, computeMac(key, payload));
}

}  // namespace

template <class Archive>
inline bool is_loading(const Archive& ar)
{
  return std::is_same_v<typename Archive::is_loading, boost::mpl::true_>;
}

template <class Archive>
void JobMessage::serialize(Archive& ar, const unsigned int version)
{
  (ar) & msg_type_;
  (ar) & job_type_;
  (ar) & desc_;
  if (!is_loading(ar)) {
    std::string eop = kEop;
    (ar) & eop;
  }
}

bool JobMessage::serializeMsg(SerializeType type,
                              JobMessage& msg,
                              std::string& str,
                              const std::string& auth_key)
{
  if (type == kWrite) {
    try {
      std::ostringstream oarchive_stream;
      boost::archive::text_oarchive archive(oarchive_stream);
      archive << msg;
      str = oarchive_stream.str();
    } catch (const boost::archive::archive_exception& e) {
      return false;
    }
    if (!auth_key.empty()) {
      str = std::string(kAuthTag) + computeMac(auth_key, str) + "\n" + str;
    }
  } else {
    std::string_view payload = str;
    if (!auth_key.empty() && !verifyAndStrip(auth_key, str, payload)) {
      return false;
    }
    try {
      std::istringstream iarchive_stream{std::string(payload)};
      boost::archive::text_iarchive archive(iarchive_stream);
      archive >> msg;
    } catch (const std::exception& e) {
      return false;
    }
  }
  return true;
}

}  // namespace dst
