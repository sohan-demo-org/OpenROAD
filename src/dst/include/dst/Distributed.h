// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "boost/asio.hpp"

namespace utl {
class Logger;
}

namespace dst {

namespace asio = boost::asio;
using asio::ip::tcp;
using Socket = asio::basic_stream_socket<tcp>;

class JobMessage;
class JobCallBack;
class Worker;

class Distributed
{
 public:
  // Name of the environment variables holding the shared secret used to
  // authenticate every DST message (HMAC-SHA256). The key can be given
  // directly or as a path to a file containing it.
  static constexpr const char* kAuthKeyEnvVar = "OPENROAD_DST_AUTH_KEY";
  static constexpr const char* kAuthKeyFileEnvVar
      = "OPENROAD_DST_AUTH_KEY_FILE";
  static constexpr std::size_t kMinAuthKeyLength = 16;
  // Upper bound on a single received message; larger ones are dropped.
  static constexpr std::size_t kMaxMessageSize = std::size_t{1} << 30;

  Distributed(utl::Logger* logger);
  ~Distributed();
  void runWorker(const char* ip, unsigned short port, bool interactive);
  void setAuthKey(const std::string& key);
  const std::string& getAuthKey();
  bool hasAuthKey() { return !getAuthKey().empty(); }
  // Directory that every file path received by a worker must be inside.
  void setSharedVolume(const std::string& dir) { shared_volume_ = dir; }
  const std::string& getSharedVolume() const { return shared_volume_; }
  void runLoadBalancer(const char* ip,
                       unsigned short port,
                       const char* workers_domain);
  void addWorkerAddress(const char* address, unsigned short port);
  bool sendJob(JobMessage& msg,
               const char* ip,
               unsigned short port,
               JobMessage& result);
  bool sendJobMultiResult(JobMessage& msg,
                          const char* ip,
                          unsigned short port,
                          JobMessage& result);
  bool sendResult(JobMessage& msg, Socket& sock);
  void addCallBack(JobCallBack* cb);
  const std::vector<JobCallBack*>& getCallBacks() const { return callbacks_; }

 private:
  struct EndPoint
  {
    std::string ip;
    uint16_t port;
    EndPoint(const std::string& ip_in, uint16_t port_in)
        : ip(ip_in), port(port_in)
    {
    }
  };
  void loadAuthKeyFromEnv();
  void checkBindAddress(const char* ip, const char* role);

  utl::Logger* logger_;
  std::once_flag auth_key_loaded_;
  std::string auth_key_;
  std::string shared_volume_;
  std::vector<EndPoint> end_points_;
  std::vector<JobCallBack*> callbacks_;
  std::vector<std::unique_ptr<Worker>> workers_;
};
}  // namespace dst
