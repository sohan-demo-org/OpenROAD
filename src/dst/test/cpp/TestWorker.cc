#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>

#include "HelperCallBack.h"
#include "Worker.h"
#include "boost/asio.hpp"
#include "boost/asio/ip/tcp.hpp"
#include "boost/bind/bind.hpp"
#include "boost/system/system_error.hpp"
#include "boost/thread/thread.hpp"
#include "dst/Distributed.h"
#include "dst/JobMessage.h"
#include "gtest/gtest.h"
#include "utl/Logger.h"

namespace dst {

namespace {

bool isSocketPermissionError(const boost::system::system_error& error)
{
  return error.code()
             == boost::system::errc::make_error_code(
                 boost::system::errc::operation_not_permitted)
         || error.code()
                == boost::system::errc::make_error_code(
                    boost::system::errc::permission_denied);
}

bool isSocketStartupRuntimeError(const std::runtime_error& error)
{
  const std::string message = error.what();
  return message.starts_with("DST-0001");
}

bool isSocketEnvironmentErrorMessage(const std::string& message)
{
  return message.find("Operation not permitted") != std::string::npos
         || message.find("Permission denied") != std::string::npos
         || message.find("DST-0001") != std::string::npos;
}

}  // namespace

TEST(test_suite, test_worker)
{
  try {
    utl::Logger* logger = new utl::Logger();
    Distributed* dist = new Distributed(logger);
    std::string local_ip = "127.0.0.1";
    unsigned short port = 1234;

    Worker* worker = new Worker(dist, logger, local_ip.c_str(), port);
    boost::thread t(boost::bind(&Worker::run, worker));

    // Checking if the worker is up and calling callbacks correctly
    dist->addCallBack(new HelperCallBack(dist));
    JobMessage msg(JobMessage::JobType::kRouting);
    JobMessage result;
    EXPECT_TRUE(dist->sendJob(msg, local_ip.c_str(), port, result));
    EXPECT_EQ(result.getJobType(), JobMessage::JobType::kSuccess);
  } catch (const boost::system::system_error& error) {
    if (isSocketPermissionError(error)) {
      GTEST_SKIP() << error.what();
    }
    throw;
  } catch (const std::runtime_error& error) {
    if (isSocketStartupRuntimeError(error)) {
      GTEST_SKIP() << error.what();
    }
    throw;
  } catch (const std::exception& error) {
    if (isSocketEnvironmentErrorMessage(error.what())) {
      GTEST_SKIP() << error.what();
    }
    throw;
  }
}

namespace {

constexpr const char* kTestKey = "0123456789abcdef-test-key";

bool sendRoutingJob(Distributed* worker_dist,
                    Distributed* leader_dist,
                    unsigned short port)
{
  utl::Logger* logger = new utl::Logger();
  const std::string local_ip = "127.0.0.1";
  Worker* worker = new Worker(worker_dist, logger, local_ip.c_str(), port);
  boost::thread t(boost::bind(&Worker::run, worker));
  t.detach();
  worker_dist->addCallBack(new HelperCallBack(worker_dist));
  JobMessage msg(JobMessage::JobType::kRouting);
  JobMessage result;
  const bool ok = leader_dist->sendJob(msg, local_ip.c_str(), port, result);
  return ok && result.getJobType() == JobMessage::JobType::kSuccess;
}

}  // namespace

TEST(test_suite, test_worker_auth_matching_key)
{
  try {
    utl::Logger* logger = new utl::Logger();
    Distributed* worker_dist = new Distributed(logger);
    Distributed* leader_dist = new Distributed(logger);
    worker_dist->setAuthKey(kTestKey);
    leader_dist->setAuthKey(kTestKey);
    EXPECT_TRUE(sendRoutingJob(worker_dist, leader_dist, 1240));
  } catch (const std::exception& error) {
    if (isSocketEnvironmentErrorMessage(error.what())) {
      GTEST_SKIP() << error.what();
    }
    throw;
  }
}

TEST(test_suite, test_worker_rejects_unauthenticated_msg)
{
  try {
    utl::Logger* logger = new utl::Logger();
    Distributed* worker_dist = new Distributed(logger);
    Distributed* leader_dist = new Distributed(logger);
    worker_dist->setAuthKey(kTestKey);
    leader_dist->setAuthKey("");
    EXPECT_FALSE(sendRoutingJob(worker_dist, leader_dist, 1241));
  } catch (const std::exception& error) {
    if (isSocketEnvironmentErrorMessage(error.what())) {
      GTEST_SKIP() << error.what();
    }
    throw;
  }
}

TEST(test_suite, test_worker_rejects_wrong_key)
{
  try {
    utl::Logger* logger = new utl::Logger();
    Distributed* worker_dist = new Distributed(logger);
    Distributed* leader_dist = new Distributed(logger);
    worker_dist->setAuthKey(kTestKey);
    leader_dist->setAuthKey("fedcba9876543210-other-key");
    EXPECT_FALSE(sendRoutingJob(worker_dist, leader_dist, 1242));
  } catch (const std::exception& error) {
    if (isSocketEnvironmentErrorMessage(error.what())) {
      GTEST_SKIP() << error.what();
    }
    throw;
  }
}

TEST(test_suite, test_worker_refuses_non_loopback_without_key)
{
  utl::Logger* logger = new utl::Logger();
  Distributed* dist = new Distributed(logger);
  dist->setAuthKey("");
  try {
    dist->runWorker("0.0.0.0", 1243, true);
    FAIL() << "worker started on a non-loopback address without a key";
  } catch (const std::runtime_error& error) {
    EXPECT_TRUE(std::string(error.what()).starts_with("DST-0024"))
        << error.what();
  }
}

TEST(test_suite, test_worker_rejects_short_key)
{
  utl::Logger* logger = new utl::Logger();
  Distributed* dist = new Distributed(logger);
  setenv(Distributed::kAuthKeyEnvVar, "short", 1);
  try {
    dist->hasAuthKey();
    FAIL() << "short key accepted";
  } catch (const std::runtime_error& error) {
    EXPECT_TRUE(std::string(error.what()).starts_with("DST-0019"))
        << error.what();
  }
  unsetenv(Distributed::kAuthKeyEnvVar);
}

}  // namespace dst
