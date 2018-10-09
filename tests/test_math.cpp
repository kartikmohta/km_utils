#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <km_utils/math.hpp>

using namespace km::types;
using namespace km::math;

TEST_CASE("powInt")
{
  CHECK(powInt(1, 100) == 1);
  CHECK(powInt(2, 0) == 1);
  CHECK(powInt(2, 1) == 2);
  CHECK(powInt(2, 2) == 4);
  CHECK(powInt(2, 3) == 8);
  CHECK(powInt(2, 4) == 16);
  CHECK(powInt(2, 5) == 32);
  CHECK(powInt(2, 6) == 64);
  CHECK(powInt(2, 7) == 128);
  CHECK(powInt(2, 8) == 256);
  CHECK(powInt(2, 9) == 512);
  CHECK(powInt(2, 10) == 1024);
  CHECK(powInt(10, 1) == 10);
  CHECK(powInt(10, 2) == 100);
  CHECK(powInt(10, 3) == 1000);
  CHECK(powInt(10, 4) == 10000);
  CHECK(powInt(10, 5) == 100000);
  CHECK(powInt(2.0, 10) == doctest::Approx(1024));
  CHECK(powInt(3.0, 10) == doctest::Approx(59049));
}

TEST_CASE("matrixPermanent")
{
  double expected_val[] = {0, 1, 2, 9, 44, 265, 1854, 14833, 133496, 1334961};
  {
    auto M = SquareMat<double, 1>::Zero().eval();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[0]));
  }
  {
    auto M = SquareMat<double, 2>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[1]));
  }
  {
    auto M = SquareMat<double, 3>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[2]));
  }
  {
    auto M = SquareMat<double, 4>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[3]));
  }
  {
    auto M = SquareMat<double, 5>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[4]));
  }
  {
    auto M = SquareMat<double, 6>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[5]));
  }
  {
    auto M = SquareMat<double, 7>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[6]));
  }
  {
    auto M = SquareMat<double, 8>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[7]));
  }
  {
    auto M = SquareMat<double, 9>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[8]));
  }
  {
    auto M = SquareMat<double, 10>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(expected_val[9]));
  }
  {
    auto M = SquareMat<double, 15>::Ones().eval();
    M.diagonal().setZero();
    CHECK(matrixPermanent(M) == doctest::Approx(481066515734.0));
  }
}
