#pragma once

#include <cassert>
#include <vector>

namespace km
{
/**
 * Generates the next combination of the integers [0,...,n-1] taken r at a time
 * in lexicographical order.
 *
 * The algorithm is from:
 * Mifsud, Charles J. "Algorithm 154: Combination in lexicographical order."
 * Communications of the ACM 6.3 (1963): 103.
 * http://dx.doi.org/10.1145/366274.366309
 */
class MifsudCombinations
{
 public:
  MifsudCombinations(unsigned int const n, unsigned int const r) : n_(n), r_(r)
  {
    assert(r <= n);
    vec_.reserve(r);
    for(unsigned int j = 0; j < r; ++j)
      vec_.push_back(j);
  }

  /**
   * Get the vector containing the combination.
   *
   * @return the vector which contains the current combination.
   */
  std::vector<unsigned int> const &get() const { return vec_; }

  /**
   * Generate the next combination
   *
   * @return true when more combinations are possible, false when done
   */
  bool next()
  {
    if(r_ == 0)
      return false;

    // Separate out the simple case for optimization
    if(vec_[r_ - 1] < n_ - 1)
    {
      vec_[r_ - 1] += 1;
      return true;
    }

    for(unsigned int i = 2; i <= r_; ++i)
    {
      unsigned int const idx = r_ - i; // Reverse index
      if(vec_[idx] < n_ - i)
      {
        vec_[idx] += 1;
        for(unsigned int j = 1; j < i; ++j)
          vec_[idx + j] = vec_[idx] + j;

        return true;
      }
    }
    return false;
  }

 private:
  std::vector<unsigned int> vec_;
  unsigned int const n_;
  unsigned int const r_;
};
}
