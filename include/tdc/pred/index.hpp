#pragma once

#include <cstdint>
#include <cstddef>

#include <tdc/vec/int_vector.hpp>
#include <tdc/math/ilog2.hpp>
#include <tdc/util/assert.hpp>
#include <tdc/util/likely.hpp>

#include "binary_search_hybrid.hpp"
#include "result.hpp"

namespace tdc {
namespace pred {

/// \brief Predecessor search using universe-based sampling.
///
/// Using this approach, we define a parameter <tt>k</tt> so that for a predecessor query, we can look up an interval of size at most <tt>2^k</tt> and
/// proceed with a \ref BinarySearchHybrid in that interval.
template<typename t_word>
class Index {
private:
    inline t_word hi(t_word x) const {
        return x >> m_lo_bits;
    }
    
    size_t m_lo_bits;
    t_word m_min, m_max;
    t_word m_key_min, m_key_max;

    vec::IntVector m_hi_idx;
    
public:
    inline Index() : m_lo_bits(0), m_min(t_word{0}-1), m_max(0), m_key_min(t_word{0}-1), m_key_max(0) {
    }

    /// \brief Constructs the index for the given keys.
    /// \param keys a pointer to the keys, that must be in ascending order
    /// \param num the number of keys
    /// \param lo_bits the number of low key bits, defining the maximum size of a search interval; lower means faster queries, but more memory usage
    Index(const t_word* keys, const size_t num, const size_t lo_bits) : m_lo_bits(lo_bits) {
    assert_sorted_ascending(keys, num);
    // build an index for high bits
    if(num == 0) {
        m_lo_bits = 0;
        m_min = t_word{0}-1;
        m_max = 0;
        m_key_min = t_word{0}-1;
        m_key_max = 0;
        return;
    }
    m_min = keys[0];
    m_max = keys[num-1];
    m_key_min = hi(m_min);
    m_key_max = hi(m_max);
    
    // std::cout << "m_min=" << m_min << ", m_key_min=" << m_key_min;
    // std::cout << ", m_max=" << m_max << ", m_key_max=" << m_key_max;
    // std::cout << " -> allocating " << (m_key_max - m_key_min + 2)
        // << " samples with bit width " << math::ilog2_ceil(num-1) << std::endl;

    m_hi_idx = vec::IntVector(size_t{m_key_max - m_key_min + 2}, std::max(size_t{1}, math::ilog2_ceil(num-1)), false);
    m_hi_idx[0] = 0; // left border of the first interval is the first entry
    assert(m_hi_idx[0] == 0);
    
    // std::cout << "keys[0] = " << keys[0] << ", hi=" << hi(keys[0]) << std::endl;
    // std::cout << "-> sample 0" << std::endl;
    
    size_t prev_key = m_key_min;
    for(size_t i = 1; i < num; i++) {
        const size_t cur_key = size_t{hi(keys[i])};
        // std::cout << "keys[" << i << "] = " << keys[i] << ", hi=" << cur_key << std::endl;
        if(cur_key > prev_key) {
            for(size_t key = prev_key + 1; key <= cur_key; key++) {
                m_hi_idx[key - size_t{m_key_min}] = i - 1;
                // std::cout << "-> sample " << i << std::endl;
                assert(m_hi_idx[0] == 0);
            }
        }
        prev_key = cur_key;
    }

    assert(prev_key == m_key_max);
    m_hi_idx[size_t{m_key_max - m_key_min + 1}] = num - 1;
    
    // std::cout << "sample: ";
    // for(size_t h = m_key_min; h <= m_key_max; h++) {
        // std::cout << m_hi_idx[h - m_key_min] << ' ';
    // }
    // std::cout << std::endl;
}

    Index(const Index& other) = default;
    Index(Index&& other) = default;
    Index& operator=(const Index& other) = default;
    Index& operator=(Index&& other) = default;

    /// \brief Finds the rank of the predecessor of the specified key.
    /// \param keys the keys that the compressed trie was constructed for
    /// \param num the number of keys
    /// \param x the key in question
    PosResult predecessor(const t_word* keys, const size_t num, const t_word x) const {
    if(tdc_unlikely(x < m_min))  return PosResult { false, 0 };
    if(tdc_unlikely(x >= m_max)) return PosResult { true, num - 1 };

    const uint64_t key = hi(x) - m_key_min;
    assert(key + 1 < m_hi_idx.size());
    const size_t q = m_hi_idx[key+1];

    if(x == keys[q]) {
        return PosResult { true, q };
    } else {
        const size_t p = m_hi_idx[key];
        return BinarySearchHybrid<t_word>::predecessor_seeded(keys, p, q, x);
    }
}
};

}} // namespace tdc::pred

