#pragma once

#include <cassert>
#include <memory>
#include <utility>

#include "bit_vector.hpp"
#include "int_vector.hpp"
#include "rank_u64.hpp"
#include "select_u64.hpp"

#include <tdc/math/imath.hpp>

namespace tdc {
namespace vec {


/// \brief A space efficient data structure for answering select queries on a \ref BitVector in constant time.
///
/// A select query, given a number \em k, finds the position of the k-th occurence of a set or unset bit, respectively, in the bit vector, starting from the beginning.
/// The data structure uses a hierarchical scheme dividing the bit vector into \em blocks and \em superblocks
/// and precomputes a number of key positions, storing them in a space efficient manner.
/// On the lowest level, the search is accelerated using the \c ctz instruction.
///
/// Note that this data structure is \em static.
/// It maintains a pointer to the underlying bit vector and will become invalid if that bit vector is changed after construction.
///
/// \tparam m_bit the bit we are interested in: 1 for finding set bits, 0 for finding unset bits
template<bool m_bit>
class BitSelect {
private:
    static constexpr uint8_t basic_rank(uint64_t v);
    static constexpr uint8_t basic_rank(uint64_t v, uint8_t x);
    static constexpr uint8_t basic_rank(uint64_t v, uint8_t a, uint8_t b);
    static constexpr uint8_t basic_select(uint64_t v, uint8_t k);
    static constexpr uint8_t basic_select(uint64_t v, uint8_t l, uint8_t k);

    std::shared_ptr<const BitVector> m_bv;

    size_t m_max;
    size_t m_block_size;
    size_t m_supblock_size;
    size_t m_blocks_per_supblock;

    IntVector m_blocks;
    IntVector m_supblocks;

public:
    /// \brief Constructs the rank data structure for the given bit vector.
    /// \param bv the bit vector
    BitSelect(std::shared_ptr<const BitVector> bv) : m_bv(bv) {
        const size_t n = m_bv->size();
        const size_t log_n = math::ilog2_ceil(n - 1);

        // construct
        m_block_size    = log_n;
        m_supblock_size = log_n * log_n;
        m_blocks_per_supblock = log_n;

        m_supblocks = IntVector(math::idiv_ceil(n, m_supblock_size), log_n);
        m_blocks = IntVector(math::idiv_ceil(n, m_block_size), log_n);

        m_max = 0;
        size_t r_sb = 0; // current bit count in superblock
        size_t r_b = 0;  // current bit count in block

        size_t cur_sb = 0;        // current superblock
        size_t cur_sb_offset = 0; // starting position of current superblock
        size_t longest_sb = 0;    // length of longest superblock

        size_t cur_b = 0; // current block

        const size_t num_blocks = m_bv->num_blocks();
        
        for(size_t i = 0; i < num_blocks; i++) {
            const auto v = m_bv->block64(i);
            
            uint8_t r;
            if(i + 1 < num_blocks) {
                // get amount flag bits in whole data element
                r = basic_rank(v);
            } else {
                // only get amount of flag bits up to end of bit vector
                const size_t m = n & 63ULL; // mod 64
                if(m > 0) {
                    r = basic_rank(v, m-1);
                } else {
                    r = basic_rank(v); // full element, but last one
                }
            }

            m_max += r;

            if(r_b + r >= m_block_size) {
                // entered new block

                // amount of bits needed to fill current block
                size_t distance_b = m_block_size - r_b;

                // stores the offset of the last bit in the current block
                uint8_t offs = 0;

                r_b += r;

                size_t distance_sum = 0;
                while(r_b >= m_block_size) {
                    // find exact position of the bit in question
                    offs = basic_select(v, offs, distance_b);
                    assert(SELECT_U64_FAIL != offs);

                    const size_t pos = (i << 6ULL) + offs; // mul 64

                    distance_sum += distance_b;
                    r_sb += distance_b;
                    if(r_sb >= m_supblock_size) {
                        // entered new superblock
                        longest_sb = std::max(longest_sb, pos - cur_sb_offset);
                        cur_sb_offset = pos;

                        m_supblocks[cur_sb++] = pos;

                        r_sb -= m_supblock_size;
                    }

                    m_blocks[cur_b++] = pos - cur_sb_offset;
                    r_b -= m_block_size;
                    distance_b = m_block_size;

                    ++offs;
                }

                assert(size_t(r) >= distance_sum);
                r_sb += r - distance_sum;
            } else {
                r_b  += r;
                r_sb += r;
            }
        }

        longest_sb = std::max(longest_sb, n - cur_sb_offset);
        const size_t w_block = math::ilog2_ceil(longest_sb);

        m_supblocks.resize(cur_sb, log_n);
        m_blocks.resize(cur_b, w_block);
    }

    /// \brief Constructs an empty, uninitialized select data structure.
    inline BitSelect()
        : m_bv(nullptr),
          m_max(0),
          m_block_size(0),
          m_supblock_size(0),
          m_blocks_per_supblock(0) {
    }

    /// \brief Copy constructor.
    /// \param other the object to copy
    inline BitSelect(const BitSelect& other) {
        *this = other;
    }

    /// \brief Move constructor.
    /// \param other the object to move
    inline BitSelect(BitSelect&& other) {
        *this = std::move(other);
    }

    /// \brief Copy assignment.
    /// \param other the object to copy
    inline BitSelect& operator=(const BitSelect& other) {
        m_bv = other.m_bv;
        m_max = other.m_max;
        m_block_size = other.m_block_size;
        m_supblock_size = other.m_supblock_size;
        m_blocks_per_supblock = other.m_blocks_per_supblock;
        m_blocks = other.m_blocks;
        m_supblocks = other.m_supblocks;
        return *this;
    }

    /// \brief Move assignment.
    /// \param other the object to move
    inline BitSelect& operator=(BitSelect&& other) {
        m_bv = std::move(other.m_bv);
        m_max = other.m_max;
        m_block_size = other.m_block_size;
        m_supblock_size = other.m_supblock_size;
        m_blocks_per_supblock = other.m_blocks_per_supblock;
        m_blocks = std::move(other.m_blocks);
        m_supblocks = std::move(other.m_supblocks);
        return *this;
    }

    /// \brief Finds the x-th occurence of \c m_bit in the bit vetor.
    /// \param x the rank of the occurence to find, must be greater than zero
    /// \return the position of the x-th occurence, or the size of the bit vector to indicate that there are no x occurences of \c m_bit
    inline size_t select(size_t x) const {
        assert(x > 0);
        if(x > m_max) return m_bv->size();

        size_t pos = 0;

        //narrow down to block
        {
            const size_t i = x / m_supblock_size;
            const size_t j = x / m_block_size;

            if(i > 0) {
                pos += m_supblocks[i-1];
                x -= i * m_supblock_size;
            }
            if(x == 0) return pos;

            // block j is the k-th block within the i-th superblock
            size_t k = j - i * m_blocks_per_supblock;
            if(k > 0) {
                pos += m_blocks[j-1];
                x   -= k * m_block_size;
            }
            if(x == 0) return pos;

            if(i > 0 || k > 0) ++pos; // offset from block boundary
        }

        // from this point forward, search directly in the bit vector
        size_t i = pos >> 6ULL; // div 64
        size_t offs  = pos & 63ULL; // mod 64

        uint64_t block = m_bv->block64(i);
        uint8_t s = basic_select(block, offs, x);
        if(s != SELECT_U64_FAIL) {
            // found in first data segment
            return pos + s - offs;
        } else {
            // linearly search in the next segments
            x -= basic_rank(block, offs, 63ULL);
            while(s == SELECT_U64_FAIL) {
                ++i;
                pos = i << 6ULL; // mul 64
                block = m_bv->block64(i);
                s = basic_select(block, x);
                if(s == SELECT_U64_FAIL) x -= basic_rank(block);
            };

            return pos + s;
        }
    }

    /// \brief Finds the x-th occurence of \c m_bit in the bit vetor.
    ///
    /// This is a convenience alias for \ref select.
    ///
    /// \param x the rank of the occurence to find, must be greater than zero
    /// \return the position of the x-th occurence, or the size of the bit vector to indicate that there are no x occurences of \c m_bit
    inline size_t operator()(size_t x) const {
        return select(x);
    }
};

/// \cond INTERNAL
template<>
inline constexpr uint8_t BitSelect<0>::basic_rank(uint64_t v) {
    return 64ULL - rank1_u64(v);
}

template<>
inline constexpr uint8_t BitSelect<0>::basic_rank(uint64_t v, uint8_t x) {
    return x + 1 - rank1_u64(v, x);
}

template<>
inline constexpr uint8_t BitSelect<0>::basic_rank(uint64_t v, uint8_t a, uint8_t b) {
    return (b-a+1) - rank1_u64(v, a, b);
}

template<>
inline constexpr uint8_t BitSelect<0>::basic_select(uint64_t v, uint8_t k) {
    return select0_u64(v, k);
}

template<>
inline constexpr uint8_t BitSelect<0>::basic_select(uint64_t v, uint8_t l, uint8_t k) {
    return select0_u64(v, l, k);
}

template<>
inline constexpr uint8_t BitSelect<1>::basic_rank(uint64_t v) {
    return rank1_u64(v);
}

template<>
inline constexpr uint8_t BitSelect<1>::basic_rank(uint64_t v, uint8_t x) {
    return rank1_u64(v, x);
}

template<>
inline constexpr uint8_t BitSelect<1>::basic_rank(uint64_t v, uint8_t a, uint8_t b) {
    return rank1_u64(v, a, b);
}

template<>
inline constexpr uint8_t BitSelect<1>::basic_select(uint64_t v, uint8_t k) {
    return select1_u64(v, k);
}

template<>
inline constexpr uint8_t BitSelect<1>::basic_select(uint64_t v, uint8_t l, uint8_t k) {
    return select1_u64(v, l, k);
}
/// \endcond

/// \brief Convenience type definition for \ref BitSelect for 0-bits.
using BitSelect0 = BitSelect<0>;

/// \brief Convenience type definition for \ref BitSelect for 1-bits.
using BitSelect1 = BitSelect<1>;

}} // namespace tdc::vec