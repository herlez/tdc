#pragma once

#include <cstddef>
#include <cstdint>

#include <tdc/pred/result.hpp>
#include <tdc/util/skip_accessor.hpp>

namespace tdc {
namespace pred {
namespace dynamic {

using Result = ::tdc::pred::Result;

class DynamicFusionNode8 {
private:
    size_t m_size;
    uint64_t m_key[8], m_index;
    uint64_t m_mask, m_branch, m_free;
    uint8_t m_bkey;

    size_t find(const uint64_t key) const;
    size_t rank(const uint64_t key) const;
    bool used(const size_t j) const;

public:
    /// \brief Constructs an empty fusion node.
    DynamicFusionNode8();

    /// \brief Selects the key with the given rank.
    /// \param i the rank in question
    uint64_t select(const size_t i) const;

    /// \brief Convenience operator for \ref select.
    /// \param i the rank in question    
    inline uint64_t operator[](const size_t i) const {
        return select(i);
    }
    
    /// \brief Finds the rank of the predecessor of the specified key in the node.
    /// \param x the key in question
    Result predecessor(const uint64_t x) const;

    /// \brief Inserts the specified key.
    /// \param key the key to insert
    void insert(const uint64_t key);

    /// \brief Removes the specified key.
    /// \param key the key to remove
    void remove(const uint64_t key);

    /// \brief Returns the current size of the fusion node.
    inline size_t size() const {
        return m_size;
    }

    // FIXME DEBUG
    void print() const;
};

}}} // namespace tdc::pred::dynamic
