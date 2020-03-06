#pragma once
    
#include <cstring>
#include <memory>
#include <type_traits>
#include <utility>

#include "item_ref.hpp"
#include <tdc/math/imath.hpp>

namespace tdc {
namespace vec {

/// \brief A static vector of items, i.e., items cannot be inserted or deleted.
///
/// This has the tiny advantage over \c std::vector that no capacity has to be maintained.
///
/// Note that the accessors will always handle \em copies of items rather than references or rvalues.
/// This vector should therefore only be used for non-complex types that have no background allocations or similar.
template<typename T>
class StaticVector {
private:
	static_assert(!std::is_same<T, bool>::value, "A StaticVector of boolean values is not supported. You'll want to use a BitVector instead.");

	friend class ItemRef<StaticVector<T>, T>;

	static constexpr size_t s_item_size = sizeof(T);
	
    static std::unique_ptr<T[]> allocate(const size_t num, const bool initialize = true) {
		T* p = new T[num];
		
		if(initialize) {
			memset(p, T(), num * s_item_size);
		}
		
		return std::unique_ptr<T[]>(p);
	}

    size_t m_size;
    std::unique_ptr<T[]> m_data;

    uint64_t get(const size_t i) const {
		return m_data[i];
	}
	
    void set(const size_t i, const T v) {
		m_data[i] = v;
	}
	
public:
	/// \brief Proxy for reading and writing a single integer.
	using ItemRef_ = ItemRef<StaticVector<T>, T>;

    /// \brief Constructs an empty vector.
    inline StaticVector() : m_size(0) {
    }
	
    /// \brief Copy constructor.
    /// \param other the vector to copy
    inline StaticVector(const StaticVector& other) {
        *this = other;
    }

    /// \brief Move constructor.
    /// \param other the vector to move
    inline StaticVector(StaticVector&& other) {
        *this = std::move(other);
    }

    /// \brief Constructs a vector with the specified length.
    /// \param size the number of items
	/// \param initialize if \c true, the items will be initialized using their default constructor
    inline StaticVector(const size_t size, const bool initialize = true) {
        m_size = size;
        m_data = allocate(size, initialize);
    }

    /// \brief Copy assignment.
    /// \param other the vector to copy
    inline StaticVector& operator=(const StaticVector& other) {
        m_size = other.m_size;
        m_data = allocate(m_size, false);
        memcpy(m_data.get(), other.m_data.get(), m_size * s_item_size);
        return *this;
    }

    /// \brief Move assignment.
    /// \param other the vector to move
    inline StaticVector& operator=(StaticVector&& other) {
        m_size = other.m_size;
        m_data = std::move(other.m_data);
        return *this;
    }
	
    /// \brief Resizes the vector with the specified new length.
    ///
    /// \param size the new number of items
    inline void resize(const size_t size) {
		StaticVector<T> new_v(size, size >= m_size); // no init if new size is smaller
		for(size_t i = 0; i < m_size; i++) {
			new_v.set(i, get(i));
		}
		
		m_size = new_v.m_size;
		m_data = std::move(new_v.m_data);
    }

    /// \brief Reads the specified item.
    /// \param i the number of the item to read
    inline T operator[](const size_t i) const {
        return get(i);
    }

    /// \brief Accesses the specified item.
    /// \param i the number of the item to access
    inline ItemRef_ operator[](const size_t i) {
        return ItemRef_(*this, i);
    }

    /// \brief The number of items in the vector.
    inline size_t size() const {
        return m_size;
    }
};

}} //namespace tdc::vec