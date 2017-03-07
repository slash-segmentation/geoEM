#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// stable_vector_iterator "wraps" a vector iterator making it not invalidate during the operation
// push_back(...), insert(end(), ...), reserve(...), and shrink_to_size(). Not "stable" (does not
// necessarily stay pointing at the same value) during operations like remove and insert.
// This does make it ~20% slower to do either * or -> (due to an additional dereference) and twice
// as large (size of 2 pointers instead of just 1).
// V is the vector type and I is the iterator type (either V::iterator or V::const_interator)
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename V, typename I>
class stable_vector_iterator
{
	template <typename V2, typename I2> friend class stable_vector_iterator;
	V* v;
	typename V::size_type i;
public:
	typedef stable_vector_iterator<V,I> type;
	typedef typename I::value_type value_type;
	typedef typename I::pointer pointer;
	typedef typename I::reference reference;
	typedef typename V::size_type size_type;
	typedef typename I::difference_type difference_type;
	typedef typename I::iterator_category iterator_category;

	inline stable_vector_iterator() : v(nullptr), i(0) { }
	inline stable_vector_iterator(V* v, size_type i) : v(v), i(i) { }
	inline stable_vector_iterator(const type& rhs) : v(rhs.v), i(rhs.i) { }
	template <typename V2, typename I2> inline stable_vector_iterator(const stable_vector_iterator<V2,I2>& rhs) : v(rhs.v), i(rhs.i) { }
	inline type& operator=(const type& rhs) { this->v = rhs.v; this->i = rhs.i; return *this; }

	inline bool operator==(const type& rhs) const { return this->v == rhs.v && this->i == rhs.i; }
	inline bool operator!=(const type& rhs) const { return this->v != rhs.v || this->i != rhs.i; }
	inline bool operator< (const type& rhs) const { return this->i <  rhs.i; }
	inline bool operator<=(const type& rhs) const { return this->i <= rhs.i; }
	inline bool operator> (const type& rhs) const { return this->i >  rhs.i; }
	inline bool operator>=(const type& rhs) const { return this->i >= rhs.i; }
	inline difference_type operator-(const type& rhs) const { return this->i - rhs.i; }
	inline reference operator*() { return this->v->operator[](i); }
	inline pointer operator->() { return &this->v->operator[](i); }
	inline reference operator*() const { return this->v->operator[](i); }
	inline pointer operator->() const { return &this->v->operator[](i); }
	inline reference operator[](difference_type n) const { return this->v->operator[](i+n); }
	inline type& operator++() { ++this->i; return *this; }
	inline type& operator--() { --this->i; return *this; }
	inline type operator++(int) { type x = *this; ++this->i; return x; }
	inline type operator--(int) { type x = *this; --this->i; return x; }
	inline type& operator+=(difference_type n) { this->i += n; return *this; }
	inline type& operator-=(difference_type n) { this->i -= n; return *this; }
	inline type operator+(difference_type n) const { return type(this->v, this->i + n); }
	inline type operator-(difference_type n) const { return type(this->v, this->i - n); }
	friend inline type operator+(difference_type n, const type& rhs) { return type(rhs.v, rhs.i + n); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// handle_iterator is an iterator that dereferences every value one additional time.
// I is an iterator of handles/iterators
// H is the handle type to return (typically I::value_type for non-const iterators)
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename I, typename H>
class handle_iterator
{
	template <typename I2, typename H2> friend class handle_iterator;
	I i;
public:
	typedef handle_iterator<I, H> type;
	typedef I iter_type;
	typedef H handle_type;
	typedef typename handle_type::value_type value_type;
	typedef handle_type pointer;
	typedef value_type& reference;
	typedef size_t      size_type;
	typedef ptrdiff_t   difference_type;
	typedef typename I::iterator_category iterator_category;
	inline handle_iterator() : i() { }
	inline handle_iterator(iter_type i) : i(i) { }
	template <class I2, class H2> inline handle_iterator(const handle_iterator<I2, H2>& i) : i(i.i) { }
	inline bool operator==(const type& rhs) const { return this->i == rhs.i; }
	inline bool operator!=(const type& rhs) const { return this->i != rhs.i; }
	inline bool operator==(const iter_type& rhs) const { return this->i == rhs; }
	inline bool operator!=(const iter_type& rhs) const { return this->i != rhs; }
	inline bool operator==(const handle_type& rhs) const { return *this->i == rhs; }
	inline bool operator!=(const handle_type& rhs) const { return *this->i != rhs; }
	//friend inline bool operator==(const iter_type& lhs, const type& rhs) { return lhs == rhs->i } // these are covered by casts
	//friend inline bool operator!=(const iter_type& lhs, const type& rhs) { return lhs != rhs->i; }
	//friend inline bool operator==(const handle_type& lhs, const type& rhs) { return lhs == *rhs->i }
	//friend inline bool operator!=(const handle_type& lhs, const type& rhs) { return lhs != *rhs->i; }
	inline bool operator< (const type& rhs) const { return this->i <  rhs.i; }
	inline bool operator<=(const type& rhs) const { return this->i <= rhs.i; }
	inline bool operator> (const type& rhs) const { return this->i >  rhs.i; }
	inline bool operator>=(const type& rhs) const { return this->i >= rhs.i; }
	inline ptrdiff_t operator-(const type& rhs) const { return this->i - rhs.i; }
	inline operator reference() const { return **this->i; }
	inline operator pointer() const { return *this->i; }
	inline reference operator*() const { return **this->i; }
	inline pointer operator->() const { return *this->i; }
	inline reference operator[](ptrdiff_t n) const { return *this->i[n]; }
	inline type& operator++() { ++this->i; return *this; }
	inline type& operator--() { --this->i; return *this; }
	inline type operator++(int) { type x = *this; ++this->i; return x; }
	inline type operator--(int) { type x = *this; --this->i; return x; }
	inline type& operator+=(ptrdiff_t n) { this->i += n; return *this; }
	inline type& operator-=(ptrdiff_t n) { this->i -= n; return *this; }
	inline type operator+(ptrdiff_t n) const { return type(this->i + n); }
	inline type operator-(ptrdiff_t n) const { return type(this->i - n); }
	friend inline type operator+(ptrdiff_t n, const type& rhs) { return type(rhs.i + n); }
};
