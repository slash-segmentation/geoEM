#pragma once

#include <vector>
#include <functional>
#include <utility>
#include <stdexcept>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>

namespace detail
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// A heap collection, which associates every value in the collection with a key and allows quick
// removal of the value with the most extreme key. The most extreme key could either be the one
// with the smallest value (min-heap, default) or the largest value (max-heap) or some custom
// definition. It also supports changing the values of keys by specifying a map type.
// 
// Template arguments:
//   K is the key type
//   V is the value type
//   Compare is a type to compare keys and defaults to std::less<K> (min-heap)
//   Mapping is a type to use to map values to keys and defaults to NoMapping which disables key updating
//   Container is the internal container used to store entries and defaults to std::vector
//
// If Mapping is not NoMapping it must take at least 2 class template arguments (for key and value)
// and must support the V type as a key and a pointer as a value. Additionally, it must provide:
//   types:     const_iterator
//   constructors: (default)
//   functions: erase(V)
//              clear()
//              std::pair<iterator,bool> insert(std::pair<V,void*>)
//              const_iterator find(V) const
//              const_iterator end() const
//              optional: reserve(size_type) [ or rehash(size_type) and double max_load_factor() ]
// std::map and std::unordered_map fill these requirements.
//
// Container must support a single template argument that can be a pointer type. It must provide
// the following:
//   types:     size_type
//   constructors: (default)
//   functions: V& operator[](size_type)
//              size_type size() const
//              reserve(size_type)
//              empty()
//              clear()
//              push_back()
//              pop_back()
//              back()
//
// Note: because C++ templates with default arguments do not map nicely to template-template
// classes you need to use template alias like the following to specify Mapping and Container:
// template<class T> using MyContainer = std::vector<T, MyAllocator>
// template<class K, class V> using MyMapping = std::unordered_map<K, V, MyHashFunction>
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, class V, class Compare, class Heap, class Entry, template<class> class Container>
class heap_base
{
    friend Entry;
    typedef heap_base<K,V,Compare,Heap,Entry,Container> Self;
public:
    typedef K key_type;
    typedef V value_type;
    typedef typename Container<Entry*>::size_type size_type;
protected:
    const Compare _comp;
    Container<Entry*> _heap;

    inline void _swap(size_type a, size_type b) { Entry *e = this->_heap[a]; (this->_heap[a] = this->_heap[b])->update(a); (this->_heap[b] = e)->update(b); }
    inline void _heapify_up(size_type pos)
    {
        while (pos > 0)
        {
            size_type parentPos = (pos - 1) / 2; // heap[i] has children heap[2*i + 1] and heap[2*i + 2] and parent heap[(i-1)/2];
            if (this->_comp(this->_heap[pos]->key, this->_heap[parentPos]->key)) { this->_swap(parentPos, pos); pos = parentPos; }
            else break;
        }
    }
    inline void _heapify_down(size_type pos)
    {
        const size_type n = this->_heap.size();
        const key_type& k = this->_heap[pos]->key;
        size_type smallest = pos;
        for(;;) // on each iteration exchange element with its smallest child
        {
            const size_type l = 2 * pos + 1, r = l + 1; // heap[i] has children heap[2*i + 1] and heap[2*i + 2] and parent heap[(i-1)/2];
            key_type key;
            if (l < n && this->_comp(this->_heap[l]->key, k)) { smallest = l; key = this->_heap[l]->key; } else { key = k; }
            if (r < n && this->_comp(this->_heap[r]->key, key)) { smallest = r; }
            else if (smallest == pos) { break; }
            this->_swap(smallest, pos);
            pos = smallest;
        }
    }
    template <class InputIterator>
    inline void _push_all(InputIterator first, InputIterator last)
    {
        this->_heap.reserve(this->_heap.size() + std::distance(first, last));
        for (size_type i = 0; first != last; ++first, ++i) { this->push_raw(first->first, first->second); }
        this->heapify();
    }
    inline void _clear_entries() { for (size_type i = 0, count = this->_heap.size(); i < count; ++i) { delete this->_heap[i]; } }

public:
    explicit inline heap_base(const Compare& comp = Compare()) : _comp(comp), _heap() { }
    inline virtual ~heap_base() { this->_clear_entries(); }
    
    // Basic container functions
    inline bool empty() const { return this->_heap.empty(); }
    inline size_type size() const { return this->_heap.size(); }
    inline virtual void reserve(size_type count) { this->_heap.reserve(count); }
    inline virtual void clear()
    {
        this->_clear_entries();
        this->_heap.clear();
    }

    // Add elements
    inline void push(const value_type& val, const key_type& key) { this->push_raw(val, key); this->_heapify_up(this->_heap.size() - 1); }
    inline void push_raw(const value_type& val, const key_type& key) { this->_heap.push_back(Entry::create(key, val, (Heap*)this)); }
    inline void heapify() { if (this->_heap.size() > 2) { for (size_type pos = (this->_heap.size() - 2) / 2; pos > 0; --pos) { _heapify_down(pos); } } if (this->_heap.size() > 0) { _heapify_down(0); } }

    // Get elements
    inline       value_type& top()       { return this->_heap[0]->val; }
    inline const value_type& top() const { return this->_heap[0]->val; }
    inline const key_type& top_key() const { return this->_heap[0]->key; }
    inline void pop() { Entry *e = this->_heap[0], *E = this->_heap.back(); this->_heap.pop_back(); e->cleanup((Heap*)this); if (e != E) { (this->_heap[0] = E)->update(0); this->_heapify_down(0); } delete e; }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template magic to determine if a Mapping has a reserve or rehash function.
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
#define TYPENAME 
#else
#define TYPENAME typename
#endif
#define CREATE_HAS_X(NAME,RT,...) \
    template<typename T> \
    class has_##NAME \
    { \
        typedef char yes[1]; typedef char no[2]; \
        template <typename V> static yes &chk(V* v, typename std::enable_if<std::is_same<decltype(v->NAME(__VA_ARGS__)),RT>::value>::type* = 0); \
        template <typename  > static no  &chk(...); \
    public: \
        const static bool value = sizeof(yes) == sizeof(chk<T>(0)); \
        typedef has_##NAME<T> type; \
    };
CREATE_HAS_X(reserve, void, TYPENAME V::size_type())
CREATE_HAS_X(rehash,  void, TYPENAME V::size_type())
namespace MPL = ::boost::mpl;
#define ENABLE_IF(...) typename ::std::enable_if< __VA_ARGS__::value >::type
template<class T> ENABLE_IF(has_reserve<T>)                                     map_reserve(T& m, typename T::size_type n) { m.reserve(n); } // has reserve
template<class T> ENABLE_IF(MPL::and_<has_rehash<T>,MPL::not_<has_reserve<T>>>) map_reserve(T& m, typename T::size_type n) { m.rehash(TYPENAME T::size_type(ceil(n/m.max_load_factor())+0.5)); } // has rehash but not reserve
template<class T> ENABLE_IF(MPL::not_<MPL::or_<has_reserve<T>,has_rehash<T>>>)  map_reserve(T&  , typename T::size_type  ) { } // has neither reserve or rehash
#undef ENABLE_IF
#undef CREATE_HAS_X
#undef TYPENAME

/////////////////////////////////////////////////
// A heap entry for a heap that has a Mapping
/////////////////////////////////////////////////
template<class H>
struct heap_with_map_entry
{
    typedef heap_with_map_entry<H> Entry;
    typedef typename H::key_type K;
    typedef typename H::value_type V;
    typedef typename H::size_type size_type;
    K key; V val; size_type index;
    inline heap_with_map_entry(const K& key, const V& val, size_type idx) : key(key), val(val), index(idx) { }
    inline static Entry* create(const K& key, const V& val, H* h)
    {
        Entry* e = new Entry(key, val, h->_heap.size());
        if (!h->_map.insert(std::make_pair(val,&e->index)).second) { delete e; throw std::out_of_range("value already exists"); }
        return e;
    }
    inline void update(size_type idx) { this->index = idx; }
    inline void cleanup(H* h) { h->_map.erase(this->val); }
};

/////////////////////////////////////////////////
// A heap entry for a heap without a Mapping
/////////////////////////////////////////////////
template<class H>
struct heap_without_map_entry
{
    typedef heap_without_map_entry<H> Entry;
    typedef typename H::key_type K;
    typedef typename H::value_type V;
    typedef typename H::size_type size_type;
    K key; V val;
    inline heap_without_map_entry(const K& key, const V& val) : key(key), val(val) { }
    inline static Entry* create(const K& key, const V& val, H* h) { return new Entry(key, val); }
    inline void update(size_type idx) { }
    inline void cleanup(H* h) { }
};

// A aliased vector with one template argument
template<class T> using vector = std::vector<T>;
} // end namespace detail

template <class, class> class NoMapping;

///////////////////////////////////////////////////////////////////////////////////////////////////
// A heap with a Mapping, adds the extra functions find, update, decrease, and increase.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, class V, class Compare = std::less<K>, template<class,class> class Mapping = NoMapping, template<class> class Container = detail::vector>
class heap :
    public detail::heap_base<K, V, Compare, heap<K,V,Compare,Mapping,Container>, detail::heap_with_map_entry<heap<K,V,Compare,Mapping,Container>>, Container>
{
    typedef heap<K,V,Compare,Mapping,Container> Self;
    typedef detail::heap_with_map_entry<Self> Entry;
    typedef detail::heap_base<K,V,Compare,Self,Entry,Container> Base;
    typedef Mapping<V, typename Base::size_type*> Map;
    friend Entry;
    Map _map;
public:
    typedef typename Base::size_type* handle;
    explicit inline heap(const Compare& comp = Compare()) : Base(comp), _map() { }
    template <class InputIterator>
    inline heap(InputIterator first, InputIterator last, const Compare& comp = Compare()) : Base(comp), _map() { this->_push_all(first, last); }
    inline virtual ~heap() { }
    inline virtual void reserve(typename Base::size_type count) { Base::reserve(count); detail::map_reserve(this->_map, count); }
    inline virtual void clear() { Base::clear(); this->_map.clear(); }

    // Finding values and updating their keys
    inline handle find(const typename Base::value_type& val) const { typename Map::const_iterator i = this->_map.find(val); return i == this->_map.end() ? nullptr : i->second; }
    inline bool update(handle h, const typename Base::key_type& key)
    {
        assert(h != nullptr);
        typename Base::size_type i = *h;
        const typename Base::key_type& old = this->_heap[i]->key;
        if (this->_comp(old, key))      { this->_heap[i]->key = key; this->_heapify_down(i); }
        else if (this->_comp(key, old)) { this->_heap[i]->key = key; this->_heapify_up(i); }
        else { return false; }
        return true;
    }
    inline bool decrease(handle h, const typename Base::key_type& key)
    {
        assert(h != nullptr);
        typename Base::size_type i = *h;
        if (!this->_comp(key, this->_heap[i]->key)) { return false; }
        this->_heap[i]->key = key;
        this->_heapify_up(i);
        return true;
    }
    inline bool increase(handle h, const typename Base::key_type& key)
    {
        assert(h != nullptr);
        typename Base::size_type i = *h;
        if (!this->_comp(this->_heap[i]->key, key)) { return false; }
        this->_heap[i]->key = key;
        this->_heapify_down(i);
        return true;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// A heap without a Mapping.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, class V, class Compare>
class heap<K, V, Compare, NoMapping, detail::vector> :
    public detail::heap_base<K, V, Compare, heap<K, V, Compare, NoMapping, detail::vector>, detail::heap_without_map_entry<heap<K, V, Compare, NoMapping, detail::vector>>, detail::vector>
{
    typedef heap<K,V,Compare,NoMapping,detail::vector> Self;
    typedef detail::heap_without_map_entry<Self> Entry;
    typedef detail::heap_base<K,V,Compare,Self,Entry,detail::vector> Base;
public:
    explicit inline heap(const Compare& comp = Compare()) : Base(comp) { }
    template <class InputIterator>
    inline heap(InputIterator first, InputIterator last, const Compare& comp = Compare()) : Base(comp) { this->_push_all(first, last); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// A heap without a Mapping but a custom Container.
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class K, class V, class Compare, template<class> class Container>
class heap<K, V, Compare, NoMapping, Container> :
    public detail::heap_base<K, V, Compare, heap<K, V, Compare, NoMapping, Container>, detail::heap_without_map_entry<heap<K, V, Compare, NoMapping, Container>>, Container>
{
    typedef heap<K,V,Compare,NoMapping,Container> Self;
    typedef detail::heap_without_map_entry<Self> Entry;
    typedef detail::heap_base<K, V, Compare, Self, Entry, Container> Base;
public:
    explicit inline heap(const Compare& comp = Compare()) : Base(comp) { }
    template <class InputIterator>
    inline heap(InputIterator first, InputIterator last, const Compare& comp = Compare()) : Base(comp) { this->_push_all(first, last); }
};
