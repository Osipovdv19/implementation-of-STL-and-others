#include <bits/stdc++.h>

template <size_t N>
class StackStorage {
private:
    char pool_[N];
    void *ptr_ = pool_;
public:
    StackStorage(): ptr_(pool_) {}

    void* getPool(size_t sz, size_t al) {
        size_t tmp = N;
        void *value = std::align(al, sz, ptr_, tmp);
        ptr_ = reinterpret_cast<char*>(ptr_) + sz;
        return value;
    }

    ~StackStorage() = default;
};

template <typename T, size_t N>
class StackAllocator {
private:
    StackStorage<N>* stg_;
public:
    typedef T value_type;

    StackAllocator() = default;

    StackAllocator(StackStorage<N>& s): stg_(&s) {}

    T* allocate(size_t sz) {
        return reinterpret_cast<T*>(stg_->getPool(sz * sizeof(T), alignof(T)));
    }

    template <typename V>
    struct rebind {
        typedef StackAllocator<V, N> other;
    };

    void deallocate(T*, size_t) const {};

    StackStorage<N>* getStg() const {
        return stg_;
    }

    template <typename V>
    StackAllocator(const StackAllocator<V, N>& al): stg_(al.getStg()) {}

    template <typename V>
    StackAllocator& operator=(const StackAllocator<V, N>& al) {
        stg_ = al.getStg();
        return *this;
    }

    bool operator==(const StackAllocator& v) {
        return stg_ == v.getStg();
    }

    bool operator!=(const StackAllocator& v) {
        return stg_ != v.getStg();
    }

    ~StackAllocator() {};
};

template <typename T, typename Allocator = std::allocator<T>>
class List {
private:
    struct BaseNode {
        BaseNode *next = nullptr;
        BaseNode *prev = nullptr;
    };

    struct Node: BaseNode {
        T value;

        Node(const T& v): value(v) {}

        template<typename... Args>
        Node(Args&&... args) : value(std::forward<Args>(args)...) {}

        Node(const T&& v): value(std::move(const_cast<T&>(v))) {}

        template<typename... Args>
        Node(const Args&... args) : value(args...) {}

        Node(T&& v) : value(std::move(const_cast<T&>(v))) {}
    };

    typedef std::allocator_traits<Allocator> Traits;
    typedef typename Traits::template rebind_alloc<Node> AllocatorNode;
    typedef typename Traits::template rebind_traits<Node> TraitsNode;

    Allocator al_;
    AllocatorNode alNode_ = al_;
    BaseNode* base_;
    size_t sz_ = 0;

public:
    static void matchNodes(BaseNode* v, BaseNode* u) {
        v->next = u;
        u->prev = v;
    }

    template <bool Flag>
    struct Iterator {
        typedef std::conditional_t<Flag, const T*, T*> pointer;
        typedef std::conditional_t<Flag, const T&, T&> reference;
        typedef std::conditional_t<Flag, const T, T> value_type;
        typedef std::ptrdiff_t difference_type;
        typedef std::bidirectional_iterator_tag iterator_category;

        BaseNode *ptr;

        Iterator(BaseNode* v): ptr(v) {}

        Iterator(const Iterator<false>& it) : ptr(it.ptr) {}

        template <bool V>
        bool operator!=(const Iterator<V>& it) const {
            return ptr != it.ptr;
        }

        template <bool V>
        bool operator==(const Iterator<V>& it) const {
            return ptr == it.ptr;
        }

        pointer operator->() const {
            return &(static_cast<std::conditional_t<Flag, const Node*, Node*>>(ptr)->value);
        }

        Iterator& operator--() {
            ptr = ptr->prev;
            return *this;
        }

        Iterator& operator++() {
            ptr = ptr->next;
            return *this;
        }

        const Iterator operator++(int) {
            Iterator temp = *this;
            ptr = ptr->next;
            return temp;
        }

        const Iterator operator--(int) {
            Iterator temp = *this;
            ptr = ptr->prev;
            return temp;
        }

        reference operator*() {
            return reinterpret_cast<Node*>(ptr)->value;
        }
    };

    typedef Iterator<true> const_iterator;
    typedef Iterator<false> iterator;
    typedef std::reverse_iterator<Iterator<true>> const_reverse_iterator;
    typedef std::reverse_iterator<Iterator<false>> reverse_iterator;

    iterator begin() const {
        return iterator(base_->next);
    }

    const_iterator cbegin() const {
        return const_iterator(base_->next);
    }

    reverse_iterator rbegin() const {
        return reverse_iterator(base_);
    }

    const_reverse_iterator crbegin() const {
        return const_reverse_iterator(base_);
    }

    iterator end() const {
        return iterator(base_);
    }

    const_iterator cend() const {
        return const_iterator(base_);
    }

    reverse_iterator rend() const {
        return reverse_iterator(base_->next);
    }

    const_reverse_iterator crend() const {
        return const_reverse_iterator(base_->next);
    }

    Allocator get_allocator() const {
        return al_;
    }

    size_t size() const {
        return sz_;
    }

    template <typename U>
    void remove(U* v) {
        TraitsNode::destroy(alNode_, reinterpret_cast<Node*>(v));;
        TraitsNode::deallocate(alNode_, reinterpret_cast<Node*>(v), 1);;
    }

    void push_back(const T& val) {
        Node* node = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, node, val);
        sz_++;
        BaseNode* temp = base_->prev;
        matchNodes(temp, node);
        matchNodes(node, base_);
    }

    void pop_back() {
        BaseNode* v = base_->prev->prev;
        BaseNode* u = base_->prev;
        matchNodes(v, base_);
        sz_--;
        remove(u);
    }

    void push_front(const T& val) {
        Node* node = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, node, val);
        sz_++;
        BaseNode* temp = base_->next;
        matchNodes(base_, node);
        matchNodes(node, temp);
    }

    void pop_front() {
        BaseNode* v = base_->next->next;
        BaseNode* u = base_->next;
        matchNodes(base_, v);
        sz_--;
        remove(u);
    }

    List& operator=(const List<T, Allocator>& v) {
        if (std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value)
            al_ = v.get_allocator();
        size_t tempSz = sz_;
        alNode_ = al_;
        BaseNode* tempBase = v.base_;
        for (size_t i = 0; i < v.size(); ++i) {
            tempBase = tempBase->next;
            try {
                push_back(reinterpret_cast<Node*>(tempBase)->value);
            }
            catch (...) {
                for (size_t j = 0; j < i; ++j) {
                    pop_back();
                }
            }
        }
        for (size_t i = 0; i < tempSz; ++i)
            pop_front();
        return *this;
    }

    void refreshBase() {
        base_ = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, base_);
        base_->next = base_->prev = base_;
    }

    List& operator=(List&& v) {
        base_ = std::move(v.base_);
        alNode_ = std::move(v.alNode_);
        al_ = std::move(v.al_);
        sz_ = v.sz_;
        v.sz_ = 0;
        v.base_ = nullptr;
        return *this;
    }

    List(const Allocator& v = Allocator()): al_(v) {
        refreshBase();
    }

    void constructor(BaseNode* tempV = NULL, size_t siz = 0, bool fl = false) {
        refreshBase();
        BaseNode *temp = base_;
        for (size_t i = 0; i < siz; i++) {
            Node *node;
            if (fl)
                tempV = tempV->next;
            try {
                node = TraitsNode::allocate(alNode_, 1);
                if (fl)
                    TraitsNode::construct(alNode_, node, reinterpret_cast<Node *>(tempV)->value);
                else
                    TraitsNode::construct(alNode_, node);
            }
            catch (...) {
                temp = base_->next;
                for (size_t j = 0; j < static_cast<size_t>(i) - 2; j++) {
                    temp = temp->next;
                    remove(temp->prev);
                }
                remove(node);
                remove(temp);
                throw;
            }
            matchNodes(temp, node);
            temp = node;
            if (!fl)
                base_->prev = node;
        }
        matchNodes(temp, base_);
    }
    List(const List<T, Allocator>& v) {
        alNode_ = al_ = std::allocator_traits<Allocator>::select_on_container_copy_construction(v.get_allocator());
        constructor(v.base_, v.size(), true);
        sz_ = v.size();
    }

    List(size_t siz, const Allocator& v = Allocator()): al_(v) {
        constructor(NULL, siz, false);
        sz_ = siz;
    }

    List(size_t siz, const T& u, const Allocator& v = Allocator()) {
        refreshBase();
        BaseNode* temp = base_;
        for (size_t i = 0; i < siz; ++i) {
            Node* node = TraitsNode::allocate(alNode_, 1);
            TraitsNode::construct(alNode_, node, u);
            matchNodes(temp, node);
            temp = node;
        }
        sz_ = siz;
        matchNodes(temp, base_);
    }

    template <bool Flag>
    void erase(Iterator<Flag> iter) {
        BaseNode* v = iter.ptr->prev;
        BaseNode* u = iter.ptr->next;
        matchNodes(v, u);
        sz_--;
        remove(iter.ptr);
    }

    void clear() {
        while (sz_) {
            erase(this->begin());
        }
    }

    template <bool Flag>
    void insert(Iterator<Flag> iter, T&& V) {
        Node* node = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, node, std::move(const_cast<T&>(V)));
        BaseNode* v = iter.ptr->prev;
        BaseNode* u = iter.ptr;
        sz_++;
        matchNodes(v, node);
        matchNodes(node, u);
    }

    template<bool Flag, typename... Args>
    void emplace(const Iterator<Flag>& it, Args&&... args) {
        Node* v = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, v, std::forward<Args>(args)...);
        BaseNode* tmp = reinterpret_cast<BaseNode*>(v);
        matchNodes(it.ptr->prev, tmp);
        matchNodes(tmp, it.ptr);
        sz_++;
    }

    ~List() {
        clear();
        TraitsNode::destroy(alNode_, base_);
        TraitsNode::deallocate(alNode_, reinterpret_cast<Node*>(base_), 1);
    }
};

template<
        typename Key,
        typename Value,
        typename Hash = std::hash<Key>,
        typename Equal = std::equal_to<Key>,
        typename Allocator = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap {
public:
    typedef std::pair<const Key, Value> NodeType;

private:
    static size_t getHash(const Key& v) {
        return Hash{}(v);
    }

    static size_t checkEqual(const Key& v, const Key& u) {
        return Equal{}(v, u);
    }

    struct Node {
        size_t hash_;
        NodeType key_;

        Node(const NodeType& v) : key_(v) {
            hash_ = getHash(key_.first);
        }

        Node(NodeType&& v) : key_(
                std::move(const_cast<Key&>(v.first)),
                std::move(const_cast<Value&>(v.second))) {
            hash_ = getHash(key_.first);
        }

        Node(const Key& key, const Value& val) : key_(key, val) {
            hash_ = getHash(key_.first);
        }

        Node(Key&& key, Value&& val) : key_(std::move(key), std::move(val)) {
            hash_ = getHash(key_.first);
        }

        Node(Node&& v) : hash_(v.hash_), key_(
                std::move(const_cast<Key&>(v.key_.first)),
                std::move(v.key_.second)) {
        }
    };

    typedef std::allocator_traits<Allocator> Traits;
    typedef typename Traits::template rebind_alloc<Node> AllocatorNode;
    typedef std::allocator_traits<AllocatorNode> TraitsNode;
    typedef typename List<Node, AllocatorNode>::iterator BucketIterator;

    double mxLoadFct_ = 0.911;
    size_t sz_ = 0;
    AllocatorNode alNode_ = Allocator();
    List<Node, AllocatorNode> lst_;
    std::vector<BucketIterator> buckets_;

public:
    template <bool Flag>
    struct Iterator {
        typedef std::conditional_t<Flag, const NodeType*, NodeType*> pointer;
        typedef std::conditional_t<Flag, const NodeType&, NodeType&> reference;
        typedef std::conditional_t<Flag, const NodeType, NodeType> value_type;
        typedef std::ptrdiff_t difference_type;
        typedef std::forward_iterator_tag iterator_category;

        BucketIterator it_;

        Iterator(const BucketIterator& v) : it_(v) {}

        reference operator*() {
            return (it_->key_);
        }

        pointer operator->() {
            return &(it_->key_);
        }

        Iterator& operator++() {
            ++it_;
            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp = *this;
            ++it_;
            return tmp;
        }

        template<bool FlagOther>
        bool operator==(const Iterator<FlagOther>& v) const {
            return it_ == v.it_;
        }

        template<bool FlagOther>
        bool operator!=(const Iterator<FlagOther>& v) const {
            return !(*this == v);
        }
    };

    typedef Iterator<false> iterator;
    typedef Iterator<true> const_iterator;

    iterator begin() const {
        return iterator(lst_.begin());
    }

    const_iterator cbegin() const {
        return const_iterator(lst_.begin());
    }

    iterator end() const {
        return iterator(lst_.end());
    }

    const_iterator cend() const {
        return const_iterator(lst_.end());
    }

    size_t size() const {
        return lst_.size();
    }

    UnorderedMap() {
        buckets_.assign(1, lst_.end());
    }

    UnorderedMap(const UnorderedMap& v){
        buckets_.assign(1, lst_.end());
        mxLoadFct_ = v.mxLoadFct_;
        for (auto& node : v.lst_)
            insert(node.key_);
    }

    UnorderedMap(UnorderedMap&& v) {
        lst_ = std::move(v.lst_);
        sz_ = v.sz_;
        mxLoadFct_ = v.mxLoadFct_;
        buckets_ = std::move(v.buckets_);
        for (auto& key : buckets_)
            if (v.lst_.end() == key) key = lst_.end();
    }

    UnorderedMap& operator=(UnorderedMap&& v) {
        lst_ = std::move(v.lst_);
        sz_ = v.sz_;
        mxLoadFct_ = v.mxLoadFct_;
        v.sz_ = 0;
        buckets_ = std::move(v.buckets_);
        for (auto& key : buckets_)
            if (v.lst_.end() == key) key = lst_.end();
        return *this;
    }

    void remove(Node* v) {
        TraitsNode::destroy(alNode_, v);
        TraitsNode::deallocate(alNode_, v, 1);
    }

    template<typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        if (mxLoadFct_ < load_factor())
            reserve(buckets_.size() + 1);
        Node* node = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, node, std::forward<Args>(args)...);
        NodeType key(
                std::move(const_cast<Key&>((node->key_).first)),
                std::move((node->key_).second));
        remove(node);
        size_t id = (getHash(key.first)) % buckets_.size();
        auto it = find(key.first);
        if (it != end())
            return {it, false};
        lst_.emplace(
                buckets_[id],
                std::move(const_cast<Key&>(key.first)),
                std::move(key.second));
        if (buckets_[id] == lst_.end())
            sz_++;
        return {iterator(--buckets_[id]), true};
    }

    std::pair<iterator, bool> insert(NodeType&& v) {
        return emplace(std::move(const_cast<Key&>(v.first)), std::move(v.second));
    }

    std::pair<iterator, bool> insert(const NodeType& v) {
        return emplace(std::move(NodeType(v)));
    }

    template<class InputIterator>
    void insert(InputIterator v, InputIterator u) {
        for (auto iter = v; iter != u; ++iter) insert(*iter);
    }

    void erase(iterator it) {
        size_t id = (it.it_->hash_) % buckets_.size();
        if (iterator(buckets_[id]) != it) {
            lst_.erase(it.it_);
        }
        else {
            lst_.erase((it++).it_);
            iterator tmp = it;
            if (tmp == end() || tmp.it_->hash_ % buckets_.size() != id) {
                buckets_[id] = lst_.end();
                sz_--;
            }
            else
                buckets_[id] = tmp.it_;
        }
    }

    void erase(iterator v, iterator u) {
        for (auto it = v; it != u; ) erase(it++);
    }

    auto findHelp(const Key& key) {
        size_t id = (getHash(key)) % buckets_.size();
        for (auto keyb = buckets_[id]; keyb != lst_.end(); ++keyb) {
            if (checkEqual(keyb->key_.first, key))
                return iterator(keyb);
            else if (keyb->hash_ % buckets_.size() != id) {
                break;
            }
        }
        return end();
    }

    iterator find(const Key& key) {
        return findHelp(key);
    }

    const_iterator find(const Key& key) const {
        return findHelp(key);
    }

    Value& operator[](const Key& key) {
        auto ans = find(key);
        if (ans == end()) ans = emplace(key, Value()).first;
        return ans->second;
    }

    Value& at(const Key& key) {
        auto ans = find(key);
        if (ans == end()) throw "Error with method - at";
        return ans->second;
    }

    size_t max_size() const {
        return 1e9;
    }

    double load_factor() const {
        return double(sz_) / buckets_.size();
    }

    double max_load_factor() const {
        return mxLoadFct_;
    }

    void max_load_factor(double v) {
        mxLoadFct_ = v;
    }

    void rehash(size_t v) {
        buckets_.resize(v, lst_.end());
        for (auto &bucket : buckets_) bucket = lst_.end();
        sz_ = 0;
        size_t id;
        for (auto key = lst_.begin(); key != lst_.end(); ) {
            id = key->hash_ % buckets_.size();
            auto tmp = key++;
            auto i1 = buckets_[id];
            auto i2 = tmp;
            lst_.matchNodes(i2.ptr->prev, i2.ptr->next);
            auto node1 = i1.ptr->next;
            auto node2 = i2.ptr;
            lst_.matchNodes(node1->prev, node2);
            lst_.matchNodes(node2, node1);
            if (buckets_[id] == lst_.end())
                sz_++;
            buckets_[id] = tmp;
        }
    }

    void reserve(size_t v) {
        rehash(std::ceil(double(v) / mxLoadFct_));
    }

    ~UnorderedMap() {}
};
