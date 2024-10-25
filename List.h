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
        void *value = std::align(al, sz, ptr_, tmp); ///
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

        Node() = default;

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
    BaseNode base_tmp;
    BaseNode* base_ = &base_tmp;
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
    void insert(Iterator<Flag> iter, const T& V) {
        Node* node = TraitsNode::allocate(alNode_, 1);
        TraitsNode::construct(alNode_, node, V);
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