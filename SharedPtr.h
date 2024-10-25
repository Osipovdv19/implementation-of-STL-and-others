#include <memory>
#include <iostream>

struct BaseControlBlock {
    size_t weakCnt = 0;
    size_t sharedCnt = 0;

    virtual void remove() = 0;
    virtual void deallocate() = 0;

    virtual void* getPtr() = 0;
    virtual ~BaseControlBlock() = default;
};

template<typename Type, typename Deleter, typename Allocator>
struct ControlBlock: BaseControlBlock {
public:
    Type* ptr;
    Deleter deleter;

    typedef std::allocator_traits<Allocator> Traits;
    typedef typename Traits::template rebind_alloc<ControlBlock> BlockAllocator;

    BlockAllocator alloc;

    ControlBlock(Type* type, Deleter deleter, Allocator allocator)
            : ptr(type), deleter(deleter), alloc(allocator) {}

    void remove() override {
        deleter(ptr); //
        ptr = nullptr;
    }

    void deallocate() override {
        alloc.deallocate(this, 1);
    }

    void* getPtr() override { return ptr; }

    ~ControlBlock() {
        deallocate();
        remove();
    }
};

template<typename Type, typename Allocator>
struct SharedControlBlock: BaseControlBlock {
    Type object;

    typedef std::allocator_traits<Allocator> Traits;
    typedef typename Traits::template rebind_alloc<SharedControlBlock> BlockAllocator;

    BlockAllocator alloc;

    void remove() override {
        alloc.destroy(&object);
    }

    template<typename... Args>
    SharedControlBlock(Allocator allocator, Args&&... args)
            : object(std::forward<Args>(args)...), alloc(allocator) {}

    void deallocate() override {
        alloc.deallocate(this, 1);
    }

    void* getPtr() override { return &object; }

    ~SharedControlBlock() {
        deallocate();
        remove();
    }
};

template<class Type>
class SharedPtr {
private:
    template<typename V, typename Alloc, typename... Args>
    friend SharedPtr<V> allocateShared(Alloc alloc, Args&&... args);

    template<class V>
    friend class SharedPtr;

    template<class V>
    friend class WeakPtr;

    BaseControlBlock* block_ = nullptr;

    template<typename Allocator, typename... Args>
    SharedPtr(Allocator all, Args&&... args) {
        typedef std::allocator_traits<Allocator> Traits;
        typedef typename Traits::template rebind_alloc<SharedControlBlock<Type, Allocator>> BlockAllocator;
        BlockAllocator bl = all;
        auto ptr = std::allocator_traits<BlockAllocator>::allocate(bl, 1);
        bl.construct(ptr, bl, std::forward<Args>(args)...);
        ptr->sharedCnt = 1;
        block_ = ptr;
    }

    SharedPtr(BaseControlBlock* block) : block_(block) {
        if (block) {
            block->sharedCnt++;
        }
    }

public:

    template <typename V, typename Deleter = std::default_delete<V>, typename Allocator = std::allocator<V>> ///
    SharedPtr(V* ptr, Deleter del, Allocator all) {
        typedef std::allocator_traits<Allocator> Traits;
        typedef typename Traits::template rebind_alloc<ControlBlock<V, Deleter, Allocator>> BlockAllocator;
        typename Traits::template rebind_alloc<ControlBlock<V, Deleter, Allocator>> bl = all;
        block_ = std::allocator_traits<BlockAllocator>::allocate(bl, 1); ///
        new (block_) ControlBlock<V, Deleter, Allocator>(ptr, del, std::move(bl));
        block_->sharedCnt = 1;
    }

    SharedPtr() = default;

    template <typename V>
    explicit SharedPtr(V* ptr): SharedPtr(ptr, std::default_delete<V>(), std::allocator<V>()) {} ///

    template <typename V, typename Deleter>
    SharedPtr(V* ptr, Deleter del): SharedPtr(ptr, del, std::allocator<V>()) {}

    SharedPtr(SharedPtr<Type>&& u) : block_(u.block_) {
        u.block_ = nullptr;
    }

    SharedPtr(const SharedPtr<Type>& u) : SharedPtr(u.block_) {}

    template<typename V>
    SharedPtr(SharedPtr<V>&& u) : block_(u.block_) {
        u.block_ = nullptr;
    }

    template<typename V>
    SharedPtr(const SharedPtr<V>& u) : SharedPtr(u.block_) {}

    template<typename V>
    void swap(SharedPtr<V>& u) {
        std::swap(u.block_, block_);
    }

    Type* operator->() const {
        if (!block_)
            return nullptr;
        auto t = block_->getPtr();
        return reinterpret_cast<Type*>(t);
    }

    Type& operator*() const {
        auto t = block_->getPtr();
        return *reinterpret_cast<Type*>(t);
    }

    SharedPtr& operator=(SharedPtr<Type>&& u) {
        SharedPtr v = std::move(u);
        swap(v);
        return *this;
    }

    SharedPtr& operator=(const SharedPtr<Type>& u) {
        SharedPtr v = u;
        swap(v);
        return *this;
    }

    template<typename V>
    SharedPtr& operator=(SharedPtr<V>&& u) {
        SharedPtr v = std::move(u);
        swap(v);
        return *this;
    }

    template<typename V>
    SharedPtr& operator=(const SharedPtr<V>& u) {
        SharedPtr v = u;
        swap(v);
        return *this;
    }

    size_t use_count() const {
        return block_->sharedCnt;
    }

    template<typename V>
    void reset(V* ptr) {
        SharedPtr<Type> v(ptr);
        swap(v);
    }

    void reset() {
        *this = SharedPtr<Type>();
    }

    Type* get() const { ///
        return block_ ? reinterpret_cast<Type*>(block_->getPtr()) : nullptr;
    }

    ~SharedPtr() {
        if (!block_) return;
        --block_->sharedCnt;
        if (!block_->sharedCnt) block_->remove();
        else return;
        if (!block_->weakCnt) block_->deallocate();
    }
};

template<typename T, typename Allocator, typename... Args>
SharedPtr<T> allocateShared(Allocator allocator, Args&&... args) {
    return SharedPtr<T>(allocator, std::forward<Args>(args)...);
}

template<typename T, typename... Args>
SharedPtr<T> makeShared(Args&&... args) {
    return allocateShared<T, std::allocator<T>, Args...>(std::allocator<T>(), std::forward<Args>(args)...);
}

template<class T>
class WeakPtr {
private:
    template<class V>
    friend class WeakPtr;

    template<typename V>
    friend class SharedPtr;

    BaseControlBlock* block_ = nullptr;

    WeakPtr(BaseControlBlock* u) : block_(u) {
        if (block_) {
            block_->weakCnt++;
        }
    }

public:

    WeakPtr(WeakPtr<T>&& u) : block_(u.block_) {
        u.block_ = nullptr;
    }

    WeakPtr(const WeakPtr<T>& u) : WeakPtr(u.block_) {}

    template<typename V>
    WeakPtr(WeakPtr<V>&& u) : block_(u.block_) {
        u.block_ = nullptr;
    }

    template<typename V>
    WeakPtr(const WeakPtr<V>& u) : WeakPtr(u.block_) {}

    WeakPtr() = default;

    template<typename V>
    void swap(WeakPtr<V>& u) {
        std::swap(u.block_, block_);
    }

    template<typename V>
    WeakPtr(const SharedPtr<V>& u) : WeakPtr(u.block_) {}

    WeakPtr(const SharedPtr<T>& u) : WeakPtr(u.block_) {}

    size_t use_count() const {
        return block_->sharedCnt;
    }

    bool expired() const {
        return !block_->sharedCnt;
    }

    WeakPtr& operator=(WeakPtr<T>&& u) {
        WeakPtr v = std::move(u);
        swap(v);
        return *this;
    }

    WeakPtr& operator=(const WeakPtr<T>& u) {
        WeakPtr v = u;
        swap(v);
        return *this;
    }

    template<typename V>
    WeakPtr& operator=(WeakPtr<V>&& u) {
        WeakPtr v = std::move(u);
        swap(v);
        return *this;
    }

    template<typename V>
    WeakPtr& operator=(const WeakPtr<V>& u) {
        WeakPtr v = u;
        swap(v);
        return *this;
    }

    SharedPtr<T> lock() const {
        return SharedPtr<T>(!expired() ? block_ : nullptr);
    }

    ~WeakPtr() {
        if (!block_) return;
        --block_->weakCnt;
        if (!block_->weakCnt && !block_->sharedCnt) block_->deallocate();
    }
};