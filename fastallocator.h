#include <iostream>
#include <type_traits>
#include <vector>

template <size_t ChunkSize>
class FixedAllocator {
private:
    struct Block {
        char data[ChunkSize];
        Block* next;
    };

    size_t sz = 1;

    Block* free_block = nullptr;
    std::vector<Block*> blocks;

    FixedAllocator() {
        if (sz < 512) {
            sz *= 2;
        }

        free_block = new Block[sz];

        blocks.push_back(free_block);

        for (size_t i = 0; i < sz - 1; ++i) {
            free_block[i].next = &free_block[i + 1];
        }

        free_block[sz - 1].next = nullptr;
    }

public:
    FixedAllocator(const FixedAllocator&) = delete;
    FixedAllocator& operator=(const FixedAllocator&) = delete;

    static FixedAllocator& get_instance () {
        static FixedAllocator instance;
        return instance;
    }

    void* allocate() {
        if (free_block == nullptr) {
            if (sz < 512) {
                sz *= 2;
            }

            free_block = new Block[sz];

            blocks.push_back(free_block);

            for (size_t i = 0; i < sz - 1; ++i) {
                free_block[i].next = &free_block[i + 1];
            }

            free_block[sz - 1].next = nullptr;
        }

        Block* res = free_block;
        free_block = free_block->next;

        return static_cast<void*>(res);
    }

    void deallocate(void* ptr) {
        Block* block = static_cast<Block*>(ptr);
        block->next = free_block;
        free_block = block;
    }

    ~FixedAllocator() {
        for (size_t i = 0; i < blocks.size(); ++i) {
            delete blocks[i];
        }
    }
};

template <typename T>
class FastAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;

    template <class U>
    struct rebind {
        using other = FastAllocator<U>;
    };

    struct propagate_on_container_copy_assignment: std::false_type {};

public:
    FastAllocator() = default;

    template <class U>
    FastAllocator<T>(const FastAllocator<U>&) {}

    template <class U>
    FastAllocator<T>& operator=(const FastAllocator<U>&) {
        return *this;
    }

    FastAllocator<T>& select_on_container_copy_construction() {
        return *this;
    }

    const FastAllocator<T>& select_on_container_copy_construction() const {
        return *this;
    }

    T* allocate(size_t n) {
        if (n == 1) {
            return reinterpret_cast<T*>(FixedAllocator<sizeof(T)>::get_instance().allocate());
        } else {
            return reinterpret_cast<T*>(::operator new(n * sizeof(T)));
        }
    }

    void deallocate(T* ptr, size_t n) {
        if (n == 1) {
            FixedAllocator<sizeof(T)>::get_instance().deallocate(reinterpret_cast<void*>(ptr));
        } else {
            ::operator delete(ptr);
        }
    }

    template <typename ...Args>
    void construct(T* ptr, const Args& ...args) {
        new (ptr) T(args...);
    }

    void destroy(T* ptr) {
        ptr->~T();
    }

    bool operator==(const FastAllocator<T>&) {
        return true;
    }

    bool operator!=(const FastAllocator<T>& another) {
        return !(*this == another);
    }

    ~FastAllocator() {}
};

template <typename T, typename Allocator = std::allocator<T>>
class List {
private:
    struct Node {
        T value;
        Node* next = nullptr;
        Node* previous = nullptr;

        Node() = default;
        Node(const T& value, Node* next, Node* previous): value(value), next(next), previous(previous) {}
    };

    size_t sz = 0;
    Node* head = nullptr;
    Node* fake = nullptr;

    using RebindAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;
    RebindAllocator alloc;
    using AllocTraits = typename std::allocator_traits<Allocator>::template rebind_traits<Node>;

    void swap(List<T, Allocator>& copy) {
        std::swap(sz, copy.sz);
        std::swap(head, copy.head);
        std::swap(fake, copy.fake);
    }

public:
    template<bool IsConst>
    struct common_iterator: public std::iterator<std::bidirectional_iterator_tag, T> {
    private:
        friend class List<T, Allocator>;

        using conditional_node_ptr = typename std::conditional<IsConst, const Node*, Node*>::type;
        using conditional_link = typename std::conditional<IsConst, const T&, T&>::type;
        using conditional_ptr = typename std::conditional<IsConst, const T*, T*>::type;

        conditional_node_ptr ptr = fake;

    public:
        explicit common_iterator(conditional_node_ptr ptr): ptr(const_cast<conditional_node_ptr>(ptr)) {}

        common_iterator(const common_iterator<false>& another): ptr(const_cast<conditional_node_ptr>(another.ptr)) {}

        common_iterator& operator++() {
            ptr = ptr->next;
            return *this;
        }

        common_iterator operator++(int) {
            common_iterator copy = *this;

            ptr = ptr->next;
            return copy;
        }

        common_iterator& operator--() {
            ptr = ptr->previous;
            return *this;
        }

        common_iterator operator--(int) {
            common_iterator copy = *this;

            ptr = ptr->previous;
            return copy;
        }

        bool operator==(const common_iterator& another) const {
            return ptr == another.ptr;
        }

        bool operator!=(const common_iterator& another) const {
            return !(*this == another);
        }

        conditional_link operator*() const {
            return (ptr->value);
        }

        conditional_ptr operator->() const {
            return &(ptr->value);
        }
    };

    using iterator = common_iterator<false>;
    using const_iterator = common_iterator<true>;

    iterator begin() {
        return iterator(head);
    }

    const_iterator begin() const {
        return const_iterator(head);
    }

    const_iterator cbegin() const {
        return const_iterator(head);
    }

    iterator end() {
        return iterator(fake);
    }

    const_iterator end() const {
        return const_iterator(fake);
    }

    const_iterator cend() const {
        return const_iterator(fake);
    }

    template<bool IsConst>
    struct common_reverse_iterator: public std::iterator<std::bidirectional_iterator_tag, T> {
    private:
        friend class List<T, Allocator>;

        using conditional_node_ptr = typename std::conditional<IsConst, const Node*, Node*>::type;
        using conditional_link = typename std::conditional<IsConst, const T&, T&>::type;
        using conditional_ptr = typename std::conditional<IsConst, const T*, T*>::type;
        using conditional_iterator = typename std::conditional<IsConst, const_iterator, iterator>::type;

        conditional_node_ptr ptr = fake;

    public:
        explicit common_reverse_iterator(Node* ptr): ptr(const_cast<conditional_node_ptr>(ptr)) {}

        common_reverse_iterator(const common_reverse_iterator<false>& another): ptr(const_cast<conditional_node_ptr>(another.ptr)) {}

        common_reverse_iterator& operator--() {
            ptr = ptr->next;
            return *this;
        }

        common_reverse_iterator operator--(int) {
            common_reverse_iterator copy = *this;

            ptr = ptr->next;
            return copy;
        }

        common_reverse_iterator& operator++() {
            ptr = ptr->previous;
            return *this;
        }

        common_reverse_iterator operator++(int) {
            common_reverse_iterator copy = *this;

            ptr = ptr->previous;
            return copy;
        }

        iterator base() {
            return iterator(const_cast<Node*>(ptr)->next);
        }

        const_iterator base() const {
            return const_iterator(const_cast<Node*>(ptr)->next);
        }

        bool operator==(const common_reverse_iterator& another) const {
            return ptr == another.ptr;
        }

        bool operator!=(const common_reverse_iterator& another) const {
            return !(*this == another);
        }

        conditional_link operator*() const {
            return (ptr->value);
        }

        conditional_ptr operator->() const {
            return &(ptr->value);
        }
    };

    using reverse_iterator = common_reverse_iterator<false>;
    using const_reverse_iterator = common_reverse_iterator<true>;

    reverse_iterator rbegin() {
        return reverse_iterator(fake->previous);
    }

    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(fake->previous);
    }

    reverse_iterator rend() {
        return reverse_iterator(head->previous);
    }

    const_reverse_iterator rend() const {
        return const_reverse_iterator(head->previous);
    }

    const_reverse_iterator crbegin() const {
        return const_reverse_iterator(fake->previous);
    }

    const_reverse_iterator crend() const {
        return const_reverse_iterator(head->previous);
    }

    explicit List(const Allocator& alloc = Allocator()) {
        this->alloc = alloc;

        fake = AllocTraits::allocate(this->alloc, 1);
        fake->next = fake;
        fake->previous = fake;

        head = fake;
    }

    List(size_t count, const Allocator& alloc = Allocator()) {
        this->alloc = alloc;

        fake = AllocTraits::allocate(this->alloc, 1);
        fake->next = fake;
        fake->previous = fake;

        head = fake;

        Node* prev = fake;

        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(this->alloc, 1);
            AllocTraits::construct(this->alloc, new_node);
            new_node->next = fake;
            new_node->previous = prev;
            prev = new_node;

            fake->previous->next = new_node;
            fake->previous = new_node;

            if (head == fake) {
                head = new_node;
                fake->next = new_node;
            }
        }

        sz = count;
    }

    List(size_t count, const T& value, const Allocator& alloc = Allocator()) {
        this->alloc = alloc;

        fake = AllocTraits::allocate(this->alloc, 1);
        fake->next = fake;
        fake->previous = fake;

        head = fake;

        Node* prev = fake;

        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(this->alloc, 1);
            AllocTraits::construct(this->alloc, new_node, value, fake, prev);
            prev = new_node;

            fake->previous->next = new_node;
            fake->previous = new_node;

            if (head == fake) {
                head = new_node;
                fake->next = new_node;
            }
        }

        sz = count;
    }

    List<T, Allocator>(const List<T, Allocator>& another) {
        alloc = AllocTraits::select_on_container_copy_construction(another.alloc);

        fake = AllocTraits::allocate(this->alloc, 1);
        fake->next = fake;
        fake->previous = fake;

        head = fake;

        Node* ptr = another.head;
        Node* prev_node = head;
        while (ptr != another.fake){
            Node* new_node = AllocTraits::allocate(alloc, 1);
            AllocTraits::construct(alloc, new_node, (*ptr).value, fake, prev_node);
            ptr = ptr->next;

            fake->previous->next = new_node;
            fake->previous = new_node;

            if (head == fake) {
                head = new_node;
                fake->next = new_node;
            }

            prev_node = new_node;
        }

        sz = another.sz;
    }

    List<T, Allocator>& operator=(const List<T, Allocator>& another) {
        List<T, Allocator> copy = another;
        swap(copy);

        if (AllocTraits::propagate_on_container_copy_assignment::value && alloc != another.alloc) {
            alloc = another.alloc;
        }

        return *this;
    }

    size_t size() const {
        return sz;
    }

    RebindAllocator get_allocator() const {
        return alloc;
    }

    void push_front(const T& value) {
        Node* new_node = AllocTraits::allocate(alloc, 1);
        AllocTraits::construct(alloc, new_node, value, head, fake);

        fake->next = new_node;
        head->previous = new_node;
        head = new_node;

        ++sz;
    }

    void push_back(const T& value) {
        Node* new_node = AllocTraits::allocate(alloc, 1);
        AllocTraits::construct(alloc, new_node, value, fake, fake->previous);

        fake->previous->next = new_node;
        fake->previous = new_node;

        if (head == fake) {
            head = new_node;
            fake->next = new_node;
        }

        ++sz;
    }

    void pop_front() {
        head->next->previous = fake;
        Node* for_pop = head;
        head = head->next;

        AllocTraits::destroy(alloc, for_pop);
        AllocTraits::deallocate(alloc, for_pop, 1);
        --sz;
    }

    void pop_back() {
        Node* for_pop = fake->previous;
        fake->previous = for_pop->previous;
        for_pop->previous->next = fake;

        if (head == for_pop) {
            fake->next = fake;
            fake->previous = fake;
            head = fake;
        }

        AllocTraits::destroy(alloc, for_pop);
        AllocTraits::deallocate(alloc, for_pop, 1);
        --sz;
    }

    void insert(iterator it, const T& value) {
        Node* new_node = AllocTraits::allocate(alloc, 1);
        AllocTraits::construct(alloc, new_node, value, const_cast<Node*>(it.ptr), it.ptr->previous);
        it.ptr->previous->next = new_node;
        it.ptr->previous = new_node;

        if (head == it.ptr) {
            head = new_node;
        }

        ++sz;
    }

    void insert(const_iterator it, const T& value) {
        Node* new_node = AllocTraits::allocate(alloc, 1);
        AllocTraits::construct(alloc, new_node, value, const_cast<Node*>(it.ptr), it.ptr->previous);
        it.ptr->previous->next = new_node;
        const_cast<Node*>(it.ptr)->previous = new_node;

        if (head == it.ptr) {
            head = new_node;
        }

        ++sz;
    }

    void erase(iterator it) {
        it.ptr->next->previous = it.ptr->previous;
        it.ptr->previous->next = it.ptr->next;

        if (head == it.ptr) {
            head = it.ptr->next;
        }

        AllocTraits::destroy(alloc, const_cast<Node*>(it.ptr));
        AllocTraits::deallocate(alloc, const_cast<Node*>(it.ptr), 1);

        --sz;
    }

    void erase(const_iterator it) {
        it.ptr->next->previous = it.ptr->previous;
        it.ptr->previous->next = it.ptr->next;

        if (head == it.ptr) {
            head = it.ptr->next;
        }

        AllocTraits::destroy(alloc, const_cast<Node*>(it.ptr));
        AllocTraits::deallocate(alloc, const_cast<Node*>(it.ptr), 1);

        --sz;
    }

    ~List() {
        Node* ptr = head;

        while (ptr != fake) {
            Node* next = ptr->next;
            AllocTraits::destroy(alloc, ptr);
            AllocTraits::deallocate(alloc, ptr, 1);
            ptr = next;
        }

        AllocTraits::destroy(alloc, fake);
        AllocTraits::deallocate(alloc, fake, 1);
    }

    void test() {
        Node* ptr = head;

        while (ptr != fake) {
            std::cout << ptr->value << ' ';
            ptr = ptr->next;
        }

        std::cout << '\n';
    }
};
