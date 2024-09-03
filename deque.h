#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

template<typename T>
class Deque {
private:
    static const int chunk_size = 200;

    int sz = 0, capacity = 0;
    int shift = 0;

    std::pair<int, int> begin_pos, end_pos;

    T** data;

    void swap(Deque<T>& copy) {
        std::swap(sz, copy.sz);
        std::swap(capacity, copy.capacity);
        std::swap(shift, copy.shift);
        std::swap(begin_pos, copy.begin_pos);
        std::swap(end_pos, copy.end_pos);
        std::swap(data, copy.data);
    }

public:
    template<bool IsConst>
    struct common_iterator: public std::iterator<std::random_access_iterator_tag, typename std::conditional<IsConst, const T, T>::type> {
    private:
        using conditional_chunks_ptr = typename std::conditional<IsConst, const T**, T**>::type;
        using conditional_ptr = typename std::conditional<IsConst, const T*, T*>::type;
        using conditional_link = typename std::conditional<IsConst, const T&, T&>::type;
        using conditional_pos = typename std::conditional<IsConst, const int, int>::type;

        conditional_chunks_ptr chunks_ptr;
        conditional_pos chunk_pos;
        conditional_pos pos;
    public:
        common_iterator(T** chunks_ptr, conditional_pos chunk_pos, conditional_pos pos): chunks_ptr(const_cast<conditional_chunks_ptr>(chunks_ptr)), chunk_pos(chunk_pos), pos(pos) {}
        explicit common_iterator(common_iterator<false>& another): chunks_ptr(const_cast<conditional_chunks_ptr>(another.chunks_ptr)), chunk_pos(another.chunk_pos), pos(another.pos) {}

        common_iterator& operator++() {
            ++pos;
            if (pos == chunk_size) {
                pos = 0;
                ++chunk_pos;
            }

            return *this;
        }

        common_iterator& operator++(int) {
            common_iterator copy = *this;

            ++pos;
            if (pos == chunk_size) {
                pos = 0;
                ++chunk_pos;
            }

            return copy;
        }

        common_iterator& operator--() {
            --pos;
            if (pos < 0) {
                pos = chunk_size - 1;
                --chunk_pos;
            }

            return *this;
        }

        common_iterator& operator--(int) {
            common_iterator copy = *this;

            --pos;
            if (pos < 0) {
                pos = chunk_size - 1;
                --chunk_pos;
            }

            return copy;
        }

        common_iterator& operator+=(const int n) {
            if (n >= 0) {
                chunk_pos += n / chunk_size;
                pos += n % chunk_size;

                if (pos >= chunk_size) {
                    pos = pos % chunk_size;
                    ++chunk_pos;
                }
            } else {
                chunk_pos -= (-n) / chunk_size;
                pos -= (-n) % chunk_size;

                if (pos < 0) {
                    pos = chunk_size + pos;
                    --chunk_pos;
                }
            }

            return *this;
        }

        common_iterator operator+(const int n) const {
            common_iterator copy = *this;
            copy += n;
            return copy;
        }

        common_iterator& operator-=(const int n) {
            if (n >= 0) {
                chunk_pos -= n / chunk_size;
                pos -= n % chunk_size;

                if (pos < 0) {
                    pos = chunk_size + pos;
                    --chunk_pos;
                }
            } else {
                chunk_pos += (-n) / chunk_size;
                pos += (-n) % chunk_size;

                if (pos >= chunk_size) {
                    pos = pos % chunk_size;
                    ++chunk_pos;
                }
            }

            return *this;
        }

        common_iterator& operator-(const int n) const {
            common_iterator<false> copy (const_cast<T**>(chunks_ptr), chunk_pos, pos);
            copy -= n;
            return copy;
        }

        int operator-(const common_iterator& another) const {
            return (chunk_pos - another.chunk_pos - 1) * chunk_size + pos + chunk_size - another.pos;
        }

        bool operator==(const common_iterator& another) const {
            return ((chunks_ptr == another.chunks_ptr) && (chunk_pos == another.chunk_pos) && (pos == another.pos)) ||
                   ((chunk_pos == another.chunk_pos - 1) && (pos == chunk_size) && (another.pos == 0)) || 
                   ((chunk_pos - 1 == another.chunk_pos) && (pos == 0) && (another.pos == chunk_size));
        }

        bool operator!=(const common_iterator& another) const {
            return !(*this == another);
        }

        bool operator>(const common_iterator& another) const {
            return (*this - another) > 0;
        }

        bool operator<(const common_iterator& another) const {
            return !(*this > another);
        }

        bool operator>=(const common_iterator& another) const {
            return (*this > another) || (*this == another);
        }

        bool operator<=(const common_iterator& another) const {
            return (*this < another) || (*this == another);
        }

        conditional_link operator*() const {
            return chunks_ptr[chunk_pos][pos];
        }

        conditional_ptr operator->() const {
            return &chunks_ptr[chunk_pos][pos];
        }
    };

    using iterator = common_iterator<false>;
    using const_iterator = common_iterator<true>;

    iterator begin() {
        return iterator(data, begin_pos.first + shift, begin_pos.second);
    }

    const_iterator begin() const {
        return const_iterator(begin());
    }

    iterator end() {
        return iterator(data, end_pos.first + shift, end_pos.second);
    }

    const_iterator end() const {
        return const_iterator(end());
    }

    const_iterator cbegin() const {
        return const_iterator(begin());
    }

    const_iterator cend() const {
        return const_iterator(end());
    }

    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    reverse_iterator rbegin() const {
        return reverse_iterator(iterator(data, end_pos.first + shift, end_pos.second));
    }

    reverse_iterator rend() const {
        return reverse_iterator(iterator(data, begin_pos.first + shift, begin_pos.second));
    }

    const_reverse_iterator crbegin() const {
        return const_reverse_iterator(const_iterator(data, end_pos.first + shift, end_pos.second - 1));
    }

    const_reverse_iterator crend() const {
        return const_reverse_iterator(const_iterator(data, begin_pos.first + shift, begin_pos.second - 1));
    }

    Deque<T>() = default;

    Deque<T>(const Deque<T>& another) {
        T** old_data = data;

        data = new T*[another.capacity];

        int i = 0;
        try {
            for (; i < another.capacity; ++i) {
                data[i] = reinterpret_cast<T*>(new int8_t[chunk_size * sizeof (T)]);
                std::uninitialized_copy(another.data[i], another.data[i] + chunk_size, data[i]);
            }
        }  catch (...) {
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k < chunk_size; ++k) {
                    data[j][k].~T();
                }

                delete[] reinterpret_cast<int8_t*>(data[j]);
            }
            delete[] data;
            data = old_data;

            throw;
        }

        for (int i = 0; i < capacity; ++i) {
            for (int j = 0; j < chunk_size; ++j) {
                data[i][j].~T();
            }

            delete[] reinterpret_cast<int8_t*>(old_data[i]);
        }

        capacity = another.capacity;
        sz = another.sz;
        begin_pos = another.begin_pos;
        end_pos = another.end_pos;
        shift = another.shift;
    }

    Deque<T>(int n) {
        capacity = 2 * (n / chunk_size + (n % chunk_size != 0));
        sz = n;

        data = new T*[capacity];

        int i = 0;
        try {
            for (; i < capacity; ++i) {
                data[i] = reinterpret_cast<T*>(new int8_t[chunk_size * sizeof (T)]);
            }

            begin_pos = end_pos = std::make_pair(capacity / 2 - n / chunk_size / 2 - (n % (chunk_size * 2) != 0), chunk_size - n / 2 % chunk_size - 1);

            for (int i = 0; i < n; ++i) {
                std::uninitialized_default_construct_n(&data[end_pos.first + shift][end_pos.second], 1);
                if (end_pos.second == chunk_size - 1) {
                    ++end_pos.first;
                }
                end_pos.second = (end_pos.second + 1) % chunk_size;
            }
        } catch (...) {
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k < chunk_size; ++k) {
                    data[j][k].~T();
                }

                delete[] reinterpret_cast<int8_t*>(data[j]);
            }
            delete[] data;
            capacity = 0;
            sz = 0;
            begin_pos = end_pos = std::make_pair(0, 0);

            throw;
        }
    }

    Deque<T>(int n, const T& value) {
        capacity = 2 * (n / chunk_size + (n % chunk_size != 0));
        sz = n;

        data = new T*[capacity];

        int i = 0;
        try {
            for (; i < capacity; ++i) {
                data[i] = reinterpret_cast<T*>(new int8_t[chunk_size * sizeof (T)]);
            }

            begin_pos = end_pos = std::make_pair(capacity / 2 - n / chunk_size / 2 - (n % (chunk_size * 2) != 0), chunk_size - n / 2 % chunk_size - 1);

            for (int i = 0; i < n; ++i) {
                memcpy(&data[end_pos.first + shift][end_pos.second], &value, sizeof(T));
                if (end_pos.second == chunk_size - 1) {
                    ++end_pos.first;
                }
                end_pos.second = (end_pos.second + 1) % chunk_size;
            }
        } catch (...) {
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k < chunk_size; ++k) {
                    data[j][k].~T();
                }

                delete[] reinterpret_cast<int8_t*>(data[j]);
            }
            delete[] data;
            capacity = 0;
            sz = 0;
            begin_pos = end_pos = std::make_pair(0, 0);

            throw;
        }
    }

    Deque<T>& operator=(const Deque<T>& another) {
        Deque<T> copy = another;
        swap(copy);
        return *this;
    }

    T& operator[](int index) {
        return data[begin_pos.first + shift + (begin_pos.second + index) / chunk_size][(begin_pos.second + index) % chunk_size];
    }

    const T& operator[](int index) const {
        return data[begin_pos.first + shift + (begin_pos.second + index) / chunk_size][(begin_pos.second + index) % chunk_size];
    }

    T& at(int index) {
        if (index < 0 || index >= sz) {
            throw std::out_of_range("");
        }
        return data[begin_pos.first + shift + (begin_pos.second + index) / chunk_size][(begin_pos.second + index) % chunk_size];
    }

    const T& at(int index) const {
        if (index < 0 || index >= sz) {
            throw std::out_of_range("");
        }
        return data[begin_pos.first + shift + (begin_pos.second + index) / chunk_size][(begin_pos.second + index) % chunk_size];
    }

    size_t size() const {
        return sz;
    }

    void push_back(const T& value) {
        std::pair<int, int> prev_end_pos = end_pos;
        if (end_pos.second >= chunk_size) {
            ++end_pos.first;
            end_pos.second = 0;
        }

        if (end_pos.first + shift >= capacity) {
            T** old_data = data;

            data = new T*[2 * std::max(capacity, 1)];

            int i = 0;
            try {
                for (; i < capacity; ++i) {
                    data[i] = old_data[i];
                }
                for (; i < 2 * std::max(capacity, 1); ++i) {
                    data[i] = reinterpret_cast<T*>(new int8_t[chunk_size * sizeof (T)]);
                }
            }  catch (...) {
                for (int j = capacity; j < i; ++j) {
                    delete[] reinterpret_cast<int8_t*>(data[j]);
                }

                data = old_data;
                end_pos = prev_end_pos;

                throw;
            }

            delete[] old_data;

            capacity = 2 * std::max(capacity, 1);
        }

        try {
            memcpy(&data[end_pos.first + shift][end_pos.second], &value, sizeof(T));
            ++sz;
            ++end_pos.second;
        }  catch (...) {
            end_pos = prev_end_pos;

            throw;
        }
    }

    void pop_back() {
        --end_pos.second;
        (data[end_pos.first + shift] + end_pos.second)->~T();
        --sz;
        if (end_pos.second == 0) {
            --end_pos.first;
            end_pos.second = chunk_size;
        }

        if (end_pos.first + shift < capacity / 4) {
            T** old_data = data;

            data = new T*[capacity / 2];

            for (int i = 0; i < capacity / 2; ++i) {
                data[i] = old_data[i];
            }

            for (int i = capacity / 2; i < capacity; ++i) {
                delete[] reinterpret_cast<int8_t*>(old_data[i]);
            }
            delete[] old_data;

            capacity /= 2;
        }
    }

    void push_front(const T& value) {
        std::pair<int, int> prev_begin_pos = begin_pos;
        --begin_pos.second;

        if (begin_pos.second < 0) {
            --begin_pos.first;
            begin_pos.second = chunk_size - 1;
        }

        if (capacity == 0) {
            end_pos = std::make_pair(-1, 200);
        }

        if (begin_pos.first + shift < 0) {
            T** old_data = data;

            data = new T*[2 * std::max(capacity, 1)];

            int i = 2 * capacity - 1;

            try {
                for (; i >= capacity; --i) {
                    data[i] = old_data[i - capacity];
                }
                if (capacity == 0) {
                    i = 1;
                }
                for (; i >= 0; --i) {
                    data[i] = reinterpret_cast<T*>(new int8_t[chunk_size * sizeof (T)]);
                }
            }  catch (...) {
                if (capacity != 0) {
                    for (int j = capacity - 1; j > i; --j) {
                        delete[] reinterpret_cast<int8_t*>(data[j]);
                    }
                } else {
                    for (int j = 1; j > i; --j) {
                        delete[] reinterpret_cast<int8_t*>(data[j]);
                    }
                }
                data = old_data;
                begin_pos = prev_begin_pos;

                throw;
            }

            delete[] old_data;

            if (capacity == 0) {
                shift = 2;
            } else {
                shift = shift + capacity;
            }

            capacity = 2 * std::max(capacity, 1);
        }

        try {
            memcpy(&data[begin_pos.first + shift][begin_pos.second], &value, sizeof(T));
            ++sz;
        }  catch (...) {
            begin_pos = prev_begin_pos;

            throw;
        }
    }

    void pop_front() {
        (data[begin_pos.first + shift] + begin_pos.second)->~T();
        --sz;
        ++begin_pos.second;
        if (begin_pos.second >= chunk_size) {
            ++begin_pos.first;
            begin_pos.second = 0;
        }

        if (capacity - begin_pos.first - shift < capacity / 4) {
            T** old_data = data;

            data = new T*[capacity / 2 + (capacity % 2 != 0)];

            for (int i = capacity / 2; i < capacity; ++i) {
                data[i - capacity / 2] = old_data[i];
            }

            for (int i = 0; i < capacity / 2; ++i) {
                delete[] reinterpret_cast<int8_t*>(old_data[i]);
            }
            delete[] old_data;

            shift -= capacity / 2;
            capacity /= 2;
        }
    }

    void insert(iterator it, const T& value) {
        T cur_insert_value = value;
        while (it != this->end()) {
            T tmp = *it;
            *it = cur_insert_value;
            cur_insert_value = tmp;

            ++it;
        }

        this->push_back(cur_insert_value);
    }

    void erase(iterator it) {
        while (it + 1 != this->end()) {
            *it = *(it + 1);
            ++it;
        }

        this->pop_back();
    }

    void test() {
        for (int i = 0; i < capacity; ++i) {
            for (int j = 0; j < chunk_size; ++j) {
                std::cout << data[i][j] << ' ';
            }
            std::cout << '\n';
        }
    }

    ~Deque() {
        for (int i = 0; i < capacity; ++i) {
            for (int j = 0; j < chunk_size; ++j) {
                data[i][j].~T();
            }
            delete[] reinterpret_cast<int8_t*>(data[i]);
        }
        delete[] data;
    }
};
