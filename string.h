#include <iostream>
#include <cstring>

class String {
private:
    size_t capacity = 2, size = 0;
    char* string = nullptr;

    friend bool operator==(const String&, const String&);

    void swap(String& s) {
        std::swap(capacity, s.capacity);
        std::swap(size, s.size);
        std::swap(string, s.string);
    }
public:
    String(const char* source_string) {
        size_t source_string_len = strlen(source_string);
        capacity = (2 > 2 * source_string_len) ? 2 : 2 * source_string_len;
        size = source_string_len;
        string = new char[(2 > 2 * source_string_len) ? 2 : 2 * source_string_len];

        memcpy(string, source_string, size);
    }

    String(size_t size, char c = '\0'): capacity((2 > 2 * size) ? 2 : size), size(size), string(new char[(2 > 2 * size) ? 2 : size]) {
        memset(string, c, size);
    }

    String(char c): capacity(2), size(1), string(new char[2]) {
        string[0] = c;
    }

    String() {}

    String(const String& s): String(s.size, '\0') {
        memcpy(string, s.string, s.size);
    }

    String& operator=(const String& s) {
        String copy = s;
        swap(copy);
        return *this;
    }

    const char& operator[](size_t index) const {
        return string[index];
    }

    char& operator[](size_t index) {
        return string[index];
    }

    size_t length() const {
        return size;
    }

    void push_back(char c) {
        if (size == capacity) {
            char* old_string = string;

            string = new char[(capacity == 0) ? 2 : capacity * 2];
            memcpy(string, old_string, capacity);
            capacity = (capacity == 0) ? 2 : capacity * 2;

            delete[] old_string;
        }

        string[size] = c;
        ++size;
    }

    char pop_back() {
        --size;
        char back_c = string[size];

        if (size < capacity / 4) {
            char* old_string = string;

            string = new char[capacity / 2];
            memcpy(string, old_string, size);
            capacity = capacity / 2;

            delete[] old_string;
        }

        return back_c;
    }

    const char& front() const {
        return string[0];
    }

    char& front() {
        return string[0];
    }

    const char& back() const {
        return string[size - 1];
    }

    char& back() {
        return string[size - 1];
    }

    String& operator+=(const String& s) {
        if (size + s.size < capacity) {
            memcpy(string + size, s.string, s.size);
            size += s.size;
        } else {
            char* old_string = string;
            string = new char[2 * (size + s.size)];

            memcpy(string, old_string, size);
            memcpy(string + size, s.string, s.size);

            delete[] old_string;

            size += s.size;
            capacity = (2 > 2 * (size + s.size)) ? 2 : 2 * (size + s.size);
        }

        return *this;
    }

    String& operator+=(const char& c) {
        this->push_back(c);

        return *this;
    }

    size_t find(const String& substring) const {
        if (size >= substring.length()) {
            for (size_t i = 0; i <= size - substring.length(); ++i) {
                if (memcmp(string + i, substring.string, substring.length()) == 0) {
                    return i;
                }
            }
        }

        return size;
    }

    size_t rfind(const String& substring) const {
        if (size >= substring.length()) {
            for (size_t i = size - substring.length() + 1; i > 0; --i) {
                if (memcmp(string + i - 1, substring.string, substring.length()) == 0) {
                    return i;
                }
            }
        }

        return size;
    }

    String substr(size_t start, size_t count) const {
        char s[count + 1];

        memcpy(s, string + start, count);
        s[count] = '\0';

        return s;
    }

    bool empty() const {
        return size == 0;
    }

    void clear() {
        delete[] string;

        capacity = 2;
        size = 0;
        string = new char[2];
    }

    ~String() {
        delete[] string;
    }
};

String operator+(const String& string1, const String& string2) {
    String copy = string1;
    copy += string2;
    return copy;
}

bool operator==(const String& string1, const String& string2) {
    if (string1.length() != string2.length()) {
        return false;
    } else {
        return memcmp(string1.string, string2.string, string1.length()) == 0;
    }
}

std::ostream& operator<<(std::ostream& out, const String& str) {
    for (size_t i = 0; i < str.length(); ++i) {
        out << str[i];
    }

    return out;
}

std::istream& operator>>(std::istream& in, String& str) {
    str.clear();

    char c;

    while (in.get(c)) {
        if (!isspace(c)) {
            str.push_back(c);
        } else {
            break;
        }
    }

    return in;
}
