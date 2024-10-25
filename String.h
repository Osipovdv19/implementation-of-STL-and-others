#include <cstring>
#include <iostream>

class String {
private:
    char* buffer_ = nullptr;
    size_t sz_;
    size_t capacity_;

public:
    String(): sz_(0), capacity_(1) {
        buffer_ = new char[capacity_];
    }

    String(size_t tmp_sz, char value): sz_(tmp_sz), capacity_(1) {
        get_capacity();
        buffer_ = new char[capacity_];
        memset(buffer_, value, sz_);
    }

    String(const String &s) {
        sz_ = s.sz_;
        capacity_ = s.capacity_;
        buffer_ = new char[capacity_];
        memcpy(buffer_, s.buffer_, sz_);
    }

    String(const char* s): sz_(strlen(s)), capacity_(1) {
        get_capacity();
        buffer_ = new char[capacity_];
        memcpy(buffer_, s, sz_);
    }

    void get_capacity() {
        while (capacity_ <= sz_) {
            capacity_ *= 2;
        }
    }

    void swap(String& s) {
        std::swap(buffer_, s.buffer_);
        std::swap(sz_, s.sz_);
        std::swap(capacity_, s.capacity_);
    }

    String& operator=(const String& s) {
        String copy = s;
        swap(copy);
        return *this;
    }

    char operator[](size_t pos) const {
        return buffer_[pos];
    }

    char& operator[](size_t pos) {
        return buffer_[pos];
    }

    size_t length() const {
        return sz_;
    }

    void recalc() {
        char *s = new char[capacity_];
        memcpy(s, buffer_, sz_);
        delete[] buffer_;
        buffer_ = s;
    }

    void push_back(char value) {
        if (capacity_ <= sz_ + 1) {
            capacity_ = capacity_ * 2;
            recalc();
        }
        sz_++;
        buffer_[sz_ - 1] = value;
    }

    void pop_back() {
        sz_--;
    }

    char& front() {
        return buffer_[0];
    }

    char front() const {
        return buffer_[0];
    }

    char& back() {
        return buffer_[sz_ - 1];
    }

    char back() const {
        return buffer_[sz_ - 1];
    }

    String& operator+=(char value) {
        push_back(value);
        return *this;
    }

    String& operator+=(const String& s) {
        get_capacity();
        recalc();
        memcpy(buffer_ + sz_, s.buffer_, s.sz_);
        sz_ += s.sz_;
        return *this;
    }

    String operator+(char b) {
        String copy = *this;
        return copy += b;
    }

    size_t find(const String& s) const {
        if (sz_ < s.sz_) {
            return sz_;
        }
        for (size_t i = 0; i < sz_; ++i) {
            if (sz_ < i + s.sz_) {
                break;
            }
            if (!memcmp(buffer_ + i, s.buffer_, s.sz_)) {
                return i;
            }
        }
        return sz_;
    }

    size_t rfind(const String& s) const {
        if (sz_ < s.sz_) {
            return sz_;
        }
        for (int i = static_cast<int>(sz_) - 1; i >= 0; --i) {
            if (sz_ >= i + s.sz_ && !memcmp(buffer_ + i, s.buffer_, s.sz_)) {
                return i;
            }
        }
        return sz_;
    }

    String substr(size_t start, size_t count) const {
        String s(count, '0');
        memcpy(s.buffer_, buffer_ + start, count);
        return s;
    }

    bool empty() const {
        return !sz_;
    }

    void clear() {
        sz_ = 0;
        capacity_ = 1;
        delete[] buffer_;
        buffer_ = new char[capacity_];
    }

    ~String() {
        delete[] buffer_;
    }
};

String operator+(const String &a, const String &b) {
    String copy_a = a;
    copy_a += b;
    return copy_a;
}

String operator+(char a, const String& b) {
    String copy_a(1, a);
    copy_a += b;
    return copy_a;
}

std::istream &operator>>(std::istream &in, String& s) {
    char value;
    s.clear();
    while (in.get(value) && !isspace(value)) {
        s += value;
    }
    return in;
}

std::ostream &operator<<(std::ostream &out, const String& s) {
    for (size_t i = 0; i < s.length(); ++i) {
        out << s[i];
    }
    return out;
}

bool operator==(const String& a, const String& b) {
    if (a.length() != b.length()) {
        return false;
    }
    for (size_t i = 0; i < a.length(); ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}