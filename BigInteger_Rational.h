#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <algorithm>

using namespace std;

typedef complex<long double> cld;
typedef int ll;

const ll base = 1e4;
const ll base_digits = 4;
const long double pi = acos(-1);

class BigInteger {
private:
    ll sign_; //1 - not negative, -1 - negative
    vector<ll> array_;

public:
    BigInteger() {
        sign_ = 1;
        array_ = {};
    }

    void readInt(ll number) {
        clear();
        sign_ = 1;
        if (number < 0) {
            sign_ = -1;
            number *= -1;
        }
        if (number) {
            while (number) {
                array_.push_back(number%base);
                number /= base;
            }
        }
        removeZeros();
    }

    void readString(string s) {
        clear();
        ll pos = 0;
        sign_ = 1;
        if (s[0] == '-') {
            sign_ = -1;
            pos++;
        }
        if (s[0] == '+') {
            pos++;
        }
        for (ll i = s.size() - 1; i >= pos; i -= base_digits) {
            ll tmp = 0;
            for (ll j = max(pos, i - base_digits + 1); j <= i; ++j) {
                tmp = tmp * 10 + s[j] - '0';
            }
            array_.push_back(tmp);
        }
        removeZeros();
    }

    string toString() const {
        string s = "";
        if (isZero()) {
            s = '0';
            return s;
        }
        if (sign_ == -1) {
            s = '-';
        }
        bool flag = false;
        for (ll j = getSize() - 1, number; j >= 0; j--) {
            number = array_[j];
            for (ll i = base / 10, digit; i >= 1; i /= 10) {
                digit = number / i %10;
                if (!flag && !digit) {
                    continue;
                }
                flag = true;
                s += char('0' + digit);
            }
        }
        return s;
    }

    BigInteger(ll number) {
        readInt(number);
    }

    BigInteger& operator=(ll number) {
        readInt(number);
        return *this;
    }

    bool operator<(const BigInteger &v) const {
        if (sign_ != v.sign_) {
            return sign_ < v.sign_;
        }
        if (getSize() != v.getSize()) {
            return sign_*getSize() < v.getSign()*v.getSize();
        }
        for (ll i = getSize() - 1; i >= 0; i--) {
            if (array_[i] != v.getDigits(i)) {
                return sign_*array_[i] < sign_*v.getDigits(i);
            }
        }
        return false;
    }

    bool operator>(const BigInteger &v) const {
        return v < *this;
    }

    bool operator==(const BigInteger &v) const {
        return !(v < *this) && !(*this < v);
    }

    bool operator>=(const BigInteger &v) const {
        return !(*this < v);
    }

    bool operator<=(const BigInteger &v) const {
        return !(*this > v);
    }

    bool operator!=(const BigInteger &v) const {
        return (v < *this) || (*this < v);
    }

    ll compareAbs(const BigInteger &v) {
        if (getSize() != v.getSize()) {
            return getSize() < v.getSize() ? -1 : 1;
        }
        for (ll i = getSize() - 1; i >= 0; i--) {
            if (array_[i] != v.getDigits(i)) {
                return array_[i] < v.getDigits(i) ? -1 : 1;
            }
        }
        return 0;
    }

    void addWithoutSign(const BigInteger &v) {
        while (getSize() < v.getSize()) {
            array_.push_back(0);
        }
        for (ll i = 0, cur_sum = 0; i < max(getSize(), v.getSize()) || cur_sum; ++i) {
            if (i == getSize()) {
                array_.push_back(0);
            }
            array_[i] += cur_sum + v.getDigits(i);
            cur_sum = 0;
            if (array_[i] >= base) {
                array_[i] -= base;
                cur_sum = 1;
            }
        }
        removeZeros();
    }

    void subWithoutSign(const BigInteger &v) {
        for (ll i = 0, cur_sum = 0; i < v.getSize() || cur_sum; ++i) {
            array_[i] += cur_sum - v.getDigits(i);
            cur_sum = 0;
            if (array_[i] < 0) {
                array_[i] += base;
                cur_sum = -1;
            }
        }
        removeZeros();
    }

    void sum(const BigInteger &v) {
        if (sign_ == v.getSign()) {
            addWithoutSign(v);
        }
        else {
            if (compareAbs(v) >= 0) {
                subWithoutSign(v);
            }
            else {
                BigInteger tmp = v;
                swap(tmp, *this);
                subWithoutSign(tmp);
            }
        }
    }

    BigInteger& operator+=(const BigInteger &v) {
        sum(v);
        return *this;
    }

    void sub(const BigInteger &v) {
        if (sign_ != v.getSign()) {
            addWithoutSign(v);
        }
        else {
            if (compareAbs(v) >= 0) {
                subWithoutSign(v);
            }
            else {
                BigInteger tmp = v;
                swap(*this, tmp);
                subWithoutSign(tmp);
                sign_ *= -1;
            }
        }
    }

    BigInteger& operator-=(const BigInteger &v) {
        sub(v);
        return *this;
    }

    void recalcBase(vector<ll> &mas, ll from, ll to) {
        if (from == to) {
            return;
        }
        vector<long long> cur(max(from, to) + 1, 0);
        vector<ll> ans;
        cur[0] = 1;
        for (ll i = 0; i < (ll)cur.size() - 1; ++i) {
            cur[i + 1] = cur[i] * 10;
        }
        ll digits = 0;
        long long tmp = 0;
        for (ll i = 0; i < (ll)mas.size(); ++i) {
            tmp += mas[i] * cur[digits];
            digits += from;
            while (to <= digits) {
                ans.push_back((long long)(tmp % cur[to]));
                tmp /= cur[to];
                digits -= to;
            }
        }
        if (tmp) {
            ans.push_back(tmp);
        }
        while (!ans.empty() && !ans.back()) {
            ans.pop_back();
        }
        mas = ans;
    }

    void getFft(vector<cld> &a, bool invert) const {
        ll n = (ll)a.size();
        for (ll i = 1, j = 0; i < n; ++i) {
            ll bit = n >> 1;
            for (; j >= bit; bit >>= 1)
                j -= bit;
            j += bit;
            if (i < j) {
                swap(a[i], a[j]);
            }
        }
        for (ll len = 2; len <= n; len <<= 1) {
            long double ang = 2 * pi / (long double)len * (invert ? -1 : 1);
            cld wlen(cos(ang), sin(ang));
            for (ll i = 0; i < n; i += len) {
                cld w(1);
                for (ll j = 0; j < len / 2; ++j) {
                    cld u = a[i + j],  v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len/2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) {
            for (ll i = 0; i < n; ++i) {
                a[i] /= n;
            }
        }
    }

    void fft(const vector<ll> &a, const vector<ll> &b, vector<ll> &ans) const {
        ll len = 1;
        vector<cld> cur_a(a.begin(), a.end()), cur_b(b.begin(), b.end());
        while ((ll)max(a.size(), b.size()) > len) {
            len *= 2;
        }
        len *= 2;
        cur_a.resize(len);
        cur_b.resize(len);
        getFft(cur_a, 0);
        getFft(cur_b, 0);
        for (ll i = 0; i < len; ++i) {
            cur_a[i] *= cur_b[i];
        }
        getFft(cur_a, 1);
        ans.resize(len);
        for (size_t i = 0, carry = 0; i < ans.size(); ++i) {
            long long cur = (long long)(cur_a[i].real() + 0.5) + carry;
            carry = cur / 10000;
            if (carry > 0 && i + 1 >= ans.size()) {
                ans.push_back(0);
            }
            ans[i] = cur % 10000;
        }
    }

    BigInteger& operator*=(const BigInteger &v) {
        BigInteger copy_v = v;
        sign_ = v.sign_*sign_;
        fft(array_, copy_v.array_ ,array_);
        removeZeros();
        return *this;
    }


    void operator*=(ll val) {
        if (val < 0) {
            sign_ = -sign_;
            val = -val;
        }
        long long tmp = 0;
        for (ll i = 0; i < getSize() || tmp; ++i) {
            if (i == (ll)array_.size()) {
                array_.push_back(0);
            }
            long long res = (long long)array_[i] * (long long) val + tmp;
            tmp = (ll) (res / base);
            array_[i] = (ll) (res % base);
        }
        removeZeros();
    }

    BigInteger operator*(ll val) {
        BigInteger copy = *this;
        copy *= val;
        return copy;
    }

    void operator/=(ll val) {
        if (isZero())
            return;
        if (val < 0) {
            sign_ *= -1;
            val *= -1;
        }
        long long cur = 0;
        for (ll i = (ll)getSize() - 1; i >= 0; --i) {
            long long tmp = (long long)base * (long long)cur + array_[i];
            cur = (long long)tmp % (long long)val;
            array_[i] = (long long)tmp / (long long)val;
        }
    }

    BigInteger operator/(ll val) {
        BigInteger copy = *this;
        copy /= val;
        return copy;
    }

    BigInteger operator-(const BigInteger &v) {
        BigInteger copy = *this;
        copy -= v;
        return copy;
    }

    BigInteger operator/(const BigInteger &v) {
        BigInteger copy = *this;
        copy /= v;
        return copy;
    }

    BigInteger operator*(const BigInteger &v) {
        BigInteger cop = *this;
        cop *= v;
        return cop;
    }

    pair<BigInteger, BigInteger> getDivMod(const BigInteger &v) {
        long long val = base / (v.array_.back() + 1);
        BigInteger a = abs() * val, b = v.abs() * val;
        BigInteger q(0), r(0);
        q.array_.resize(a.getSize());
        for (ll i = a.getSize() - 1; i >= 0; i--) {
            r *= base;
            r += a.array_[i];
            long long fst = r.getSize() <= b.getSize() ? 0 : r.getDigits(b.getSize());
            long long snd = r.getSize() <= b.getSize() - 1 ? 0 : r.getDigits(b.getSize() - 1);
            long long d = ((long long) base * fst + snd) / b.array_.back();
            r -= b * d;
            while (r < 0) {
                r += b;
                --d;
            }
            q.array_[i] = d;
        }
        q.sign_ = sign_ * v.getSign();
        r.sign_ = sign_;
        q.removeZeros();
        r.removeZeros();
        pair<BigInteger, BigInteger> ans = make_pair(q, r / val);
        if (ans.second < 0) {
            ans.second += v;
        }
        return ans;
    }

    BigInteger& operator/=(const BigInteger &v) {
        if (!isZero()) {
            *this = getDivMod(v).first;
            removeZeros();
        }
        return *this;
    }

    BigInteger& operator%=(const BigInteger &v) {
        if (!isZero()) {
            *this = (*this - (*this / v) * v);
        }
        return *this;
    }

    BigInteger& operator++() {
        sum(BigInteger(1));
        return *this;
    }

    BigInteger& operator--() {
        sub(BigInteger(1));
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger ans = *this;
        ++(*this);
        return ans;
    }

    BigInteger operator--(int) {
        BigInteger ans = *this;
        --(*this);
        return ans;
    }

    ll getSize() const {
        return array_.size();
    }

    ll getSign() const {
        return sign_;
    }

    void setSign(ll value) {
        sign_ = value;
    }

    ll getDigits(ll i) const {
        return i < getSize() ? array_[i] : 0;
    }

    void clear() {
        sign_ = 1;
        array_.clear();
    }

    BigInteger abs() const {
        BigInteger ans = *this;
        ans.sign_ = 1;
        return ans;
    }

    void removeZeros() {
        while (!array_.empty() && !array_.back()) {
            array_.pop_back();
        }
        if (array_.empty()) {
            sign_ = 1;
        }
    }

    bool isZero() const {
        return array_.empty();
    }

    explicit operator bool() const {
        return !isZero();
    }

    void swapSize() {
        sign_ *= -1;
    }
};

BigInteger operator-(const BigInteger &a) {
    BigInteger ans = a;
    if (!a.isZero()) {
        ans.swapSize();
    }
    return ans;
}

BigInteger operator*(const BigInteger &a, const BigInteger &b) {
    BigInteger copy_a = a;
    copy_a *= b;
    return copy_a;
}

BigInteger operator+(const BigInteger &a, const BigInteger &b) {
    BigInteger copy_a = a;
    copy_a += b;
    return copy_a;
}

BigInteger operator-(const BigInteger &a, const BigInteger &b) {
    BigInteger copy_a = a;
    copy_a -= b;
    return copy_a;
}

BigInteger operator/(const BigInteger &a, const BigInteger &b) {
    BigInteger copy_a = a;
    copy_a /= b;
    return copy_a;
}

BigInteger operator%(const BigInteger &a, const BigInteger &b) {
    BigInteger copy_a = a;
    copy_a %= b;
    return copy_a;
}

istream& operator>>(std::istream &in, BigInteger& v) {
    string s;
    in >> s;
    v.readString(s);
    return in;
}

ostream& operator<<(std::ostream &out, const BigInteger& v) {
    out << v.toString();
    return out;
}

BigInteger gcd(const BigInteger &a, const BigInteger &b) {
    BigInteger copy_a = a, copy_b = b;
    while (!copy_a.isZero()) {
        copy_b %= copy_a;
        swap(copy_a, copy_b);
    }
    return copy_b;
}

class Rational {
private:
    BigInteger a_, b_;

public:
    Rational(): a_(0), b_(1) {}

    Rational(ll value): a_(value), b_(1) {}

    Rational(const BigInteger& n): a_(n), b_(1) {}

    Rational& operator=(const Rational &v) {
        a_ = v.a_;
        b_ = v.b_;
        return *this;
    }

    Rational& operator+=(const Rational &v) {
        a_ = a_ * v.b_ + b_ * v.a_;
        b_ *= v.b_;
        return *this;
    }

    Rational& operator-=(const Rational &v) {
        a_ = a_ * v.b_ - b_ * v.a_;
        b_ *= v.b_;
        return *this;
    }

    Rational& operator*=(const Rational &v) {
        a_ *= v.a_;
        b_ *= v.b_;
        return *this;
    }

    Rational& operator/=(const Rational &v) {
        a_ *= v.b_;
        b_ *= v.a_;
        return *this;
    }

    string toString() {
        string ans = "";
        getOk();
        ans += a_.toString();
        if (b_ != 1) {
            ans += '/';
            ans += b_.toString();
        }
        return ans;
    }

    string asDecimal(size_t precision = 0) {
        BigInteger v = a_;
        string ans = "", tmp = "";
        for (size_t i = 0; i < precision; ++i) {
            v *= 10;
        }
        v /= b_;
        ans += v.toString();
        reverse(ans.begin(), ans.end());
        if (ans.back() == '-') {
            tmp += '-';
            ans.pop_back();
        }
        while (ans.size() <= precision) {
            ans += '0';
        }
        if (precision > 0) {
            ans.insert(ans.begin() + precision, '.');
        }
        while (ans.size() > 1 && ans[ans.size() - 2] != '.' && ans.back() == '0') {
            ans.pop_back();
        }
        ans += tmp;
        reverse(ans.begin(), ans.end());
        return ans;
    }

    void getOk() {
        BigInteger g = gcd(a_, b_);
        a_ /= g;
        b_ /= g;
        getOkSign();
    }

    void getOkSign() {
        if (b_.getSign() < 0) {
            a_.setSign(a_.getSign() * (-1));
            b_.setSign(1);
        }
    }

    bool operator<(const Rational &u) const{
        return a_ * u.b_ < b_ * u.a_;
    }

    Rational operator-() const{
        Rational ans = *this;
        ans.a_ = -ans.a_;
        return ans;
    }
};

bool operator>(const Rational &v, const Rational &u) {
    return u < v;
}

bool operator<=(const Rational &v, const Rational &u) {
    return !(v > u);
}

bool operator>=(const Rational &v, const Rational &u) {
    return !(v < u);
}

bool operator==(const Rational &v, const Rational &u) {
    return v >= u && v <= u;
}

bool operator!=(const Rational &v, const Rational &u) {
    return !(v == u);
}

Rational operator+(const Rational &v, const Rational &u) {
    Rational copy_v = v;
    copy_v += u;
    return copy_v;
}

Rational operator-(const Rational &v, const Rational &u) {
    Rational copy_v = v;
    copy_v -= u;
    return copy_v;
}

Rational operator*(const Rational &v, const Rational &u) {
    Rational copy_v = v;
    copy_v *= u;
    return copy_v;
}

Rational operator/(const Rational &v, const Rational &u) {
    Rational copy_v = v;
    copy_v /= u;
    return copy_v;
}