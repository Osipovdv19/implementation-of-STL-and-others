#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <algorithm>
#include <cassert>

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


    BigInteger& operator*=(ll val) {
        if (val < 0) {
            sign_ = -sign_;
            val = -val;
        }
        if (val == 1) {
            return *this;
        }
        long long tmp = 0;
        for (ll i = 0; i < getSize() || tmp; ++i) {
            if (i == (ll) array_.size()) {
                array_.push_back(0);
            }
            long long res = (long long)array_[i] * (long long) val + tmp;
            tmp = (ll) (res / base);
            array_[i] = (ll) (res % base);
        }
        removeZeros();
        return *this;
    }

    BigInteger operator*(ll val) {
        BigInteger copy = *this;
        copy *= val;
        return copy;
    }

    BigInteger& operator/=(ll val) {
        if (isZero() || val == 1) {
            return *this;
        }
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
        removeZeros();
        return *this;
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
        if (v.getSize() < 3) {
            copy /= v.sign_ * (v.getDigits(0) + v.getDigits(1) * base);
        }
        else {
            copy /= v;
        }
        return copy;
    }

    BigInteger operator*(const BigInteger &v) {
        BigInteger copy = *this;
        if (v.getSize() < 3) {
            copy *= v.sign_ * (v.getDigits(0) + v.getDigits(1) * base);
        }
        else {
            copy *= v;
        }
        return copy; //
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
            ans.second.removeZeros();
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

    bool isEven() const {
        return getDigits(0) % 2 == 0;
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

BigInteger gcd(const BigInteger& v, const BigInteger& u) {
    BigInteger a = v;
    BigInteger b = u;
    BigInteger ans = 1;
    if (a < 0) {
        a = -a;
    }
    while (true) {
        if (a == 1 || b == 1) {
            break;
        }
        else if (a.isZero()) {
            ans *= b;
            break;
        }
        else if (b.isZero()) {
            ans *= a;
            break;
        }
        else if (a.isEven() && b.isEven()) {
            ans *= 2;
            a /= 2;
            b /= 2;
        }
        else if (a.isEven()) {
            a /= 2;
        }
        else if (b.isEven()) {
            b /= 2;
        }
        else if (a <= b) {
            b -= a;
            b /= 2;
        }
        else {
            a -= b;
            a /= 2;
        }
    }
    return ans;
}

class Rational {
private:
    BigInteger a;
    BigInteger b;

public:
    Rational(): a(0), b(1) {}

    Rational(ll value): a(value), b(1) {}

    Rational(const BigInteger& n): a(n), b(1) {}

    Rational& operator=(const Rational &v) {
        a = v.a;
        b = v.b;
        return *this;
    }

    Rational& operator+=(const Rational &v) {
        a = a * v.b + b * v.a;
        b *= v.b;
        getOk();
        return *this;
    }

    Rational& operator-=(const Rational &v) {
        a = a * v.b - b * v.a;
        b *= v.b;
        getOk();
        return *this;
    }

    Rational& operator*=(const Rational &v) {
        a *= v.a;
        b *= v.b;
        getOk();
        return *this;
    }

    Rational& operator/=(const Rational &v) {
        a *= v.b;
        b *= v.a;
        getOk();
        return *this;
    }

    Rational& operator/=(ll val) {
        b *= val;
        getOk();
        return *this;
    }

    Rational operator/(ll val) {
        Rational copy = *this;
        copy /= val;
        return copy;
    }

    Rational& operator*=(ll val) {
        a *= val;
        getOk();
        return *this;
    }

    Rational operator*(ll val) {
        Rational copy = *this;
        copy *= val;
        return copy;
    }

    string toString() {
        string ans = "";
        getOk();
        ans += a.toString();
        if (b != 1) {
            ans += '/';
            ans += b.toString();
        }
        return ans;
    }

    string asDecimal(size_t precision = 0) const {
        BigInteger v = a;
        string ans = "";
        string tmp = "";
        for (size_t i = 0; i < precision; ++i) {
            v *= 10;
        }
        v /= b;
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
        getOkSign();
        BigInteger g = gcd(a, b);
        a /= g;
        b /= g;
        getOkSign();
    }

    void getOkSign() {
        if (b.getSign() < 0) {
            a.setSign(a.getSign() * (-1));
            b.setSign(1);
            a.removeZeros();
            b.removeZeros();
        }
    }

    bool operator<(const Rational &u) const{
        return a * u.b < b * u.a;
    }

    Rational operator-() const{
        Rational ans = *this;
        ans.a = -ans.a;
        return ans;
    }

    BigInteger& getA() {
        return a;
    }

    BigInteger& getB() {
        return b;
    }

    explicit operator double() const {
        return stod(asDecimal(8));
    }
};

istream& operator>>(istream& in, Rational& v) {
    in >> v.getA();
    v.getB() = 1;
    return in;
}

ostream& operator<<(ostream& out, Rational& v) {
//    out << v.asDecimal(11);
    v.getOk();
//    out << v.getA() << '/' << v.getB() << ' ';
    out << v.asDecimal(11) << ' ';
    return out;
}

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

template<ll N, ll K>
struct CheckIsPrime {
    static const bool value = N % K == 0 ? false : CheckIsPrime<N, K - 1>::value;
};

template<ll N>
struct CheckIsPrime<N, 2> {
    static const bool value = N % 2 != 0;
};

template <ll N>
struct IsPrime {
    static const bool value = CheckIsPrime<N, min(N - 1, (ll)sqrt(N) + 10)>::value;
};

template<>
struct IsPrime<0> {
    static const bool value = false;
};

template<>
struct IsPrime<1> {
    static const bool value = false;
};

template<>
struct IsPrime<2> {
    static const bool value = true;
};

template <size_t N>
class Residue {
private:
    ll val_;
public:
    Residue(ll val = 0): val_(val) {
        relax();
    };

    static int binPow(ll a, ll k, ll M) {
        ll ans = 1;
        a %= M;
        while (k) {
            if (k & 1) {
                ans = ans * a % M;
            }
            a = a * a % M;
            k >>= 1;
        }
        return ans;
    }

    void relax() {
        ll m = N;
        if (val_ >= 0) {
            val_ = val_ % m;
        }
        else {
            val_ = -val_;
            val_ = (m - val_ % m) % m;
        }
    }

    Residue& operator+=(const Residue &u) {
        val_ += u.val_;
        relax();
        return *this;
    }

    Residue& operator*=(const Residue &u) {
        val_ *= u.val_;
        relax();
        return *this;
    }

    Residue& operator-=(const Residue &u) {
        val_ -= u.val_;
        relax();
        return *this;
    }

    Residue& operator/=(const Residue &u) {
        static_assert(IsPrime<N>::value, "Error with /=");
        val_ *= binPow(u.val_, (ll)N - 2, N);
        relax();
        return *this;
    }

    explicit operator int() const {
        return val_;
    }
};

template <size_t N>
bool operator==(const Residue<N> &a, const Residue<N> &b) {
    return int(a) == int(b);
}

template <size_t N>
bool operator!=(const Residue<N> &a, const Residue<N> &b) {
    return int(a) != int(b);
}

template <size_t N>
Residue<N> operator+(const Residue<N> &a, const Residue<N> &b) {
    Residue<N> copy = a;
    copy += b;
    return copy;
}

template <size_t N>
Residue<N> operator-(const Residue<N> &a, const Residue<N> &b) {
    Residue<N> copy = a;
    copy -= b;
    return copy;
}

template <size_t N>
Residue<N> operator/(const Residue<N> &a, const Residue<N> &b) {
    Residue<N> copy = a;
    copy /= b;
    return copy;
}

template <size_t N>
Residue<N> operator*(const Residue<N> &a, const Residue<N> &b) {
    Residue<N> copy = a;
    copy *= b;
    return copy;
}

template <size_t N>
ostream& operator<<(ostream&out, const Residue<N> &a) {
    out << int(a);
    return out;
}

template<size_t N, size_t M, typename Field = Rational>
class Matrix {
private:
    vector<vector<Field>> mt_;

public:
    Matrix(Field val = 1) {
        mt_.assign(N, vector<Field>(M, Field(0)));
        if (N != M) {
            return;
        }
        for (size_t i = 0; i < N; ++i) {
            mt_[i][i] = val;
        }
    }

    Matrix(const vector<vector<Field>> &lst) : mt_(lst) {
    }

    Matrix(const initializer_list<initializer_list<Field>> &lst) {
        mt_.resize(N, vector<Field>(M));
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                mt_[i][j] = 0;
            }
        }
        ll i = 0, j;
        for (auto v : lst) {
            j = 0;
            for (auto u : v) {
                mt_[i][j] = u;
                j++;
            }
            i++;
        }
    }

    vector<Field>& operator[](ll i) {
        return mt_[i];
    }

    const vector<Field>& operator[](ll i) const {
        return mt_[i];
    }

    bool operator==(const Matrix<N, M, Field> &u) const {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                if (mt_[i][j] != u[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const Matrix<N, M, Field> &u) const {
        return !(*this == u);
    }

    Matrix<N, M, Field>& operator=(const Matrix<N, M, Field> &u) {
        mt_ = u.mt_;
        return *this;
    }

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field> u) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                mt_[i][j] += u[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field> &u) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                mt_[i][j] -= u[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field>& operator*=(const Matrix<N, M, Field> &u) {
        *this = *this * u;
        return *this;
    }

    Matrix<N, M, Field>& operator*=(Field u) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                mt_[i][j] *= u;
            }
        }
        return *this;
    }

    template<size_t P>
    Matrix<N, P, Field> operator*(const Matrix<M, P, Field> &u) const {
        Matrix<N, P, Field> ans(Field(0));
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < P; ++j) {
                for (size_t k = 0; k < M; ++k) {
                    ans[i][j] += mt_[i][k] * u[k][j];
                }
            }
        }
        return ans;
    }

    Matrix<N, M, Field> operator*(Field val) const{
        Matrix<N, M, Field> ans = *this;
        return ans *= val;
    }

    Matrix<N, M, Field> operator+(const Matrix<N, M, Field> &u) const{
        Matrix<N, M, Field> ans = *this;
        return ans += u;
    }

    Matrix<N, M, Field> operator-(const Matrix<N, M, Field> &u) const{
        Matrix<N, M, Field> ans = *this;
        return ans -= u;
    }

    Matrix<M, N, Field> transposed() const{
        Matrix<M, N, Field> ans; //
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                ans[j][i] = mt_[i][j];
            }
        }
        return ans;
    }

    pair<Field, Matrix<N, M, Field>> getGauss() const {
        Matrix<N, M, Field> v = *this;
        Field ans(1);
        size_t k;
        for (size_t j = 0; j < M; ++j) {
            k = 0;
            for (size_t i = j; i < N; ++i)
                if (v[i][j] != Field(0)) {
                    k = i + 1;
                    break;
                }
            if (!k) {
                ans = 0;
                continue;
            }
            k--;
            swap(v[j], v[k]);
            if (j != k) {
                ans *= -1;
            }
            ans *= v[j][j];
            for (size_t i = 0; i < M; ++i) {
                if (j == i) {
                    continue;
                }
                v[j][i] /= v[j][j];
            }
            v[j][j] /= v[j][j];
            for (size_t i = 0; i < N; ++i) {
                Field u = v[i][j];
                if (i == j || u == Field(0)) {
                    continue;
                }
                for (size_t t = 0; t < M; ++t) {
                    v[i][t] -= v[j][t] * u;
                }
            }
        }
        return {ans, v};
    }

    Field det() const {
        static_assert(N == M);
        return getGauss().first;
    }

    size_t rank() const{
        size_t ans = N;
        for (auto &v : getGauss().second.mt_) {
            ll cur = 0;
            for (auto u : v) {
                cur += (u == Field(0));
            }
            ans -= (cur == M);
        }
        return ans;
    }

    Field trace() const {
        static_assert(N == M, "trace/n");
        Field ans = 0;
        for (size_t i = 0; i < N; ++i) {
            ans += mt_[i][i];
        }
        return ans;
    }

    void invert() {
        static_assert(N == M);
        Matrix<N, M * 2, Field> ans;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                ans[i][j] = mt_[i][j];
                ans[i][j + N] = 0;
            }
            ans[i][i + N] = 1;
        }
        ans = ans.getGauss().second;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                mt_[i][j] = ans[i][j + M];
            }
        }
    }

    Matrix<N, M, Field> inverted() const {
        static_assert(N == M);
        Matrix<N, M, Field> ans = *this;
        ans.invert();
        return ans;
    }

    vector<Field> getRow(size_t i) const {
        return mt_[i];
    }

    vector<Field> getColumn(size_t i) const {
        vector<Field> ans;
        for (size_t j = 0; j < M; ++j) {
            ans.push_back(mt_[j][i]);
        }
        return ans;
    }
};

template<size_t N, size_t M, typename Field=Rational>
Matrix<N, M, Field> operator*(Field val, const Matrix<N, M, Field> &u){
    return u * val;
}

template<size_t N, size_t M, typename Field=Rational>
istream& operator>>(istream& in, Matrix<N, M, Field>& v) {
    for (ll i = 0, num; i < N; ++i) {
        for (ll j = 0; j < M; ++j) {
            in >> num;
            v[i][j] = num;
        }
    }
    return in;
}

template<size_t N, size_t M, typename Field=Rational>
ostream& operator<<(ostream& out, Matrix<N, M, Field>& v) {
    for (ll i = 0; i < N; ++i) {
        for (ll j = 0; j < M; ++j) {
            out << v[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

template<size_t N, typename Field=Rational>
using SquareMatrix = Matrix<N, N, Field>;