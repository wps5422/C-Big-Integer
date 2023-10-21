#include <algorithm>
#include <array>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono>
#include <cassert>

using namespace std;

struct bigint {
    using i64 = long long;
    static constexpr int base = 1e9, base_digits = 9;
    static constexpr i64 base64 = 1e9;
    vector<int> words;
    int sign = 1;

    bigint() = default;

    bigint(i64 v) { *this = v; }

    bigint &operator=(i64 v) {
        sign = v < 0 ? -1 : 1;
        v *= sign;
        words.clear();
        for (; v; v /= base) words.push_back(int(v % base64));
        return *this;
    }

    explicit bigint(const string &s) { read(s); }

    bigint &operator+=(const bigint &other) {
        const int other_words_size = int((other.words).size());
        if (sign == other.sign) {
            for (int i = 0, carry = 0; i < other_words_size || carry; ++i) {
                if (i == int((words).size())) words.push_back(0);
                words[i] += carry + (i < other_words_size ? other.words[i] : 0);
                carry = words[i] >= base;
                if (carry) words[i] -= base;
            }
        } else if (other != 0) *this -= -other;
        return *this;
    }

    friend bigint operator+(bigint a, const bigint &b) { return a += b; }

    bigint &operator-=(const bigint &other) {
        if (sign == other.sign) {
            if ((sign == 1 && *this >= other) || (sign == -1 && *this <= other)) {
                for (int i = 0, carry = 0; i < int((other.words).size()) || carry; ++i) {
                    words[i] -= carry + (i < int((other.words).size()) ? other.words[i] : 0);
                    carry = words[i] < 0;
                    if (carry) words[i] += base;
                }
                trim();
            } else {
                *this = other - *this;
                this->sign = -this->sign;
            }
        } else *this += -other;
        return *this;
    }

    friend bigint operator-(bigint a, const bigint &b) { return a -= b; }

    bigint &operator*=(int v) {
        if (v < 0) sign = -sign, v = -v;
        for (int i = 0, carry = 0; i < int((words).size()) || carry; ++i) {
            if (i == int((words).size())) words.push_back(0);
            i64 cur = (i64) words[i] * v + carry;
            carry = int(cur / base);
            words[i] = int(cur % base);
        }
        trim();
        return *this;
    }

    bigint operator*(int v) const { return bigint(*this) *= v; }

    friend pair<bigint, bigint> divmod(const bigint &a1, const bigint &b1) {
        int norm = base / (b1.words.back() + 1);
        bigint a = a1.abs() * norm, b = b1.abs() * norm, q, r;
        q.words.resize(a.words.size());
        for (int i = int(a.words.size()) - 1; i >= 0; --i) {
            r *= base;
            r += a.words[i];
            i64 s1 = b.words.size() < r.words.size() ? r.words[b.words.size()] : 0;
            i64 s2 = (b.words).size() - 1 < (r.words).size() ? r.words[(b.words).size() - 1] : 0;
            int d = int((s1 * base64 + s2) / b.words.back());
            r -= b * d;
            while (r < 0) r += b, --d;
            q.words[i] = d;
        }
        q.sign = a1.sign * b1.sign;
        r.sign = a1.sign;
        q.trim();
        r.trim();
        return {q, r / norm};
    }

    bigint operator/(const bigint &v) const { return divmod(*this, v).first; }

    bigint operator%(const bigint &v) const { return divmod(*this, v).second; }

    bigint &operator/=(int v) {
        if (v < 0) sign = -sign, v = -v;
        for (i64 i = int((words).size()) - 1, rem = 0; i >= 0; --i) {
            i64 cur = words[i] + rem * base64;
            words[i] = int(cur / v);
            rem = cur % v;
        }
        trim();
        return *this;
    }

    bigint operator/(int v) const { return bigint(*this) /= v; }

    int operator%(int v) const {
        if (v < 0) v = -v;
        i64 m = 0;
        for (int i = int(words.size()) - 1; i >= 0; --i) {
            m = (words[i] + m * base64) % v;
        }
        return int(m * sign);
    }

    bigint &operator*=(const bigint &v) { return *this = *this * v; }

    bigint &operator/=(const bigint &v) { return *this = *this / v; }

    bool operator<(const bigint &v) const {
        if (sign != v.sign) return sign < v.sign;
        if (int((words).size()) != int((v.words).size()))
            return int((words).size()) * sign < int((v.words).size()) * v.sign;
        for (int i = (int((words).size())) - 1; i >= 0; --i)
            if (words[i] != v.words[i]) return words[i] * sign < v.words[i] * sign;
        return false;
    }

    bool operator>(const bigint &v) const { return v < *this; }

    bool operator<=(const bigint &v) const { return !(v < *this); }

    bool operator>=(const bigint &v) const { return !(*this < v); }

    bool operator==(const bigint &v) const { return !(*this < v) && !(v < *this); }

    bool operator!=(const bigint &v) const { return *this < v || v < *this; }

    void trim() {
        while (!words.empty() && words.back() == 0) words.pop_back();
        if (!words.empty()) sign = 1;
    }

    friend bigint operator-(bigint v) {
        if (!v.words.empty()) v.sign = -v.sign;
        return v;
    }

    bigint abs() const { return sign == 1 ? *this : -*this; }

    void read(const string &s) {
        sign = 1;
        words.clear();
        int pos = 0;
        while (pos < int((s).size()) && (s[pos] == '-' || s[pos] == '+')) {
            if (s[pos] == '-') sign = -sign;
            ++pos;
        }
        for (int i = int((s).size()) - 1; i >= pos; i -= base_digits) {
            int x = 0;
            for (int j = max(pos, i - base_digits + 1); j <= i; j++)
                x = x * 10 + s[j] - '0';
            words.push_back(x);
        }
        trim();
    }

    friend istream &operator>>(istream &is, bigint &v) {
        string s;
        is >> s;
        v.read(s);
        return is;
    }

    friend ostream &operator<<(ostream &os, const bigint &v) {
        if (v.sign == -1) os << '-';
        os << (v.words.empty() ? 0 : v.words.back());
        for (int i = (int((v.words).size()) - 1) - 1; i >= 0; --i) {
            os << setw(base_digits) << setfill('0') << v.words[i];
        }
        return os; // pad with zeroes
    }

    string to_string() const {
        string res;
        if (sign == -1) res += '-';
        res += ::to_string(words.empty() ? 0 : words.back());
        for (int i = (int((words).size()) - 1) - 1; i >= 0; --i) {
            string tmp = std::to_string(words[i]);
            res += string(base_digits - int((tmp).size()), '0') + tmp;
        }
        return res;
    }

    static vector<int> convert_base(const vector<int> &a, int old_digits, int new_digits) {
        vector<i64> p(max(old_digits, new_digits) + 1);
        p[0] = 1;
        for (int i = (1); i < (int((p).size())); ++i) p[i] = p[i - 1] * 10;
        vector<int> res;
        i64 cur = 0;
        int cur_digits = 0;
        for (int v: a) {
            cur += v * p[cur_digits];
            cur_digits += old_digits;
            while (cur_digits >= new_digits) {
                res.push_back(int(cur % p[new_digits]));
                cur /= p[new_digits];
                cur_digits -= new_digits;
            }
        }
        res.push_back(int(cur));
        while (!res.empty() && res.back() == 0) res.pop_back();
        return res;
    }

    static constexpr int karatsuba_threshold = 32;

    static vector<i64> karatMul(const vector<i64> &a, const vector<i64> &b) {
        int n = int((a).size());
        vector<i64> res(2 * n);
        if (n <= karatsuba_threshold) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    res[i + j] += a[i] * b[j];
                }
            }
            return res;
        }
        int k = n / 2;
        vector<i64> a1(begin(a), begin(a) + k), a2(k + begin(a), end(a));
        vector<i64> b1(begin(b), begin(b) + k), b2(k + begin(b), end(b));
        const vector<i64> lhs = karatMul(a1, b1), rhs = karatMul(a2, b2);
        const int l_size = int(lhs.size()), r_size = int(rhs.size());
        for (int i = 0; i < (k); ++i) {
            a2[i] += a1[i], b2[i] += b1[i];
        }
        vector<i64> r = karatMul(a2, b2);
        for (int i = 0; i < l_size; ++i) {
            r[i] -= lhs[i];
        }
        for (int i = 0; i < r_size; ++i) {
            r[i] -= rhs[i];
        }
        for (int i = 0; i < (int((r).size())); ++i) {
            res[i + k] += r[i];
        }
        for (int i = 0; i < l_size; ++i) {
            res[i] += lhs[i];
        }
        for (int i = 0; i < r_size; ++i) {
            res[i + n] += rhs[i];
        }
        return res;
    }

    bigint operator*(const bigint &v) const {
        if (min(int((words).size()), int((v.words).size())) < 150) {
            return mul_simple(v);
        }
        bigint res;
        res.sign = sign * v.sign; // should work as long as # of digits isn't too large (> LLONG_MAX/10^{12})
        vector<int> a6 = convert_base(this->words, base_digits, 6); // blocks of 10^6 instead of 10^9
        vector<int> b6 = convert_base(v.words, base_digits, 6);
        vector<i64> a(begin(a6), end(a6)), b(begin(b6), end(b6));
        while (int((a).size()) < int((b).size())) {
            a.push_back(0);
        }
        while (int((b).size()) < int((a).size())) {
            b.push_back(0);
        }
        while (a.size() & (a.size() - 1)) {
            a.push_back(0), b.push_back(0);
        }
        vector<i64> c = karatMul(a, b);
        i64 cur = 0;
        for (int i = 0; i < (int((c).size())); ++i) { // process carries
            cur += c[i];
            res.words.push_back(int(cur % base64));
            cur /= 1000000;
        }
        res.words = convert_base(res.words, 6, base_digits);
        res.trim();
        return res;
    }

    bigint mul_simple(const bigint &v) const {
        bigint res;
        res.sign = sign * v.sign;
        res.words.resize(int((words).size()) + int((v.words).size()));
        for (int i = 0; i < (int((words).size())); ++i)
            if (words[i]) {
                i64 cur = 0;
                for (int j = 0; j < v.words.size() || cur; ++j) {
                    cur += res.words[i + j] + (i64) words[i] * (j < v.words.size() ? v.words[j] : 0);
                    res.words[i + j] = int(cur % base);
                    cur /= base;
                }
            }
        res.trim();
        return res;
    }
};

string hash_string(string &ss) {
    int res = 0;
    for (char c: ss) {
        res = res * 10 + c - '0';
    }
    return "Length: " + to_string(ss.size()) + " Hash: " + to_string(res);
}

void test_mul() {
    {
        string a = "1928374659133333333333333333333333333";
        string b = "192837465913348888888888888888888888333333333333333333333333";
        bigint x(a);
        bigint y(b);
        bigint z = x * y;
        string res = "371862882598789949280919851851851850716253441473328148148147505356595103888888888888888888888889";
        assert(z == bigint(res));
    }
    {
        bigint fac = 1;
        for (int i = 1; i <= 100; ++i) {
            fac = fac * i;
        }
        string fact = "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000";
        string fact_hash = hash_string(fact);
        string tmp = fac.to_string();
        string fac_hash = hash_string(tmp);
        cerr << fact_hash << ' ' << fac_hash << '\n';
        assert(fact_hash == fac_hash);
    }
}


