#pragma once
#ifdef _MSC_VER
#include <sstream>
#endif

#include <climits>
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <chrono>

#define CHECK_INT_TYPE(type) \
    static_assert(std::is_same_v<type, int8_t> || \
                  std::is_same_v<type, int16_t> || \
                  std::is_same_v<type, int32_t> || \
                  std::is_same_v<type, int64_t>, \
                  "type must be int8, int16, int32, or int64")

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"

class Integer {
public:
    using word = uint64_t;
    static constexpr word base_digit = 10;
    std::vector<word> words;
    bool is_negative = false;

    Integer() : is_negative(false) {}

    Integer(size_t n, word w, bool negative = false) : words(n, w), is_negative(negative) {}

    [[maybe_unused]] Integer(const word *a, const word *b, bool negative = false) : words(a, b),
                                                                                    is_negative(negative) {}

    Integer(const Integer &a) {
        words = a.words;
        is_negative = a.is_negative;
    }

    Integer &operator=(const Integer &a) = default;


#pragma clang diagnostic push
#pragma ide diagnostic ignored "google-explicit-constructor"
    static constexpr char MIN_NUMERIC_CHARACTER = '0';

    template<class type>
    Integer(type number) {
        CHECK_INT_TYPE(type);
        if (number < 0) {
            is_negative = true;
            number *= -1;
        }
        while (number != 0) {
            words.push_back(number);
            for (int j = 0; j < 64; j++) number >>= 1;
        }
    }

#pragma clang diagnostic pop

    void read(const std::string &s) {
        words.clear();
        this->is_negative = false;
        uint32_t i = 0, end = s.size();
        if (s[i] == '-') {
            ++i;
            is_negative = true;
        }
        for (; i != end; ++i) {
            mul_word(base_digit);
            add_word(s[i] - MIN_NUMERIC_CHARACTER);
        }
    }

    explicit Integer(const std::string &s) : is_negative(false) { read(s); }

    inline size_t size() const { return words.size(); }

    inline word &operator[](size_t i) { return words[i]; }

    const word &operator[](size_t i) const { return words[i]; }

    inline bool is_zero() const { return words.empty(); }

    Integer &set_negative(bool b) {
        this->is_negative = b;
        return *this;
    }

    Integer &truncate() {
        while (!words.empty() && words.back() == 0) words.pop_back();
        return *this;
    }

    size_t bit_size() const {
        if (size() == 0) return 0;
        size_t last = size() - 1;
        size_t result = word_bit_size((*this)[last]) + last * 64;
        return result;
    }

    static int cmp_abs(const Integer &a, const Integer &b) {
        if (a.words == b.words) return 0;
        return std::lexicographical_compare(a.words.rbegin(), a.words.rend(), b.words.rbegin(), b.words.rend()) ? -1
                                                                                                                : 1;
    }

    static inline int cmp(const Integer &a, const Integer &b) {
        if (a.size() == 0 && b.size() == 0) return 0;
        else if (!a.is_negative && !b.is_negative) return +cmp_abs(a, b);
        else if (a.is_negative && b.is_negative) return -cmp_abs(a, b);
        else return a.is_negative && !b.is_negative ? -1 : +1;
    }

    static size_t word_bit_size(word a) {
        for (int i = 64 - 1; i >= 0; i--) if ((a >> i) & 1) return i + 1;
        return 0;
    }

    static word add_carry(word *a, word b) { return (*a += b) < b; }

    static word sub_carry(word *a, word b) {
        word tmp = *a;
        return (*a = tmp - b) > tmp;
    }

    static word word_mul_hi(word a, word b) {
        word a_hi = a >> 32;
        word a_lo = a & UINT_MAX;
        word b_hi = b >> 32;
        word b_lo = b & UINT_MAX;
        word tmp = ((a_lo * b_lo) >> 32) + a_hi * b_lo;
        tmp = (tmp >> 32) + ((a_lo * b_hi + (tmp & UINT_MAX)) >> 32);
        return tmp + a_hi * b_hi;
    }

    static Integer &add_unsigned_overwrite(Integer &a, const Integer &b) {
        size_t i, na = a.size(), nb = b.size(), n = std::max(na, nb);
        a.words.resize(n);
        word carry = 0;
        for (i = 0; i < nb; i++) {
            carry = add_carry(&a[i], carry);
            carry += add_carry(&a[i], b[i]);
        }
        for (; i < n && carry; i++) carry = add_carry(&a[i], carry);
        if (carry) a.words.push_back(carry);
        return a.truncate();
    }

    static Integer &sub_unsigned_overwrite(Integer &a, const Integer &b) {
        size_t i, na = a.size(), nb = b.size();
        word carry = 0;
        for (i = 0; i < nb; i++) {
            carry = sub_carry(&a[i], carry);
            carry += sub_carry(&a[i], b[i]);
        }
        for (; i < na && carry; i++) carry = sub_carry(&a[i], carry);
        return a.truncate();
    }

    static Integer native_multiple(const Integer &a, const Integer &b) {
        size_t na = a.size(), nb = b.size(), nc = na + nb + 1;
        Integer c(nc, 0, a.is_negative ^ b.is_negative), carries(nc, 0);
        for (size_t ia = 0; ia < na; ia++) {
            for (size_t ib = 0; ib < nb; ib++) {
                size_t i = ia + ib, j = i + 1;
                carries[i + 1] += add_carry(&c[i], a[ia] * b[ib]);
                carries[j + 1] += add_carry(&c[j], word_mul_hi(a[ia], b[ib]));
            }
        }
        return add_unsigned_overwrite(c, carries).truncate();
    }

    static Integer karatsuba_multiple(const Integer &a, const Integer &b) {
        size_t na = a.size(), nb = b.size(), n = std::max(na, nb), m2 = n / 2 + (n & 1);
        Integer a_parts[2], b_parts[2];
        split(a, a_parts, 2, m2);
        split(b, b_parts, 2, m2);
        m2 *= 64;
        Integer z0 = a_parts[0] * b_parts[0];
        Integer z1 = (a_parts[0] + a_parts[1]) * (b_parts[0] + b_parts[1]);
        Integer z2 = a_parts[1] * b_parts[1];
        Integer result = z2;
        result <<= m2;
        result += z1 - z2 - z0;
        result <<= m2;
        result += z0;
        return result;
    }

    static Integer multiple(const Integer &a, const Integer &b) {
        return (a.size() > 20 && b.size() > 20 ? karatsuba_multiple(a, b) : native_multiple(a, b));
    }

    static Integer add_signed(const Integer &a, bool a_negative, const Integer &b, bool b_negative) {
        if (a_negative == b_negative) return add_unsigned(a, b).set_negative(a_negative);
        if (cmp_abs(a, b) >= 0) return sub_unsigned(a, b).set_negative(a_negative);
        return sub_unsigned(b, a).set_negative(b_negative);
    }

    Integer &operator>>=(size_t n_bits) {
        if (n_bits == 0) return *this;
        size_t n_words = n_bits / 64;
        if (n_words >= size()) {
            words.resize(0);
            return *this;
        }
        n_bits %= 64;
        if (n_bits == 0) {
            for (size_t i = 0; i < size() - n_words; i++) {
                (*this)[i] = (*this)[i + n_words];
            }
        } else {
            word hi, lo = (*this)[n_words];
            for (size_t i = 0; i < size() - n_words - 1; i++) {
                hi = (*this)[i + n_words + 1];
                (*this)[i] = (hi << (64 - n_bits)) | (lo >> n_bits);
                lo = hi;
            }
            (*this)[size() - n_words - 1] = lo >> n_bits;
        }
        words.resize(size() - n_words);
        return truncate();
    }

    Integer &operator<<=(size_t n_bits) {
        if (n_bits == 0) return *this;
        size_t n_words = n_bits / 64;
        n_bits %= 64;
        size_t old_size = size();
        size_t n = old_size + n_words + (n_bits != 0);
        words.resize(n);
        if (n_bits == 0) {
            for (size_t i = n; i-- > n_words;) (*this)[i] = (*this)[i - n_words];
        } else {
            word lo, hi = 0;
            for (size_t i = n - 1; i > n_words; i--) {
                lo = (*this)[i - n_words - 1];
                (*this)[i] = (hi << n_bits) | (lo >> (64 - n_bits));
                hi = lo;
            }
            (*this)[n_words] = hi << n_bits;
        }
        for (size_t i = 0; i < n_words; i++) (*this)[i] = 0;
        return truncate();
    }

    static void div_mod(const Integer &generator, Integer denominator, Integer &quotient, Integer &remainder) {
        quotient = 0;
        remainder = generator;
        if (cmp_abs(remainder, denominator) >= 0) {
            int n = int(generator.bit_size() - denominator.bit_size());
            denominator <<= n;
            for (; n >= 0; n--) {
                if (cmp_abs(remainder, denominator) >= 0) {
                    sub_unsigned_overwrite(remainder, denominator);
                    quotient.set_bit(n);
                }
                denominator >>= 1;
            }
        }
        quotient.set_negative(generator.is_negative ^ denominator.is_negative);
        remainder.set_negative(generator.is_negative);
    }

    static void div_mod_half_word(const Integer &generator, Integer &quotient, word &remainder) {
        remainder = 0;
        Integer dst(generator.size(), 0);
        for (size_t i = generator.size(); i-- > 0;) {
            word dst_word = 0;
            word src_word = generator[i];
            word parts[2] = {src_word >> 32, src_word & UINT_MAX};
            for (int j: {0, 1}) {
                remainder <<= 32;
                remainder |= parts[j];
                word div_word = remainder / base_digit;
                word mod_word = remainder % base_digit;
                remainder = mod_word;
                dst_word <<= 32;
                dst_word |= div_word;
            }
            dst[i] = dst_word;
        }
        quotient = dst.truncate().set_negative(generator.is_negative);
    }

    static void split(const Integer &a, Integer *parts, size_t n_parts, size_t n) {
        size_t i = 0;
        for (size_t k = 0; k < n_parts; k++) {
            Integer &part = parts[k];
            part.words.resize(n);
            for (size_t j = 0; j < n && i < a.size(); j++) part[j] = a[i++];
            part = part.truncate();
        }
    }

    static Integer div(const Integer &generator, const Integer &denominator) {
        Integer quotient, remainder;
        div_mod(generator, denominator, quotient, remainder);
        return quotient;
    }

    static Integer mod(const Integer &generator, const Integer &denominator) {
        Integer quotient, remainder;
        div_mod(generator, denominator, quotient, remainder);
        return remainder;
    }

    static Integer add_unsigned(const Integer &a, const Integer &b) {
        Integer result(a);
        return add_unsigned_overwrite(result, b);
    }

    static Integer sub_unsigned(const Integer &a, const Integer &b) {
        Integer result(a);
        return sub_unsigned_overwrite(result, b);
    }

    static Integer add(const Integer &a, const Integer &b) {
        Integer result = add_signed(a, a.is_negative, b, b.is_negative);
        return result;
    }

    static Integer sub(const Integer &a, const Integer &b) {
        Integer result = add_signed(a, a.is_negative, b, !b.is_negative);
        return result;
    }

    Integer &set_bit(int i) {
        size_t i_word = i / 64, i_bit = i % 64;
        if (size() <= i_word) words.resize(i_word + 1);
        (*this)[i_word] |= ((word) 1) << i_bit;
        return *this;
    }

    Integer &mul_word(word b) {
        word carry = 0;
        for (size_t i = 0; i < size(); i++) {
            word a = (*this)[i];
            word tmp = a * b;
            carry = add_carry(&tmp, carry);
            carry += word_mul_hi(a, b);
            (*this)[i] = tmp;
        }
        if (carry) words.push_back(carry);
        return truncate();
    }

    Integer &add_word(word carry) {
        if (words.empty()) words.resize(1);
        for (size_t i = 0; i < words.size() && carry; i++) carry = add_carry(&(*this)[i], carry);
        if (carry) words.push_back(carry);
        return truncate();
    }

    [[maybe_unused]]  std::string to_string() const {
        if (words.empty()) return "0";
        std::string result;
        Integer tmp(*this);
        while (tmp.size() != 0) {
            word remainder;
            div_mod_half_word(tmp, tmp, remainder);
            result.push_back(char(MIN_NUMERIC_CHARACTER + remainder));
        }
        if (is_negative) result.push_back('-');
        std::reverse(result.begin(), result.end());
        return result;
    }

    Integer &operator+=(const Integer &b) { return *this = add(*this, b); }

    Integer &operator-=(const Integer &b) { return *this = sub(*this, b); }

    Integer &operator*=(const Integer &b) { return *this = multiple(*this, b); }

    Integer &operator/=(const Integer &b) {
        assert(!b.is_zero());
        return *this = div(*this, b);
    }

    Integer &operator%=(const Integer &b) {
        assert(!b.is_zero());
        return *this = mod(*this, b);
    }

    bool operator==(const Integer &b) const { return cmp(*this, b) == 0; }

    bool operator!=(const Integer &b) const { return cmp(*this, b) != 0; }

    bool operator<=(const Integer &b) const { return cmp(*this, b) <= 0; }

    bool operator>=(const Integer &b) const { return cmp(*this, b) >= 0; }

    bool operator<(const Integer &b) const { return cmp(*this, b) < 0; }

    bool operator>(const Integer &b) const { return cmp(*this, b) > 0; }

    Integer operator+(const Integer &b) const { return add(*this, b); }

    Integer operator-(const Integer &b) const { return sub(*this, b); }

    Integer operator*(const Integer &b) const { return multiple(*this, b); }

    Integer operator/(const Integer &b) const {
        assert(!b.is_zero());
        return div(*this, b);
    }

    Integer operator%(const Integer &b) const {
        assert(!b.is_zero());
        return mod(*this, b);
    }

    Integer operator-() const { return Integer(*this).set_negative(!is_negative); }

    Integer operator>>(size_t n_bits) const { return Integer(*this) >>= n_bits; }

    Integer operator<<(size_t n_bits) const { return Integer(*this) <<= n_bits; }

    friend std::istream &operator>>(std::istream &stream, Integer &v) {
        std::string s;
        stream >> s;
        v.read(s);
        return stream;
    }

    friend std::ostream &operator<<(std::ostream &stream, const Integer &v) {
        std::cout << v.to_string();
        return stream;
    }

    template<class type>
    void operator*=(type number) {
        CHECK_INT_TYPE(type);
        if (number < 0) {
            is_negative = !is_negative;
            number *= -1;
        }
        mul_word((word) number);
    }
};

#pragma clang diagnostic pop
