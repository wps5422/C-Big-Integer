#pragma once

#include <cmath>
#include <cstdint>
#include <climits>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cstdlib>



class Integer {
public:
    using word = unsigned long long int;
    static constexpr int prime_candidates[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    std::vector<word> words;
    bool is_negative = false;

    static word char_to_word(char c) {
        if (c >= '0' && c <= '9') return c - '0';
        else if (c >= 'a' && c <= 'z') return c - 87;
        else if (c >= 'A' && c <= 'Z') return c - 55;
        return ULLONG_MAX;
    }

    static word word_gcd(word a, word b) {
        while (true) {
            if (a == 0) return b;
            b %= a;
            if (b == 0) return a;
            a %= b;
        }
    }

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
    Integer(int64_t i) : is_negative(i < 0) {
        auto u = std::abs(i);
        while (u != 0) {
            words.push_back(u);
            for (int j = 0; j < 64; j++) u >>= 1;
        }
    }
#pragma clang diagnostic pop

    [[maybe_unused]] explicit Integer(const char *c, word base = 10, char **end_ptr = nullptr) : is_negative(false) {
        if (*c == '-') {
            c++;
            is_negative = true;
        }
        // read digits
        for (; *c; c++) {
            mul_word(base);
            word b = char_to_word(*c);
            if (b >= base) break;
            add_word(b);
        }
        if (end_ptr) *end_ptr = (char *) c;
    }

    [[nodiscard]] size_t size() const { return words.size(); }

    word &operator[](size_t i) { return words[i]; }

    const word &operator[](size_t i) const { return words[i]; }

    Integer &set_negative(bool b) {
        this->is_negative = b;
        return *this;
    }

    Integer &truncate() {
        while (size() > 0 && words.back() == 0) words.pop_back();
        return *this;
    }

    [[nodiscard]] size_t bit_size() const {
        if (size() == 0) return 0;
        size_t last = size() - 1;
        size_t result = word_bit_size((*this)[last]) + last * 64;
        return result;
    }

    [[nodiscard]] size_t count_trailing_zeros() const {
        for (size_t i = 0; i < size(); i++) {
            word w = (*this)[i];
            if (w) return word_count_trailing_zeros(w) + i * 64;
        }
        return 0;
    }

    static int cmp_abs(const Integer &a, const Integer &b) {
        size_t na = a.size(), nb = b.size();
        if (na != nb) {
            return na < nb ? -1 : +1;
        }
        for (size_t i = na; i-- > 0;) {
            word wa = a[i], wb = b[i];
            if (wa != wb) {
                return wa < wb ? -1 : +1;
            }
        }
        return 0;
    }

    static int cmp(const Integer &a, const Integer &b) {
        if (a.size() == 0 && b.size() == 0) return 0;
        if (!a.is_negative && !b.is_negative) return +cmp_abs(a, b);
        if (a.is_negative && b.is_negative) return -cmp_abs(a, b);
        return a.is_negative && !b.is_negative ? -1 : +1;
    }

    static size_t word_bit_size(word a) {
        for (int i = 64 - 1; i >= 0; i--) if ((a >> i) & 1) return i + 1;
        return 0;
    }

    static size_t word_count_trailing_zeros(word a) {
        for (size_t i = 0; i < 64; i++) if ((a >> i) & 1) return i;
        return 64;
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
        //assert(cmp_abs(a, b) >= 0);
        size_t i, na = a.size(), nb = b.size();
        word carry = 0;
        for (i = 0; i < nb; i++) {
            carry = sub_carry(&a[i], carry);
            carry += sub_carry(&a[i], b[i]);
        }
        for (; i < na && carry; i++) carry = sub_carry(&a[i], carry);
        //assert(!carry);
        return a.truncate();
    }

    static Integer mul_long(const Integer &a, const Integer &b) {
        size_t na = a.size(), nb = b.size(), nc = na + nb + 1;
        Integer c(nc, 0, a.is_negative ^ b.is_negative), carries(nc, 0);
        for (size_t ia = 0; ia < na; ia++) {
            for (size_t ib = 0; ib < nb; ib++) {
                size_t i = ia + ib, j = i + 1;
                // WARNING: Might overflow if word size is chosen too small
                carries[i + 1] += add_carry(&c[i], a[ia] * b[ib]);
                carries[j + 1] += add_carry(&c[j], word_mul_hi(a[ia], b[ib]));
            }
        }
        return add_unsigned_overwrite(c, carries).truncate();
    }

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
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
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
    static Integer multiple(const Integer &a, const Integer &b) {
        return (a.size() > 20 && b.size() > 20 ? karatsuba_multiple(a, b) : mul_long(a, b));
    }
#pragma clang diagnostic pop

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
            for (size_t i = n; i-- > n_words;) {
                (*this)[i] = (*this)[i - n_words];
            }
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

    static void div_mod_half_word(const Integer &generator, word denominator, Integer &quotient, word &remainder) {
        remainder = 0;
        Integer dst(generator.size(), 0);

        for (size_t i = generator.size(); i-- > 0;) {
            word dst_word = 0;
            word src_word = generator[i];
            word parts[2];
            parts[0] = src_word >> 64 / 2;
            parts[1] = src_word & UINT_MAX;
            {
                remainder <<= 64 / 2;
                remainder |= parts[0];

                word div_word = remainder / denominator;
                word mod_word = remainder % denominator;
                remainder = mod_word;

                dst_word <<= 64 / 2;
                dst_word |= div_word;
            }
            {
                remainder <<= 64 / 2;
                remainder |= parts[1];

                word div_word = remainder / denominator;
                word mod_word = remainder % denominator;
                remainder = mod_word;

                dst_word <<= 64 / 2;
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

    [[maybe_unused]] static Integer gcd(const Integer &a0, const Integer &b0) {

        if (a0.size() == 1 && b0.size() == 1) {
            return {1, word_gcd(a0[0], b0[0])};
        }

        Integer a(a0), b(b0);
        a.is_negative = b.is_negative = false;

        if (a.size() == 0) return b0;
        if (b.size() == 0) return a0;

        size_t n = std::min(a.count_trailing_zeros(), b.count_trailing_zeros());

        a >>= n;
        b >>= n;

        do {
            b >>= b.count_trailing_zeros();
            if (cmp_abs(a, b) > 0) a.words.swap(b.words);
            sub_unsigned_overwrite(b, a);
        } while (b.size() > 0);

        a <<= n;

        return a;
    }

    typedef void (*random_func)(uint8_t *bytes, size_t n_bytes);

    static Integer random_bits(size_t n_bits, random_func func) {
        if (n_bits == 0) return 0;
        size_t partial_bits = n_bits % 64;
        size_t n_words = n_bits / 64 + (partial_bits > 0);
        size_t n_bytes = n_words * sizeof(word);
        Integer result(n_words, 0);
        auto *bytes = (uint8_t *) &result[0];
        func(bytes, n_bytes);
        if (partial_bits) {
            size_t too_many_bits = 64 - partial_bits;
            result.words.back() >>= too_many_bits;
        }
        return result;
    }

    static Integer random_inclusive(const Integer &inclusive, random_func func) {
        size_t n_bits = inclusive.bit_size();
        while (true) {
            Integer result = random_bits(n_bits, func);
            if (result <= inclusive) return result;
        }
    }

    static Integer random_exclusive(const Integer &exclusive, random_func func) {
        size_t n_bits = exclusive.bit_size();
        while (true) {
            Integer result = random_bits(n_bits, func);
            if (result < exclusive) return result;
        }
    }

    [[maybe_unused]] static Integer
    random_second_exclusive(const Integer &inclusive_min_val, const Integer &exclusive_max_val, random_func func) {
        return inclusive_min_val + random_exclusive(exclusive_max_val - inclusive_min_val, func);
    }

    [[maybe_unused]] static Integer
    random_both_inclusive(const Integer &inclusive_min_val, const Integer &inclusive_max_val, random_func func) {
        return inclusive_min_val + random_inclusive(inclusive_max_val - inclusive_min_val, func);
    }

    Integer &set_bit(size_t i) {
        size_t i_word = i / 64;
        size_t i_bit = i % 64;
        if (size() <= i_word) words.resize(i_word + 1);
        (*this)[i_word] |= ((word) 1) << i_bit;
        return *this;
    }

    [[nodiscard]] word get_bit(size_t i) const {
        size_t i_word = i / 64;
        size_t i_bit = i % 64;
        if (i_word >= size()) return 0;
        return ((*this)[i_word] >> i_bit) & 1;
    }

    [[maybe_unused]] void clr_bit(size_t i) {
        size_t i_word = i / 64;
        size_t i_bit = i % 64;
        if (i_word >= size()) return;
        word mask = 1;
        mask <<= i_bit;
        (*this)[i_word] &= ~mask;
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

    Integer &add_word(word carry, size_t i0 = 0) {
        if (i0 >= size()) words.resize(i0 + 1);
        for (size_t i = i0; i < size() && carry; i++) {
            carry = add_carry(&(*this)[i], carry);
        }
        if (carry) words.push_back(carry);
        return truncate();
    }

    [[maybe_unused]] void export_to_array(
            std::vector<char> &text,
            word base = 10,
            const char *alphabet = "0123456789abcdefghijklmnopqrstuvwxyz"
    ) const {
        if (words.empty()) {
            text.push_back('0');
        } else {
            Integer tmp(*this);
            while (tmp.size() > 0) {
                word remainder;
                div_mod_half_word(tmp, base, tmp, remainder);
                text.push_back(alphabet[remainder]);
            }
            if (is_negative) text.push_back('-');
            std::reverse(text.begin(), text.end());
        }
        text.push_back('\0');
    }

    [[maybe_unused]] [[nodiscard]] double to_double() const {
        if (size() == 0) return 0.0;
        double d = 0.0, base = ::pow(2.0, 64);
        for (size_t i = size(); i-- > 0;) {
            d *= base;
            d += (double) (*this)[i]; // ???
        }
        return is_negative ? -d : d;
    }

    [[maybe_unused]] [[nodiscard]] std::string to_string() const {
        if (words.empty()) return "0";
        std::string result;
        Integer tmp(*this);
        while (tmp.size() != 0) {
            word remainder;
            div_mod_half_word(tmp, 10, tmp, remainder);
            result.push_back(char('0' + remainder));
        }
        if (is_negative) result.push_back('-');
        std::reverse(result.begin(), result.end());
        return result;
    }
    [[maybe_unused]] void print() const {
        if (words.empty()){
            putchar('0');
            return;
        }
        std::string result;
        Integer tmp(*this);
        while (tmp.size() != 0) {
            word remainder;
            div_mod_half_word(tmp, 10, tmp, remainder);
            result.push_back(char('0' + remainder));
        }
        if (is_negative){
            putchar('-');
        }
        for(int i = (int)result.size() - 1; i >= 0; --i){
            putchar(result[i] + '0');
        }
    }

    [[maybe_unused]] [[nodiscard]] Integer pow(size_t exponent) const {
        Integer result(1), p(*this);
        for (; exponent; exponent >>= 1) {
            if (exponent & 1) {
                result *= p;
                --exponent;
            }
            p *= p;
        }
        return result;
    }

    [[maybe_unused]] [[nodiscard]] Integer mod_pow(Integer exponent, const Integer &modulus) const {
        Integer result(1), base = (*this) % modulus;
        for (; exponent.size() > 0; exponent >>= 1) {
            if (exponent.get_bit(0)) {
                result = (result * base) % modulus;
            }
            base = (base * base) % modulus;
        }
        return result;
    }

    static Integer power_modulo(Integer base, Integer exponent, const Integer &modulo) {
        Integer result = 1;
        while (exponent != 0) {
            if (exponent % 2 == 1) {
                result *= base;
                result %= modulo;
            }
            exponent >>= 1;
            base *= base;
            base %= modulo;
        }
        return result;
    }

    [[maybe_unused]] [[nodiscard]] bool is_prime() const {
        // make sure n < 3 317 044 064 679 887 385 961 981
        if (*this < 2) return false;
        for (int prime: prime_candidates) {
            if (*this % prime == 0) return *this == prime;
        }
        Integer r = *this - 1;
        int e = 0;
        while (r % 2 == 0) {
            r >>= 1, ++e;
        }
        for (int prime: prime_candidates) {
            Integer x = power_modulo(prime, r, *this);
            for (int t = 0; t < e && x > 1; ++t) {
                Integer y = x * x % *this;
                if (y == 1 && x != *this - 1) return false;
                x = y;
            }
            if (x != 1) return false;
        }
        return true;
    }

    [[maybe_unused]] [[nodiscard]] Integer sqrt() const {
        Integer n = *this;
        int bit = (int) bit_size();
        if (bit & 1) bit ^= 1;
        Integer result = 0;
        for (; bit >= 0; bit -= 2) {
            Integer tmp = result;
            tmp.set_bit(bit);
            if (n >= tmp) {
                n -= tmp;
                result.set_bit(bit + 1);
            }
            result >>= 1;
        }
        return result;
    }


    Integer &operator++() {
        add_word(1);
        return *this;
    }

    Integer &operator+=(const Integer &b) { return *this = add(*this, b); }

    Integer &operator-=(const Integer &b) { return *this = sub(*this, b); }

    Integer &operator*=(const Integer &b) { return *this = multiple(*this, b); }

    Integer &operator/=(const Integer &b) { return *this = div(*this, b); }

    Integer &operator%=(const Integer &b) { return *this = mod(*this, b); }

    bool operator==(const Integer &b) const { return cmp(*this, b) == 0; }

    bool operator!=(const Integer &b) const { return cmp(*this, b) != 0; }

    bool operator<=(const Integer &b) const { return cmp(*this, b) <= 0; }

    bool operator>=(const Integer &b) const { return cmp(*this, b) >= 0; }

    bool operator<(const Integer &b) const { return cmp(*this, b) < 0; }

    bool operator>(const Integer &b) const { return cmp(*this, b) > 0; }

    Integer operator+(const Integer &b) const { return add(*this, b); }

    Integer operator-(const Integer &b) const { return sub(*this, b); }

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"
    Integer operator*(const Integer &b) const { return multiple(*this, b); }
#pragma clang diagnostic pop

    Integer operator/(const Integer &b) const { return div(*this, b); }

    Integer operator%(const Integer &b) const { return mod(*this, b); }

    Integer operator-() const { return Integer(*this).set_negative(!is_negative); }

    Integer operator>>(size_t n_bits) const { return Integer(*this) >>= n_bits; }

    Integer operator<<(size_t n_bits) const { return Integer(*this) <<= n_bits; }
};
