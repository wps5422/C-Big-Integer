#include<bits/stdc++.h>
#include "integer.h"

class Integer{
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
  
   [[maybe_unused]] void clr_bit(size_t i) {
        size_t i_word = i / 64;
        size_t i_bit = i % 64;
        if (i_word >= size()) return;
        word mask = 1;
        mask <<= i_bit;
        (*this)[i_word] &= ~mask;
    }
  
  static word word_gcd(word a, word b) {
        while (true) {
            if (a == 0) return b;
            b %= a;
            if (b == 0) return a;
            a %= b;
        }
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
  
    
}
