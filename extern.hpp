#include<bits/stdc++.h>
#include "integer.h"

class Integer{
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
}
