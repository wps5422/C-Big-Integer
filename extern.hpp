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
}
