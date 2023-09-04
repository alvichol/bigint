#include "big_integer.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

static const std::size_t DIGIT_LENGTH = 32;
static const uint32_t BASE = 10;
static const uint32_t LONG_BASE = 1e9;
static const std::size_t DEC_SHIFT = 9;
static const uint32_t DIG_MAX = std::numeric_limits<uint32_t>::max();
static const int64_t SMIN64 = std::numeric_limits<int64_t>::min();
static const int64_t SMAX64 = std::numeric_limits<int64_t>::max();
static const int32_t SMIN32 = std::numeric_limits<int32_t>::min();
static const int32_t SMAX32 = std::numeric_limits<int32_t>::max();
static const int64_t SHIFT_ONE = 1LL << DIGIT_LENGTH;

std::size_t big_integer::len() const noexcept {
  return digits.size();
}

bool big_integer::is_zero() const noexcept {
  return len() == 0;
}

void big_integer::trim_zero() noexcept {
  while (len() > 0 && digits[len() - 1] == 0) {
    digits.pop_back();
  }
}

bool big_integer::abs_lt(const big_integer& other) const {
  if (len() != other.len()) {
    return len() < other.len();
  }
  if (len() == 0) {
    return false;
  }
  return std::lexicographical_compare(digits.rbegin(), digits.rend(), other.digits.rbegin(), other.digits.rend());
}

void big_integer::swap(big_integer& other) noexcept {
  std::swap(digits, other.digits);
  std::swap(sign, other.sign);
}

void big_integer::adjust_size(const big_integer& other) {
  if (other.len() > len()) {
    digits.resize(other.len(), 0);
  }
}

void big_integer::position_increment(std::size_t pos) {
  digits.resize(len() + 1, 0);
  while (digits[pos] == DIG_MAX) {
    digits[pos] = 0;
    pos++;
  }
  digits[pos]++;
  trim_zero();
}

void big_integer::position_decrement(std::size_t pos) {
  while (digits[pos] == 0) {
    pos++;
  }
  digits[pos]--;
}

big_integer& big_integer::unsigned_sum(const big_integer& other) {
  adjust_size(other);
  digits.reserve(len() + 1);
  int carry = 0;
  for (std::size_t i = 0; i < other.len(); i++) {
    uint32_t temp = digits[i];
    digits[i] += carry + other.digits[i];
    carry = (digits[i] < std::max(temp, other.digits[i]) || (std::min(temp, other.digits[i]) == DIG_MAX && carry == 1));
  }
  if (carry) {
    position_increment(other.len());
  }
  return *this;
}

big_integer& big_integer::unsigned_diff(const big_integer& other, bool less) {
  digits.reserve(other.len());
  int carry = 0;
  for (std::size_t i = 0; i < (less ? len() : other.len()); i++) {
    uint32_t temp = digits[i];
    digits[i] = less ? other.digits[i] - digits[i] : digits[i] - other.digits[i];
    digits[i] -= carry;
    carry = (less ? other.digits[i] < temp : temp < other.digits[i] || (temp == other.digits[i] && carry));
  }
  if (less) {
    std::copy(other.digits.begin() + static_cast<ptrdiff_t>(len()), other.digits.end(),
              digits.begin() + static_cast<ptrdiff_t>(len()));
  }
  if (carry) {
    position_decrement(other.len());
  }
  trim_zero();
  if (is_zero()) {
    sign = false;
  }
  return *this;
}

big_integer::big_integer() = default;

big_integer::big_integer(const big_integer& other) = default;

big_integer::big_integer(int a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(long long a)
    : big_integer(static_cast<uint64_t>(a == SMIN64 ? SMAX64 : a * (a < 0 ? -1 : 1)) + (a == SMIN64)) {
  sign = (a < 0);
  trim_zero();
}

big_integer::big_integer(unsigned int a) : big_integer(static_cast<uint64_t>(a)) {}

big_integer::big_integer(unsigned long long a) {
  digits = {static_cast<uint32_t>(a % SHIFT_ONE), static_cast<uint32_t>(a >> DIGIT_LENGTH)};
  trim_zero();
}

big_integer::big_integer(long int a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(unsigned long int a) : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(const std::string& str) {
  std::size_t first = 0, last = str.size();
  if (first == last) {
    throw std::invalid_argument("empty string not a number");
  }

  if (str == "0") {
    return;
  }

  if (str[0] == '-') {
    sign = true;
    first++;
    if (first == last) {
      throw std::invalid_argument("minus isn't a number");
    }
  }

  digits = std::vector<uint32_t>(1, 0);
  for (; first < last; first += DEC_SHIFT) {
    std::size_t cnt = 0;
    int digit = std::stoi(str.substr(first, DEC_SHIFT), &cnt);
    if (cnt != std::min(last - first, DEC_SHIFT) || digit < 0) {
      throw std::invalid_argument("unknown digit");
    }

    uint32_t mul = LONG_BASE;
    for (std::size_t i = 0; i < DEC_SHIFT - cnt; i++) {
      mul /= BASE;
    }
    mul_short(mul);

    uint32_t temp = digits[0];
    digits[0] += static_cast<uint32_t>(digit);
    if (temp > DIG_MAX - static_cast<uint32_t>(digit)) {
      position_increment(1);
    }
  }
  trim_zero();
}

big_integer::~big_integer() = default;

big_integer& big_integer::operator=(const big_integer& other) {
  if (&other != this) {
    big_integer(other).swap(*this);
  }
  return *this;
}

big_integer& big_integer::signed_plus_minus(const big_integer& rhs, bool plus) {
  if (plus ? (sign == rhs.sign) : (sign != rhs.sign)) {
    return unsigned_sum(rhs);
  }
  bool less = abs_lt(rhs);
  sign = less ? plus == rhs.sign : sign;
  if (less) {
    return unsigned_diff(rhs, true);
  }
  return unsigned_diff(rhs, false);
}

big_integer& big_integer::operator+=(const big_integer& rhs) {
  return signed_plus_minus(rhs, true);
}

big_integer& big_integer::operator-=(const big_integer& rhs) {
  return signed_plus_minus(rhs, false);
}

big_integer& big_integer::operator+=(int a) {
  uint32_t digit = static_cast<uint32_t>(static_cast<int64_t>(a) * (a < 0 ? -1LL : 1LL));
  if (len() == 0) {
    digits = {digit};
    sign = (a < 0);
    return *this;
  }
  if (sign ^ (a < 0)) {
    if (len() == 1 && digits[0] < digit) {
      sign = !sign;
      digits[0] = digit - digits[0];
      return *this;
    }
    bool less = digits[0] < digit;
    digits[0] -= digit;
    if (less) {
      position_decrement(1);
    }
  } else {
    digits.reserve(len() + 1);
    bool too_much = digits[0] > DIG_MAX - digit;
    digits[0] += digit;
    if (too_much) {
      position_increment(1);
    }
  }
  return *this;
}

big_integer& big_integer::operator-=(int a) {
  if (a == SMIN32) {
    digits.reserve(len() + 1);
    *this += SMAX32;
    position_increment(0);
  } else {
    *this += -a;
  }
  return *this;
}

void big_integer::mul_short(uint32_t digit) {
  uint32_t carry = 0;
  digits.reserve(len() + 1);
  for (std::size_t i = 0; i < len(); i++) {
    uint64_t prod = digits[i];
    prod *= static_cast<uint64_t>(digit);
    prod += static_cast<uint64_t>(carry);
    digits[i] = static_cast<uint32_t>(prod % (static_cast<uint64_t>(1) << DIGIT_LENGTH));
    carry = static_cast<uint32_t>(prod >> DIGIT_LENGTH);
  }

  if (carry != 0) {
    digits.push_back(carry);
  }
}

big_integer& big_integer::operator*=(const big_integer& rhs) {
  sign ^= rhs.sign;
  std::size_t n = rhs.len(), m = len();
  digits.resize(n + m, 0);
  for (std::size_t i = 0; i < m; i++) {
    std::swap(digits[n + m - 1 - i], digits[m - 1 - i]);
  }

  for (std::size_t i = 0; i < m; i++) {
    uint32_t carry = 0;
    uint32_t digit = digits[n + i];
    digits[n + i] = 0;

    for (std::size_t j = 0; j < n; j++) {
      uint64_t prod = digit;
      prod *= static_cast<uint64_t>(rhs.digits[j]);
      prod += static_cast<uint64_t>(carry);

      auto least_bits = static_cast<uint32_t>(prod % (static_cast<uint64_t>(1) << DIGIT_LENGTH));
      auto most_bits = static_cast<uint32_t>(prod >> DIGIT_LENGTH);
      uint32_t last_val = digits[i + j];
      digits[j + i] += least_bits;

      if (digits[i + j] < std::max(last_val, least_bits)) {
        most_bits++;
      }
      carry = most_bits;
    }

    digits[i + n] = carry;
  }
  trim_zero();
  return *this;
}

template <class BinaryFunction>
big_integer& big_integer::bit_operation(const big_integer& other, BinaryFunction func) {
  adjust_size(other);
  digits.reserve(len() + 1);
  uint32_t carry1 = 1, carry2 = 1;

  for (std::size_t pos = 0; pos < len(); ++pos) {
    uint32_t d1 = digits[pos], d2 = (pos == other.len()) ? static_cast<uint32_t>(0) : other.digits[pos];
    if (sign) {
      d1 = ~d1 + carry1;
      carry1 = (d1 == 0) && carry1;
    }
    if (other.sign) {
      d2 = ~d2 + carry2;
      carry2 = (d2 == 0) && carry2;
    }
    digits[pos] = func(d1, d2);
  }

  uint32_t d1 = sign ? DIG_MAX : static_cast<uint32_t>(0);
  uint32_t d2 = other.sign ? DIG_MAX : static_cast<uint32_t>(0);
  sign = func(d1, d2);
  if (sign) {
    std::transform(digits.begin(), digits.end(), digits.begin(), std::bit_not());
    position_increment(0);
  }

  trim_zero();
  return *this;
}

big_integer& big_integer::div_mod(const big_integer& other, big_integer& remainder) {
  if (is_zero()) {
    return *this;
  }
  std::size_t n = other.len();
  if (len() < n) {
    swap(remainder);
    return *this;
  }

  int final_sign = sign ^ other.sign, initial_sign = sign;
  sign = false;

  big_integer b(other);
  b.sign = false;
  int k = 0;
  while (b.digits[n - 1] < static_cast<uint32_t>(1 << (DIGIT_LENGTH - 1 - k))) {
    ++k;
  }

  if (k != 0) {
    b <<= k;
    *this <<= k;
  }

  std::size_t m = len() - other.len();
  big_integer q;
  q.digits.resize(m + 1, 0);
  if (!abs_lt(b << static_cast<int>(m * DIGIT_LENGTH))) {
    q.digits[m] = 1;
    unsigned_diff(b << static_cast<int>(m * DIGIT_LENGTH), false);
  }

  do {
    if (m == 0) {
      break;
    }
    --m;
    auto q_star =
        (static_cast<uint64_t>(n + m < len() ? digits[n + m] : 0) * (static_cast<uint64_t>(1) << DIGIT_LENGTH) +
         static_cast<uint64_t>(n + m <= len() ? digits[n + m - 1] : 0)) /
        b.digits[n - 1];
    if (q_star > static_cast<uint64_t>(DIG_MAX)) {
      q_star = static_cast<uint64_t>(DIG_MAX);
    }
    q.digits[m] = static_cast<uint32_t>(q_star);
    big_integer temp(b);
    temp <<= static_cast<int>(m * DIGIT_LENGTH);
    big_integer diff(temp);
    diff.mul_short(q.digits[m]);
    diff.trim_zero();

    while (abs_lt(diff)) {
      --q.digits[m];
      diff -= temp;
    }

    unsigned_diff(diff, false);
  } while (m > 0);

  q.trim_zero();
  q.sign = final_sign;
  if (k > 0) {
    *this >>= k;
  }
  sign = !is_zero() && initial_sign;
  swap(remainder);
  swap(q);
  return *this;
}

big_integer& big_integer::operator/=(const big_integer& rhs) {
  big_integer remainder(0);
  return div_mod(rhs, remainder);
}

big_integer& big_integer::operator%=(const big_integer& rhs) {
  big_integer remainder(0);
  div_mod(rhs, remainder);
  swap(remainder);
  return *this;
}

big_integer& big_integer::operator&=(const big_integer& rhs) {
  return bit_operation(rhs, std::bit_and());
}

big_integer& big_integer::operator|=(const big_integer& rhs) {
  return bit_operation(rhs, std::bit_or());
}

big_integer& big_integer::operator^=(const big_integer& rhs) {
  return bit_operation(rhs, std::bit_xor());
}

big_integer& big_integer::operator<<=(int rhs) {
  if (rhs == 0) {
    return *this;
  }

  std::size_t upper = (rhs + DIGIT_LENGTH - 1) / DIGIT_LENGTH;
  std::size_t idx = len() + upper;
  std::size_t diff = DIGIT_LENGTH - rhs % DIGIT_LENGTH;
  digits.resize(idx, 0);

  do {
    idx--;
    if (diff == DIGIT_LENGTH) {
      digits[idx] = (idx >= upper) ? digits[idx - upper] : 0;
      continue;
    }
    if (idx >= upper - 1) {
      digits[idx] = (digits[idx - (upper - 1)] % (static_cast<uint32_t>(1) << diff)) << (rhs % DIGIT_LENGTH);
      if (idx >= upper) {
        digits[idx] += digits[idx - upper] >> diff;
      }
    } else {
      digits[idx] = 0;
    }
  } while (idx > 0);

  trim_zero();
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  std::size_t diff = DIGIT_LENGTH - rhs % DIGIT_LENGTH;

  for (std::size_t idx = 0; idx < len(); idx++) {
    if (idx + rhs / DIGIT_LENGTH < len()) {
      digits[idx] = digits[idx + rhs / DIGIT_LENGTH] >> (rhs % DIGIT_LENGTH);

      if (idx + rhs / DIGIT_LENGTH < len() - 1) {
        digits[idx] += (digits[idx + rhs / DIGIT_LENGTH + 1] % (static_cast<uint32_t>(1) << (rhs % DIGIT_LENGTH)))
                    << diff;
      }
    } else {
      digits[idx] = 0;
    }
  }

  trim_zero();
  if (sign) {
    --*this;
  }
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  big_integer tmp(*this);
  tmp.sign = !is_zero() && !sign;
  return tmp;
}

big_integer big_integer::operator~() const {
  big_integer tmp(*this);
  tmp.sign = !is_zero() && !sign;
  --tmp;
  return tmp;
}

big_integer& big_integer::operator++() {
  if (sign) {
    position_decrement(0);
    sign = !is_zero() && sign;
  } else {
    position_increment(0);
  }
  return *this;
}

big_integer big_integer::operator++(int) {
  big_integer temp = *this;
  ++*this;
  return temp;
}

big_integer& big_integer::operator--() {
  if (sign || is_zero()) {
    position_increment(0);
    sign = true;
  } else {
    position_decrement(0);
  }
  return *this;
}

big_integer big_integer::operator--(int) {
  big_integer temp = *this;
  --*this;
  return temp;
}

big_integer operator+(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp += b;
  return tmp;
}

big_integer operator-(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp -= b;
  return tmp;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp *= b;
  return tmp;
}

big_integer operator/(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp /= b;
  return tmp;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp %= b;
  return tmp;
}

big_integer operator&(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp &= b;
  return tmp;
}

big_integer operator|(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp |= b;
  return tmp;
}

big_integer operator^(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp ^= b;
  return tmp;
}

big_integer operator<<(const big_integer& a, int b) {
  big_integer tmp(a);
  tmp <<= b;
  return tmp;
}

big_integer operator>>(const big_integer& a, int b) {
  big_integer tmp(a);
  tmp >>= b;
  return tmp;
}

bool operator==(const big_integer& a, const big_integer& b) = default;

bool operator!=(const big_integer& a, const big_integer& b) {
  return !(a == b);
}

bool operator<(const big_integer& a, const big_integer& b) {
  if (a.sign != b.sign) {
    return a.sign;
  }
  return a.sign ? b.abs_lt(a) : a.abs_lt(b);
}

bool operator>(const big_integer& a, const big_integer& b) {
  return b < a;
}

bool operator<=(const big_integer& a, const big_integer& b) {
  return !(a > b);
}

bool operator>=(const big_integer& a, const big_integer& b) {
  return !(a < b);
}

std::string to_string(const big_integer& a) {
  big_integer tmp(a);
  bool sign = a.sign;
  tmp.sign = false;
  std::string out;
  if (tmp.is_zero()) {
    out = "0";
    return out;
  }
  while (!tmp.is_zero()) {
    std::size_t idx = tmp.len();
    uint64_t remainder = 0;
    do {
      --idx;
      remainder *= (static_cast<uint64_t>(1) << DIGIT_LENGTH);
      remainder += static_cast<uint64_t>(tmp.digits[idx]);
      tmp.digits[idx] = static_cast<uint32_t>(remainder / static_cast<uint64_t>(LONG_BASE));
      remainder %= static_cast<uint64_t>(LONG_BASE);
    } while (idx > 0);
    tmp.trim_zero();
    std::string str = std::to_string(remainder);
    std::reverse(str.begin(), str.end());
    if (!tmp.is_zero()) {
      str.resize(DEC_SHIFT, '0');
    }
    out += str;
  }
  if (sign) {
    out += "-";
  }
  std::reverse(out.begin(), out.end());
  return out;
}

std::ostream& operator<<(std::ostream& s, const big_integer& a) {
  return s << to_string(a);
}
