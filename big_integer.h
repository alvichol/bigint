#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

struct big_integer {
  big_integer();
  big_integer(const big_integer& other);
  big_integer(int a);
  big_integer(long long a);
  big_integer(unsigned int a);
  big_integer(unsigned long long a);
  big_integer(long int a);
  big_integer(unsigned long int a);
  big_integer(const std::string& str);
  ~big_integer();

  big_integer& operator=(const big_integer& other);

  big_integer& operator+=(const big_integer& rhs);
  big_integer& operator+=(int a);
  big_integer& operator-=(const big_integer& rhs);
  big_integer& operator-=(int a);
  big_integer& operator*=(const big_integer& rhs);
  big_integer& operator/=(const big_integer& rhs);
  big_integer& operator%=(const big_integer& rhs);

  big_integer& operator&=(const big_integer& rhs);
  big_integer& operator|=(const big_integer& rhs);
  big_integer& operator^=(const big_integer& rhs);

  big_integer& operator<<=(int rhs);
  big_integer& operator>>=(int rhs);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend bool operator==(const big_integer& a, const big_integer& b);
  friend bool operator!=(const big_integer& a, const big_integer& b);
  friend bool operator<(const big_integer& a, const big_integer& b);
  friend bool operator>(const big_integer& a, const big_integer& b);
  friend bool operator<=(const big_integer& a, const big_integer& b);
  friend bool operator>=(const big_integer& a, const big_integer& b);

  friend std::string to_string(const big_integer& a);

  void swap(big_integer& other) noexcept;
  std::size_t len() const noexcept;
  bool is_zero() const noexcept;
  bool abs_lt(const big_integer& other) const;
  void mul_short(uint32_t digit);

private:
  void adjust_size(const big_integer& other);
  void position_increment(std::size_t pos);
  void position_decrement(std::size_t pos);
  big_integer& unsigned_sum(const big_integer& other);
  big_integer& unsigned_diff(const big_integer& other, bool less);
  big_integer& signed_plus_minus(const big_integer& rhs, bool plus);
  void trim_zero() noexcept;
  big_integer& div_mod(const big_integer& other, big_integer& remainder);

  template <class BinaryFunction>
  big_integer& bit_operation(const big_integer& other, BinaryFunction func);

private:
  std::vector<uint32_t> digits;
  bool sign{false};
};

big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator-(const big_integer& a, const big_integer& b);
big_integer operator*(const big_integer& a, const big_integer& b);
big_integer operator/(const big_integer& a, const big_integer& b);
big_integer operator%(const big_integer& a, const big_integer& b);

big_integer operator&(const big_integer& a, const big_integer& b);
big_integer operator|(const big_integer& a, const big_integer& b);
big_integer operator^(const big_integer& a, const big_integer& b);

big_integer operator<<(const big_integer& a, int b);
big_integer operator>>(const big_integer& a, int b);

bool operator==(const big_integer& a, const big_integer& b);
bool operator!=(const big_integer& a, const big_integer& b);
bool operator<(const big_integer& a, const big_integer& b);
bool operator>(const big_integer& a, const big_integer& b);
bool operator<=(const big_integer& a, const big_integer& b);
bool operator>=(const big_integer& a, const big_integer& b);

std::string to_string(const big_integer& a);
std::ostream& operator<<(std::ostream& s, const big_integer& a);
