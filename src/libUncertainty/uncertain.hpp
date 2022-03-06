#pragma once

#include <utility>

/** @file uncertain.hpp
  * @brief The uncertain class
  * @author C.D. Clark III
  * @date 03/05/22
  */

namespace libUncertainty
{
template<typename NT, typename UT = decltype(NT() - NT())>
class uncertain
{
 public:
  using nominal_type     = NT;
  using uncertainty_type = UT;

  uncertain()                 = default;
  uncertain(const uncertain&) = default;
  uncertain(uncertain&)       = default;
  uncertain(uncertain&&)      = default;
  uncertain& operator=(const uncertain&) = default;
  uncertain& operator=(uncertain&&) = default;

  uncertain(nominal_type nom): m_storage(std::move(nom), 0) {}
  uncertain(nominal_type nom, uncertainty_type unc): m_storage(std::move(nom), std::move(unc)) {}
  uncertain(const std::array<nominal_type,2>& data): m_storage(data[0], data[1]) {}
  template<typename T, typename std::enable_if< std::is_same<decltype(std::declval<T>()[0],int()),int>::value ,int>::type = 0>
  uncertain(const T& data) : m_storage(data[0], data[1]){}



  nominal_type     nominal() const { return m_storage.first; }
  uncertainty_type uncertainty() const { return m_storage.second; }

  void nominal(nominal_type a_val) { m_storage.first = a_val; }
  void uncertainty(uncertainty_type a_val) { m_storage.second = a_val; }

  nominal_type upper() const { return m_storage.first + static_cast<nominal_type>(m_storage.second); }
  nominal_type lower() const { return m_storage.first - static_cast<nominal_type>(m_storage.second); }

  auto relative_uncertainty() const -> decltype(this->uncertainty() / this->nominal()) { return this->uncertainty() / this->nominal(); }

  friend std::ostream& operator<<(std::ostream& out, const uncertain<nominal_type, uncertainty_type>& a_val)
  {
    out << a_val.nominal() << " +/- " << a_val.uncertainty();
    return out;
  }

 private:
  std::pair<nominal_type, uncertainty_type> m_storage;
};
}  // namespace libUncertainty
