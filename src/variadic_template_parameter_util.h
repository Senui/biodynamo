#ifndef VARIADIC_TEMPLATE_PARAMETER_UTIL_H_
#define VARIADIC_TEMPLATE_PARAMETER_UTIL_H_

namespace bdm {

// Inspiration taken from:
// https://stackoverflow.com/questions/4691657/is-it-possible-to-store-a-template-parameter-pack-without-expanding-it

/// This solution is a workaround to pass around template parameter packs
/// without expanding them. This is necesary since creating a typedef or alias
/// fails. The following code example shows how to extract the parameter pack:
///
///     // Example: Create a std::tuple from a given VariadicTypedef
///     // Type definition
///     template <typename... Types>
///     struct ConvertToTuple {};
///
///     // Template specialization with automatic type deduction gives back
///     // the original parameter pack
///     template<typename... Types>
///     struct ConvertToTuple< VariadicTypedef<Types...> > {
///        typedef std::tuple<Types...> type;
///     }
///
///     // code that wants to pass on a parameter pack to a class
///     MyClass<VariadicTypedef<int, double>> class;
///
///     // client code that receives a VariadicTypedef type and wants to create
///     // a tuple can do so by:
///     typename ConvertToTuple<TVariadicTypedef>::type tuple;
template <typename... Types>
struct VariadicTypedef {};

/// Forward declaration of default template parameter requires wrapping of
/// VariadicTypedef
/// This method makes testing easier since the wrapper does not have to be
/// defined manually
template <typename... Types>
struct VariadicTypedefWrapper {
  typedef VariadicTypedef<Types...> types;  // NOLINT
};

/// Defines the atomic simulation objects which should be used in a certain
/// simulation.
/// CAUTION: Needs to be called inside namespace `::bdm`, since `AtomicTypes`
/// was forward declared in this namespace. Otherwise compilation will fail.
#define BDM_DEFINE_ATOMIC_TYPES(...)            \
  struct AtomicTypes {                          \
    typedef VariadicTypedef<__VA_ARGS__> types; \
  };

}  // namespace bdm

#endif  // VARIADIC_TEMPLATE_PARAMETER_UTIL_H_
