# BiggerInts

BiggerInts is a header-only C++17 library that generates signed and unsigned integer types of increasingly larger sizes.
All features are held in [BiggerInts.h](BiggerInts.h) in a namespace by the same name.

The following type aliases are provided:

* `uint_t<N>` - Represents an `N`-bit unsigned integer type.
`N` can be any power of 2 starting at 8. For 8, 16, 32, and 64 this corresponds to a built-in C++ type.
Any larger value creates a unique non-dynamic value type with all the relevant operators overloaded.
* `int_t<N>` - Same as `uint_t<N>` but is signed (2's complement).
* `bigint` - Represents an arbitrarily-sized signed integer type (2's complement).
This has infinite range, but requires dynamically-allocated resources.

All included types have a full set of `type_traits` and `numeric_limits` specializations for standard library integration.

So now for the obvious question: why use `uint_t<N>` or `int_t<N>` when `bigint` exists?
`bigint` is the powerhouse type, able to store any finite integral value.
However, this requires dynamic allocations, as well as additional logic at runtime for any operation on them (e.g. when adding `bigint` values with different number of bits stored internally).
`uint_t<N>` and `int_t<N>`, on the other hand have no dynamic resources and their sizes are known at compile-time, which allows the compiler/optimizer to do things like static loop unrolling.
Because of this, appropriate usage of `uint_t<N>` or `int_t<N>` is frequently an order of magnitude faster than `bigint`.
Additionally, operations on `uint_t<N>` and `int_t<N>` are `constexpr`, which means they can be used in compile-time contexts including `constexpr` functions/variables or `constinit` variables.

However, the operating term here is *"appropriate"* usage.
For instance, if the value you're storing would definitely fit in a 256-bit integer but you choose to go completely overboard and use an 8192-bit type, it will likely be slower than just using bigint.
Thus, if you have an upper bound for what values you expect to store, use the smallest type possible.
But if you have no such upper bound, `bigint` is the way to go.

The basis for the BiggerInts types is `std::uint64_t`.
Thus, you should definitely be compiling in 64-bit mode to guarantee best performance.
