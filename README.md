# BiggerInts

BiggerInts is a C++ library that procedurally generates signed and unsigned integers of any number of bits.

All features are held in [BiggerInts.h](BiggerInts.h) in a namespace by the same name. Two alias templates are defined: uint_t and int_t. They take a (64-bit unsigned) argument specifying the number of bits and generate an unsigned or signed type respectively. Everything else in the file is subject to change or removal.

Any generated type functions identically to a built-in integer type *(except that casts are explicit and may be required in some places for overload resolution)*. If the number of bits is a power of 2 greater than 7, the types created are compact (i.e. `sizeof(int_t<n>) == n / 8`). Additionaly, if the number of bits is 8, 16, 32, or 64, a built-in type is returned. The signed and unsigned variants for a particular number of bits are identical (i.e. you can safely reinterpret_cast between them to avoid a *potentially expensive* copy). Signed variants use 2's complement.

Generated types have a full set of type_traits and numeric_limits specializations.

All operations on generated types are constant expressions and flagged as noexcept *(except for division/modulus)*.

BiggerInts uses a **lot** of 64-bit arithmetic and recursive templating. To maximize speed, compile in 64-bit mode with maximum function inlining.
