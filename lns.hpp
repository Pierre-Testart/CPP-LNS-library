/*
* Created by Pierre Testart
* Student at Politecnico di Milano
*/

#ifndef LNS_HPP
#define LNS_HPP

#include <cstdint>
#include <limits>
#include <cmath>
#include <array>
#include <type_traits>
#include <memory>

namespace lns
{
    /*
    * LNS number with:
     * - I bits for the integer part of the exponent (including sign)
     * - F bits for the fractional part of the exponent
     * - an approximation level of A for optimized addition (if A == -1, optimized addition is not used)
     * With a greater value of A, lookup tables will be lighter in memory, but less accurate.
     * The precision of tables is such that the maximum difference between 2 consecutive values is no more than 2^A.
     * There is a default value for A based on the value of F.
    */
    template<int I, int F, int A>
    class lns_t;

    // predefined types for each exponent size
    using lns16_t = lns_t<6, 10, 0>;
    using lns32_t = lns_t<8, 24, 10>;
    using lns64_t = lns_t<11, 53, -1>;

    /* Tables end depending on the number of fractional bits available. It is unnecessary to generate values after this
     * index, because they will be smaller than the resolution of the LNS type, and be represented as zero.
     * Some values are missing because they are too large for a 64-bit integer.
     * */
    constexpr uint64_t tablesEndByFractionalBits[59] = {3ull, 6ull, 15ull, 37ull, 90ull, 210ull, 483ull, 1093ull,
        2441ull, 5392ull, 11807ull, 25660ull, 55415ull, 119021ull, 254425ull, 541616ull, 1148766ull, 2428604ull,
        5119350ull, 10762987ull, 22574549ull, 47246249ull, 98686800ull, 205762207ull, 428301629ull, 890157689ull,
        1847424240ull, 3829066206ull, 7926567868ull, 16390006646ull, 33853755115ull, 69854993876ull, 144004955052ull,
        296599844703ull, 610379558556ull, 1255118855612ull, 2578957188491ull, 5295353332584ull, 10865584580636ull,
        22280925009259ull, 45661361782698ull, 93521747366565ull, 191441539061758ull, 391679179875619ull,
        800950510876076ull, 1637085254162964ull, 3344538693793819ull, 6829812641122331ull, 13941109198287050ull,
        28445132593062037ull, 58016022074678060ull, 118284416061382948ull, 241074720748169069ull, 491165801584519587ull,
        1000309371251887465ull, 2036794164632064873ull, 4145060247760972985ull, 8431912594992253323ull,
        17142916109569665846ull};

    template<int I, int F, int A = (F > 40 ? -1 : (F <= 16 ? 0 : F > 32 ? F - 12 : (5 * (F - 16)) / 4))>
    class lns_t
    {
        static_assert(I >= 1, "Error: template parameter I must be at least 1");
        static_assert(F >= 0, "Error: template parameter F must be positive");
        static_assert(I + F == 16 || I + F == 32 || I + F == 64, "Error: the sum of template parameters I and F must be either 16, 32 or 64");
        static_assert(A >= -1, "Error: template parameter A must be -1 or a positive integer");

        using type = lns_t<I, F, A>;
        using exponent_t = typename std::conditional<I + F == 64, int64_t, typename std::conditional<I + F == 32, int32_t, int16_t>::type>::type;
        static constexpr exponent_t exponentOne = (exponent_t)1 << F;
        /* when signBit is set to one, it indicates a negative value
         * only one of the flag bits can be set at the same time, except for negative infinity which has both signBit and infBit */
        static constexpr uint8_t signBit = 1 << 0;
        static constexpr uint8_t zeroBit = 1 << 1;
        static constexpr uint8_t  infBit = 1 << 2;
        static constexpr uint8_t  nanBit = 1 << 3;
        // when one of the special value bits is set, the value of the exponent becomes irrelevant
        static constexpr uint8_t specialValueMask = zeroBit | infBit | nanBit;

        // The lookup tables declared below are used for faster LNS addition.

        // a larger interval must be covered for 64-bit exponents, so more subtables are necessary for a smooth variation of precision
        static constexpr int tableCount = I + F == 16 ? 8 : (I + F == 32 ? 16 : 48);

        static exponent_t phiPlusValue(double expValue)
        {
            return static_cast<exponent_t>(std::round(std::log2(1 + expValue) * exponentOne));
        }
        static exponent_t phiMinusValue(double expValue)
        {
            return static_cast<exponent_t>(std::round(std::log2(1 - expValue) * exponentOne));
        }

        static constexpr uint64_t tablesEnd = F < 59 ? tablesEndByFractionalBits[F] : std::numeric_limits<uint64_t>::max();

        // get the optimal number of ignored bits for a subtable
        static int getIgnoredBits(bool isPhiPlus, uint64_t lowerTableLimit, uint64_t upperTableLimit)
        {
            double prevTableLimit = static_cast<double>(lowerTableLimit) / exponentOne;
            // when the function is called, A is never -1, but it is still checked to compute the following variable at compile-time
            constexpr uint64_t acceptableDifference = A == -1 ? 0 : (uint64_t)1 << A;

            // decrement number of ignored bits until the difference between adjacent values is small enough
            int maxIgnoredBits = static_cast<int>(std::log2(upperTableLimit - lowerTableLimit)) - 5;
            for(int i = maxIgnoredBits; i > 0; --i)
            {
                double multiplicator = std::exp2(i - F);
                double expValue1 = std::exp2(-(prevTableLimit + multiplicator));
                auto value1 = isPhiPlus ? phiPlusValue(expValue1) : phiMinusValue(expValue1);
                double expValue2 = std::exp2(-(prevTableLimit + 2 * multiplicator));
                auto value2 = isPhiPlus ? phiPlusValue(expValue2) : phiMinusValue(expValue2);
                uint64_t tableDifference = static_cast<uint64_t>(std::abs(value1 - value2));

                if(tableDifference <= acceptableDifference)
                    return i;
            }
            return 0;
        }

        template<bool PhiPlus>
        struct LookupTable
        {
            explicit LookupTable()
            {
                if(A == -1) // lookup table are disabled
                    return;

                limits[tableCount - 1] = tablesEnd;
                for(int i = tableCount - 2; i >= 0; --i)
                {
                    limits[i] = PhiPlus ? static_cast<uint64_t>(((i + 1.0) / tableCount) * tablesEnd) : limits[i + 1] / 2;
                }

                for(int i = 0; i < tableCount; ++i)
                {
                    ignoredBits[i] = getIgnoredBits(PhiPlus, i == 0 ? 0 : limits[i - 1], limits[i]);
                    ignoredMasks[i] = ((uint64_t)1 << ignoredBits[i]) - 1;
                }

                sizes[0] = static_cast<uint64_t>(ceil(limits[0] * std::exp2(static_cast<double>(-ignoredBits[0]))));
                totalSize = sizes[0];
                for(int i = 1; i < tableCount; ++i)
                {
                    sizes[i] = static_cast<uint64_t>(ceil((limits[i] - limits[i - 1]) * std::exp2(static_cast<double>(-ignoredBits[i]))));
                    totalSize += sizes[i];
                }

                offsets[0] = 0;
                for(int i = 1; i < tableCount; ++i)
                    offsets[i] = offsets[i - 1] + sizes[i - 1];

                values = std::make_unique<exponent_t[]>(totalSize);
                values[0] = PhiPlus ? exponentOne : std::numeric_limits<exponent_t>::lowest();

                for(int tableId = 0; tableId < tableCount; ++tableId)
                {
                    double prevTableLimit = tableId == 0 ? 0.0 : static_cast<double>(limits[tableId - 1]) / exponentOne;
                    double multiplicator = std::exp2(ignoredBits[tableId] - F);

                    for(uint64_t i = tableId == 0 ? 1 : 0; i < sizes[tableId]; ++i)
                    {
                        double expValue = std::exp2(-(prevTableLimit + i * multiplicator));
                        values[offsets[tableId] + i] = PhiPlus ? phiPlusValue(expValue) : phiMinusValue(expValue);
                    }
                }
            }
            // Linear interpolation between the table values
            inline exponent_t interpolate(int tableId, uint64_t value) const
            {
                auto index = value >> ignoredBits[tableId];
                auto indexWithOffset = index + offsets[tableId];

                auto lowerValue = values[indexWithOffset];
                if(ignoredBits[tableId] == 0) // if the table does not skip any values, no need to interpolate
                    return lowerValue;

                if(index == sizes[tableId] - 1)
                {
                    // special case: the lower value is the last value in the table, so we use the slope of the previous interval

                    if(PhiPlus)
                        return lowerValue - (((value & ignoredMasks[tableId]) * (values[indexWithOffset - 1] - lowerValue)) >> ignoredBits[tableId]);
                    else
                        return lowerValue + (((value & ignoredMasks[tableId]) * (lowerValue - values[indexWithOffset - 1])) >> ignoredBits[tableId]);
                }
                else
                {
                    // interpolate normally between the lower and upper values in the table

                    if(PhiPlus)
                        return lowerValue - (((value & ignoredMasks[tableId]) * (lowerValue - values[indexWithOffset + 1])) >> ignoredBits[tableId]);
                    else
                        return lowerValue + (((value & ignoredMasks[tableId]) * (values[indexWithOffset + 1] - lowerValue)) >> ignoredBits[tableId]);
                }
            }

            inline exponent_t getValue(uint64_t exponentDifference) const
            {
                // binary search of the subtable to read from
                int minTable = 0, maxTable = tableCount - 1;
                while(minTable < maxTable)
                {
                    int table = (minTable + maxTable) / 2;
                    if(exponentDifference >= limits[table])
                        minTable = table + 1;
                    else
                        maxTable = table;
                }
                return interpolate(minTable, minTable == 0 ? exponentDifference : exponentDifference - limits[minTable - 1]);
            }

            std::unique_ptr<exponent_t[]> values;
            uint64_t limits[tableCount], sizes[tableCount], offsets[tableCount], ignoredMasks[tableCount], totalSize;
            int ignoredBits[tableCount];
        };

        static const LookupTable<true> phiPlus;
        static const LookupTable<false> phiMinus;

    public:
        static constexpr int integerBits = I;
        static constexpr int fractionalBits = F;
        static constexpr int approximationLevel = A;

        // default constructor: initialize number to zero (no need to initialize exponent)
        inline lns_t() : m_flags(zeroBit) {}

        // constructor from a number (floating-point or integer)
        template<typename Type, typename = typename std::enable_if<std::is_arithmetic<Type>::value>::type>
        inline explicit lns_t(Type value)
        {
            // first, handle special cases
            if(value == 0)
                m_flags = zeroBit;
            else if(std::numeric_limits<Type>::has_infinity && std::isinf(value)) // detects both +inf and -inf
                m_flags = infBit | (value < 0 ? signBit : (uint8_t)0);
            else if((std::numeric_limits<Type>::has_quiet_NaN || std::numeric_limits<Type>::has_signaling_NaN) && std::isnan(value))
                m_flags = nanBit;
            else
            {
                /* if the number is not a special value, set the flags to store the sign,
                 * and set the exponent as the log of the absolute value of the number */
                m_flags = static_cast<uint8_t>(value < 0);
                m_exponent = static_cast<exponent_t>(std::log2(std::fabs(value)) * exponentOne);
            }
        }

        // conversion to a floating-point or integral type
        template<typename Type, typename = typename std::enable_if<std::is_arithmetic<Type>::value>::type>
        inline explicit operator Type() const
        {
            if(std::is_integral<Type>::value)
            {
                // cast to integer by casting to floating-point and casting again

                // select floating-point type based on integer size
                using FloatType = typename std::conditional<(sizeof(Type) > 4), long double,
                        typename std::conditional<(sizeof(Type) > 2), double, float
                        >::type>::type;
                return static_cast<Type>((FloatType)*this);
            }
            else
            {
                // cast to floating-point

                // first, check if the number is one of the special values
                if(m_flags & specialValueMask)
                {
                    if(m_flags & zeroBit)
                        return 0;
                    if(m_flags & infBit)
                        return (m_flags & signBit) ? -std::numeric_limits<Type>::infinity() : std::numeric_limits<Type>::infinity();

                    // only special value left: NaN
                    return std::numeric_limits<Type>::quiet_NaN();
                }

                // the number is not a special value
                Type expVal = std::exp2(static_cast<Type>(m_exponent) / exponentOne);
                return m_flags & signBit ? -expVal : expVal;
            }
        }

        inline void operator*=(const type& other)
        {
            // first, check if one of the operands is a special value
            if((m_flags | other.m_flags) & specialValueMask)
            {
                uint8_t flags = m_flags | other.m_flags;

                // if one of the operands is NaN, or if the operands are zero and +-inf, the result is NaN
                if((flags & nanBit) || (flags & (zeroBit | infBit)) == (zeroBit | infBit))
                {
                    m_flags = nanBit;
                }
                // if one of the operands is +-inf, the result is +-inf, also depending on the other operand's sign
                else if(flags & infBit)
                {
                    m_flags = infBit | ((m_flags ^ other.m_flags) & signBit);
                }
                // only special case left: one of the operands is zero, so the result is zero
                else
                {
                    m_flags = zeroBit;
                }
            }
            else
            {
                /* the number is not a special value: the result exponent is the sum of the exponents, and
                 * the result sign is an exclusive disjunction of the operands' signs */

                // handle integer underflow when the absolute value of the result is too small to be represented
                if(m_exponent < 0 && other.m_exponent < 0 && std::numeric_limits<exponent_t>::lowest() - m_exponent > other.m_exponent)
                {
                    m_flags = zeroBit;
                }
                else
                {
                    m_exponent += other.m_exponent;
                    m_flags ^= other.m_flags;
                }
            }
        }
        inline type operator*(const type& other) const
        {
            auto number = *this;
            number *= other;
            return number;
        }
        inline void operator/=(const type& other)
        {
            // first, check if one of the operands is a special value
            if((m_flags | other.m_flags) & specialValueMask)
            {
                /* if one of the operands is NaN, or if both operands are zero, or if both operands and +-inf,
                 * the result is NaN */
                if(((m_flags | other.m_flags) & nanBit) || (m_flags & other.m_flags & (infBit | zeroBit)))
                {
                    m_flags = nanBit;
                }
                // if the denominator is zero, the result is +-inf, depending on the numerator's sign
                else if(other.m_flags & zeroBit)
                {
                    m_flags = (m_flags & signBit) | infBit;
                }
                // if the numerator is zero, or the denominator is +-inf, the result is zero
                else if((m_flags & zeroBit) || (other.m_flags & infBit))
                {
                    m_flags = zeroBit;
                }
                /* only special case left: +-inf divided by finite value
                 * the result is +-inf, depending on the sign of both operands */
                else
                {
                    m_flags ^= (other.m_flags & signBit);
                }
            }
            else
            {
                /* the number is not a special value: the result exponent is the difference of the exponents, and
                 * the result sign is an exclusive disjunction of the operands' signs */

                // handle integer underflow when the absolute value of the result is too small to be represented
                if(m_exponent < 0 && other.m_exponent > 0 && std::numeric_limits<exponent_t>::lowest() + other.m_exponent > m_exponent)
                {
                    m_flags = zeroBit;
                }
                else
                {
                    m_exponent -= other.m_exponent;
                    m_flags ^= other.m_flags;
                }
            }
        }
        inline type operator/(const type& other) const
        {
            auto number = *this;
            number /= other;
            return number;
        }
        inline type operator-() const
        {
            // same number if it is NaN or zero, otherwise flip sign bit
            return type(m_exponent, (m_flags & (nanBit | zeroBit)) ? m_flags : (m_flags ^ signBit));
        }

        inline void operator+=(const type& other)
        {
            // first, check if one of the operands is a special value
            if((m_flags | other.m_flags) & specialValueMask)
            {
                // if the left operand is zero, the result of the sum is the right operand
                if(m_flags & zeroBit)
                {
                    m_flags = other.m_flags;
                    m_exponent = other.m_exponent;
                }
                /*
                 * if the right operand is zero, the result of the sum is the left operand,
                 * so the number can be left unchanged. For this reason, we only continue the
                 * operation if the right operand is NOT zero
                 */
                else if(!(other.m_flags & zeroBit))
                {
                    /* if the left operand is not a special value, the result is the right operand, because
                     * when x is a finite value:
                     * x + inf = inf
                     * x + -inf = -inf
                     * x + NaN = NaN */
                    if(!(m_flags & specialValueMask))
                    {
                        m_flags = other.m_flags;
                        m_exponent = other.m_exponent;
                    }
                    // if one of the operands is NaN, the result is NaN
                    else if((m_flags | other.m_flags) & nanBit)
                    {
                        m_flags = nanBit;
                    }
                    // at this point, the left operand can only be +-inf

                    // if the right operand is also +-inf, with the opposite sign, the result is NaN
                    else if((other.m_flags & specialValueMask) && ((m_flags ^ other.m_flags) & signBit))
                    {
                        m_flags = nanBit;
                    }
                    /* otherwise, the sum is between +-inf and a finite value, the result is +-inf,
                     * no need to change the number */
                }
            }
            else if(m_exponent == other.m_exponent) // same absolute value
            {
                if(m_flags == other.m_flags) // same sign
                {
                    /* adding the same number is equivalent to multiplying by 2, which is done by
                     * incrementing the exponent */
                    m_exponent += exponentOne;
                }
                else
                {
                    m_flags = zeroBit; // opposite sign: the result of the sum is zero
                }
            }
            else
            {
                if(A != -1) // use optimized addition with lookup tables
                {
                    exponent_t maxExponent, minExponent;
                    uint8_t maxFlags;
                    if(m_exponent > other.m_exponent)
                    {
                        maxExponent = m_exponent;
                        minExponent = other.m_exponent;
                        maxFlags = m_flags;
                    }
                    else
                    {
                        maxExponent = other.m_exponent;
                        minExponent = m_exponent;
                        maxFlags = other.m_flags;
                    }

                    // compute difference in a way to avoid overflows
                    uint64_t exponentDifference = maxExponent >= 0 && minExponent < 0 ?
                          static_cast<uint64_t>(maxExponent) + static_cast<uint64_t>(-(minExponent + 1)) + 1 :
                                  static_cast<uint64_t>(maxExponent - minExponent);

                    if(exponentDifference >= tablesEnd)
                    {
                        // the exponent difference is too high to be in the tables, the phi value is zero
                        m_exponent = maxExponent;
                        m_flags = maxFlags;
                    }
                    else
                    {
                        if(m_flags == other.m_flags)
                        {
                            // same sign: read from table phi+
                            m_exponent = maxExponent + phiPlus.getValue(exponentDifference);
                            // no need to change the flags, the result has the same sign
                        }
                        else
                        {
                            // different signs: read from table phi-
                            m_exponent = maxExponent + phiMinus.getValue(exponentDifference);
                            m_flags = maxFlags;
                        }
                    }
                }
                else // default slower addition
                {
                    /*
                     * To compute sum of positive x and y, use:
                     * exponent(x + y) = exponentX + log2(1 + 2^(exponentY - exponentX))
                     * sign(x + y) = 0
                     * If both x and y are negative, use the same exponent formula, but:
                     * sign(x + y) = 1
                     * If the numbers have different signs:
                     * if exponentY > exponentX, switch x and y
                     * exponent(x + y) = exponentX + log2(1 - 2^(exponentY - exponentX))
                     * sign(x + y) = sign(x) // keep the sign of the biggest number
                    */

                    constexpr double inverseExponentOne = 1 / static_cast<double>(exponentOne);
                    if(m_flags == other.m_flags) // same sign
                    {
                        m_exponent += static_cast<exponent_t>(std::log2(1 + std::exp2((other.m_exponent - m_exponent) * inverseExponentOne)) * exponentOne);
                        // no need to change the sign bit
                    }
                    else // different sign
                    {
                        exponent_t maxExponent, minExponent;
                        if(m_exponent > other.m_exponent)
                        {
                            maxExponent = m_exponent;
                            minExponent = other.m_exponent;
                        }
                        else
                        {
                            maxExponent = other.m_exponent;
                            minExponent = m_exponent;
                            m_flags = other.m_flags; // the result sign is the sign of the biggest operand
                        }
                        m_exponent = maxExponent + static_cast<exponent_t>(std::log2(1 - std::exp2((minExponent - maxExponent) * inverseExponentOne)) * exponentOne);
                    }
                }
            }
        }
        inline type operator+(const type& other) const
        {
            auto number = *this;
            number += other;
            return number;
        }
        inline void operator-=(const type& other)
        {
            *this += -other;
        }
        inline type operator-(const type& other) const
        {
            auto number = *this;
            number -= other;
            return number;
        }
        inline bool operator==(const type& other) const
        {
            /* for two numbers to be equal:
             * - the flags must be the same
             * - the numbers must not be NaN, as NaN == NaN is false by convention
             * - if they are not special values, the exponents must also be the same */
            return m_flags == other.m_flags && !((m_flags | other.m_flags) & nanBit) &&
                (m_exponent == other.m_exponent || (m_flags & specialValueMask));
        }
        inline bool operator!=(const type& other) const
        {
            return !(*this == other);
        }
        inline bool operator<(const type& other) const
        {
            // check for a comparison involving NaN, which is always false by convention
            if((m_flags | other.m_flags) & nanBit)
                return false;

            // if the signs are different, the negative number is always smaller
            if((m_flags ^ other.m_flags) & signBit)
                return m_flags & signBit;

            // check if one of the operands is a special value
            if((m_flags | other.m_flags) & specialValueMask)
            {
                // comparing the same special value: x < x is false
                if(m_flags == other.m_flags)
                    return false;

                /* comparing finite value x != 0 and special value of same sign (zero being positive):
                 * - if x is negative, it is always greater than the other operand, because the only negative
                 * special value is -inf
                 * - if x is positive, the other operand can be zero or +inf. x is always greater than zero
                 * By conclusion, the only case when x < special value is when the special value is +inf */
                if(!(m_flags & specialValueMask))
                {
                    // checking with == and not & checks that the flags are set to +inf, and not +-inf
                    return other.m_flags == infBit;
                }
                // same reasoning of the other operand is the non-special value
                if(!(other.m_flags & specialValueMask))
                {
                    return m_flags != infBit;
                }

                /* only case left: comparing zero and +inf
                 * the operation returns true iff the left operand is the one to be zero */
                return m_flags & zeroBit;
            }

            // the operands are non-special values of same sign, we compare the exponents
            return m_flags ? m_exponent > other.m_exponent : m_exponent < other.m_exponent;
        }
        inline bool operator>(const type& other) const
        {
            return other < *this;
        }
        inline bool operator<=(const type& other) const
        {
            // if one of the operands is NaN, the result of the comparison must be false
            return !((m_flags | other.m_flags) & nanBit) && !(other < *this);
        }
        inline bool operator>=(const type& other) const
        {
            // if one of the operands is NaN, the result of the comparison must be false
            return !((m_flags | other.m_flags) & nanBit) && !(*this < other);
        }

        inline type square() const
        {
            /* multiply exponent by 2 and keep special value bits
             * this makes every value positive, and preserves the special values zero, inf and NaN */
            return type(m_exponent * 2, m_flags & specialValueMask);
        }

        inline type sqrt() const
        {
            // if the number is NaN or negative, its square root is NaN
            if(m_flags & (nanBit | signBit))
                return type(0, nanBit);

            // otherwise, simply divide exponent by 2
            return type(m_exponent / 2, m_flags);
        }

        inline type inverse() const
        {
            // check the special cases where the number is 0 or +-inf
            if(m_flags & (zeroBit | infBit))
            {
                // 1 / 0 = inf
                if(m_flags & zeroBit)
                {
                    return type(0, infBit);
                }
                // 1 / +-inf = 0
                if(m_flags & infBit)
                {
                    return type(0, zeroBit);
                }
            }
            /* the result exponent is the opposite of the input exponent
             * the flags are preserved to keep the same sign, and so that 1 / NaN = NaN */
            return type(-m_exponent, m_flags);
        }
        template<typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
        inline FloatType log2() const
        {
            // check if the value is special or negative
            if(m_flags)
            {
                // log2(0) = -inf
                if(m_flags & zeroBit)
                {
                    return -std::numeric_limits<FloatType>::infinity();
                }
                // log2(+inf) = +inf
                if(m_flags == infBit)
                {
                    return std::numeric_limits<FloatType>::infinity();
                }
                // log2 of NaN or a negative value is NaN
                return std::numeric_limits<FloatType>::quiet_NaN();
            }
            // the log2 of the number is simply the exponent, converted to the floating-point type
            return static_cast<FloatType>(m_exponent) / exponentOne;
        }

        inline type abs() const
        {
            /* same number without the sign bit
             * this preserves the special values zero, +inf and NaN, and turns -inf into +inf */
            return type(m_exponent, m_flags & specialValueMask);
        }

        inline bool isZero() const
        {
            return m_flags & zeroBit;
        }

        inline bool isPositive() const // strictly positive test, excluding zero
        {
            return !(m_flags & (zeroBit | signBit | nanBit));
        }

        inline bool isNegative() const // strictly negative test, excluding zero
        {
            return m_flags & signBit;
        }

        inline bool isInf() const
        {
            return m_flags & infBit;
        }

        inline bool isPosInf() const
        {
            return m_flags == infBit;
        }

        inline bool isNegInf() const
        {
            return m_flags == (infBit | signBit);
        }

        inline bool isNan() const
        {
            return m_flags & nanBit;
        }

    private:
        // private constructor used to create a LNS number directly from the flags and exponent
        inline lns_t(exponent_t exponent, uint8_t flags) : m_exponent(exponent), m_flags(flags) {}

        exponent_t m_exponent;
        uint8_t m_flags;
    };

    template<int I, int F, int A>
    const typename lns_t<I, F, A>::template LookupTable<true> lns_t<I, F, A>::phiPlus{};
    template<int I, int F, int A>
    const typename lns_t<I, F, A>::template LookupTable<false> lns_t<I, F, A>::phiMinus{};
}

#endif
