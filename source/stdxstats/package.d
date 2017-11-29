/**
 * This module provides simple range-based statistics functions.
 *
 * Currently Implemented:
 * $(UL
 *     $(LI mean)
 *     $(LI variance)
 *     $(LI standard deviation)
 * )
 */
module stdxstats;

import std.traits;
import std.typecons;
import std.range.primitives;

/**
 * Finds the mean (colloquially known as the average) of a range.
 *
 * For user-defined types an extra parameter `seed` is needed in order to correctly
 * seed the summation with the equivalent to `0`.
 *
 * The first overload of this function will return `T.init` if the range
 * is empty. However, the second overload will return `seed` on empty ranges.
 *
 * This function is $(BIGOH r.length).
 *
 * Note: On systems where `real` represents an 80bit floating point, this
 * function will be more slightly more accurate if `is(T == real)`. But,
 * it will also be notably slower.
 *
 * Params:
 *     r = An $(REF_ALTTEXT input range, isInputRange, std,range,primitives)
 *     seed = For user defined types. Should be equivalent to `0`.
 *
 * Returns:
 *     The mean of `r` when `r` is non-empty.
*/
T mean(T = double, R)(R r)
if (isInputRange!R &&
    isNumeric!(ElementType!R) &&
    !isInfinite!R)
{
    if (r.empty)
        return T.init;

    Unqual!T meanRes = 0;
    size_t i = 1;

    // Knuth & Welford mean calculation
    // division per element is slower, but more accurate
    for (; !r.empty; r.popFront())
    {
        const T delta = r.front - meanRes;
        meanRes += delta / i++;
    }

    return meanRes;
}

/// ditto
auto mean(R, T)(R r, T seed)
if (isInputRange!R &&
    !isNumeric!(ElementType!R) &&
    is(typeof(r.front + seed)) &&
    is(typeof(r.front / size_t(1))) &&
    !isInfinite!R)
{
    import std.algorithm.iteration : sum, reduce;

    // per item division vis-a-vis the previous overload is too
    // inaccurate for integer division, which the user defined
    // types might be representing
    static if (hasLength!R)
    {
        if (r.length == 0)
            return seed;

        return sum(r, seed) / r.length;
    }
    else
    {
        import std.typecons : tuple;

        if (r.empty)
            return seed;

        auto pair = reduce!((a, b) => tuple(a[0] + 1, a[1] + b))
            (tuple(size_t(0), seed), r);
        return pair[1] / pair[0];
    }
}

///
@safe @nogc pure nothrow unittest
{
    import std.math : approxEqual, isNaN;

    static immutable arr1 = [1, 2, 3];
    static immutable arr2 = [1.5, 2.5, 12.5];

    assert(arr1.mean.approxEqual(2));
    assert(arr2.mean.approxEqual(5.5));

    assert(arr1[0 .. 0].mean.isNaN);
}

@safe pure nothrow unittest
{
    import std.internal.test.dummyrange : ReferenceInputRange;
    import std.math : approxEqual;

    auto r1 = new ReferenceInputRange!int([1, 2, 3]);
    assert(r1.mean.approxEqual(2));

    auto r2 = new ReferenceInputRange!double([1.5, 2.5, 12.5]);
    assert(r2.mean.approxEqual(5.5));
}

// Test user defined types
@system pure unittest
{
    import std.bigint : BigInt;
    import std.internal.test.dummyrange : ReferenceInputRange;
    import std.math : approxEqual;

    auto bigint_arr = [BigInt("1"), BigInt("2"), BigInt("3"), BigInt("6")];
    auto bigint_arr2 = new ReferenceInputRange!BigInt([
        BigInt("1"), BigInt("2"), BigInt("3"), BigInt("6")
    ]);
    assert(bigint_arr.mean(BigInt(0)) == BigInt("3"));
    assert(bigint_arr2.mean(BigInt(0)) == BigInt("3"));

    BigInt[] bigint_arr3 = [];
    assert(bigint_arr3.mean(BigInt(0)) == BigInt(0));

    struct MyFancyDouble
    {
       double v;
       alias v this;
    }

    // both overloads
    auto d_arr = [MyFancyDouble(10), MyFancyDouble(15), MyFancyDouble(30)];
    assert(mean!(double)(cast(double[]) d_arr).approxEqual(18.333));
    assert(mean(d_arr, MyFancyDouble(0)).approxEqual(18.333));
}

/**
 * Calculates variance of a range of number-like elements. Finds either the
 * population variance or the sample variance based on the `population` argument.
 *
 * If the range has less than 3 elements, `T.init` will be returned.
 *
 * This function is $(BIGOH r.length).
 *
 * Note: On systems where `real` represents an 80bit floating point, this
 * function will be more slightly more accurate if `is(T == real)`. But,
 * it will also be notably slower.
 *
 * Params:
 *     r = An $(REF_ALTTEXT input range, isInputRange, std,range,primitives)
 *     of number-like elements
 *     population = If `true` gives the population variance and not the sample
 *     variance
 *     seed = For user defined types. Should be equivalent to `0`.
 * Returns:
 *     If `r` has three or more elements, the variance of `r`, as type `T`.
 *
 *     Otherwise, `T.init` is returned.
 */
T variance(R, T = double)(R r, Flag!"Population" population = No.Population)
if (isInputRange!R &&
    isNumeric!(ElementType!R) &&
    isNumeric!(T) &&
    !isInfinite!R)
{
    Unqual!T mean = 0;
    Unqual!T var = 0;

    return varianceImpl(r, mean, var, population);
}

/// ditto
T variance(R, T)(R r, T seed, Flag!"Population" population = No.Population)
if (isInputRange!R &&
    !isInfinite!R &&
    !isNumeric!T &&
    is(typeof(r.front + seed)) &&
    is(typeof(r.front / size_t(1))))
{
    Unqual!T mean = seed;
    Unqual!T var = seed;

    return varianceImpl(r, mean, var, population);
}

private T varianceImpl(R, T)(auto ref R r, ref T mean, ref T var, bool population)
{
    // giving the variance with one or two elements is incorrect
    static if (hasLength!R)
    {
        if (r.length < 3)
            return T.init;
    }
    else
    {
        if (r.empty)
            return T.init;
    }

    size_t i = 1;

    // Welford’s method for single pass variance
    foreach (e; r)
    {
        const T oldMean = mean;
        mean = mean + (e - mean) / i++;
        var = var + (e - mean) * (e - oldMean);
    }

    static if (!hasLength!R)
        if (i <= 3)
            return T.init;

    if (population)
        return var / (i - 1);

    return var / (i - 2);
}

///
@safe pure nothrow unittest
{
    import std.math : approxEqual, isNaN;
    auto arr1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    auto arr2 = [1, 10, 40, 15, 4, 5, 22];

    // sample variance
    assert(arr1.variance.approxEqual(9.1666));
    assert(arr2.variance.approxEqual(184.4761));

    // whole population variance
    assert(arr1.variance(Yes.Population).approxEqual(8.25));
    assert(arr2.variance(Yes.Population).approxEqual(158.1224));

    // ranges with less than three elements return T.init
    assert(arr1[0 .. 0].variance.isNaN);
    assert(arr1[0 .. 2].variance.isNaN);
}

@system pure unittest
{
    import std.bigint : BigInt;

    auto bigIntArr = [BigInt(1), BigInt(2), BigInt(3), BigInt(4), BigInt(5)];
    assert(bigIntArr.variance(BigInt(0)) == 7);
    assert(bigIntArr.variance(BigInt(0), Yes.Population) == 6);

    assert(bigIntArr[0 .. 0].variance(BigInt(0)) == BigInt.init);
    assert(bigIntArr[0 .. 2].variance(BigInt(0)) == BigInt.init);
}

@system pure unittest
{
    import std.bigint : BigInt;
    import std.internal.test.dummyrange : ReferenceInputRange;
    import std.math : approxEqual;

    auto r1 = new ReferenceInputRange!int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    auto r2 = new ReferenceInputRange!BigInt([
        BigInt(1), BigInt(2), BigInt(3), BigInt(4), BigInt(5)
    ]);

    assert(r1.variance.approxEqual(9.1666));
    assert(r2.variance(BigInt(0)) == 7);
}

// test nogc
@safe @nogc pure nothrow unittest
{
    import std.math : approxEqual;
    static immutable arr1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    assert(arr1.variance.approxEqual(9.1666));
}

/**
 * Finds the standard deviation from a $(REF_ALTTEXT input range, isInputRange, std,range,primitives)
 * of number like elements. Finds either the population standard deviation or the sample
 * standard deviation based on the `population` argument.
 *
 * If the range has less than 3 elements, `T.init` will be returned.
 *
 * This function is $(BIGOH r.length).
 *
 * Note: On systems where `real` represents an 80bit floating point, this
 * function will be more slightly more accurate if `is(T == real)`. But,
 * it will also be notably slower.
 *
 * Params:
 *     r = An $(REF_ALTTEXT input range, isInputRange, std,range,primitives)
 *     of number-like elements
 *     population = If `true` gives the population standard deviation and not the sample
 *     standard deviation
 *     seed = For user defined types. Should be equivalent to `0`.
 * Returns:
 *     If `r` has three or more elements, the standard deviation of `r`, as type `T`.
 *
 *     Otherwise, `T.init` is returned.
 */
T standardDeviation(R, T = double)(R r, Flag!"Population" population = No.Population)
if (isInputRange!R &&
    isNumeric!(ElementType!R) &&
    isFloatingPoint!(T) &&
    !isInfinite!R)
{
    import std.math : sqrt;
    return r.variance!(R, T)(population).sqrt;
}

/// ditto
T standardDeviation(R, T)(R r, T seed, Flag!"Population" population = No.Population)
if (isInputRange!R &&
    !isInfinite!R &&
    !isNumeric!T &&
    is(typeof(r.front + seed)) &&
    is(typeof(r.front / size_t(1))) &&
    is(typeof(r.front++)))
{
    return manualSqrt(r.variance(seed, population), seed);
}

///
@safe pure nothrow unittest
{
    import std.math : approxEqual, isNaN;

    int[] arr1 = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    assert(arr1.standardDeviation.approxEqual(2.7386));
    assert(arr1.standardDeviation(Yes.Population).approxEqual(2.5819));

    double[] arr2 = [
        38.94, 27.62, 51.15, 23.06, 4.89, 8.17, 12.19,
        27.30, 14.60, 23.16, 37.94, 34.56, 60.14
    ];
    assert(arr2.standardDeviation.approxEqual(16.4207));
    assert(arr2.standardDeviation(Yes.Population).approxEqual(15.7765));

    assert(arr1[0 .. 0].variance.isNaN);
    assert(arr1[0 .. 2].variance.isNaN);
}

// User defined types
@system pure unittest
{
    import std.bigint : BigInt;

    auto bigIntArr = [
        BigInt(1500), BigInt(2000), BigInt(3500), BigInt(4000), BigInt(5000),
        BigInt(6000), BigInt(7500), BigInt(8000), BigInt(9500)
    ];
    assert(bigIntArr.standardDeviation(BigInt(0)) == 2753);
    assert(bigIntArr.standardDeviation(BigInt(0), Yes.Population) == 2595);

    assert(bigIntArr[0 .. 0].standardDeviation(BigInt(0)) == BigInt.init);
    assert(bigIntArr[0 .. 2].standardDeviation(BigInt(0)) == BigInt.init);
}

// test nogc
@safe @nogc pure nothrow unittest
{
    import std.math : approxEqual;
    static immutable arr1 = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    assert(arr1.standardDeviation.approxEqual(2.7386));
}

// because std.math.sqrt doesn't work with user defined types
private T manualSqrt(T)(T x, T seed)
{
    static if (is(typeof(x == 0)))
        if (x == 0)
            return T.init;

    if (x == seed)
        return seed;

    // More accurate sqrt for int types
    static if (is(typeof(x << 1)) && is(typeof(x >> 1)))
    {
        T op = x;
        T res = seed;
        T one = 1 << 30;


        // "one" starts at the highest power of four <= than the argument.
        while (one > op)
            one >>= 2;

        while (one != seed)
        {
            if (op >= res + one)
            {
                op = op - (res + one);
                res = res +  2 * one;
            }
            res >>= 1;
            one >>= 2;
        }

        if (op > res)
            ++res;

        return res;
    }
    else static if (is(typeof(x - 1.0)) && is(typeof(x >= 1.0)) && is(typeof(x / 1.0)))
    {
        static T manualAbs(T x)
        {
            if (x < 0)
                return x * -1;
            return x;
        }

        T guess = seed;
        ++guess;

        while (manualAbs((guess * guess) / x  - 1.0) >= 0.00001)
            guess = ((x / guess) + guess) / 2;

        return guess;
    }
    else
    {
        static assert(0, "Cannot find the standardDeviation of " ~ T.stringof ~
            ", as it doesn't provide the necessary opBinary functionality.");
    }
}

pure unittest
{
    import std.bigint : BigInt;
    import std.math : approxEqual;

    assert(manualSqrt(9.0, 0.0).approxEqual(3));
    assert(manualSqrt(30.25, 0.0).approxEqual(5.5));
    assert(manualSqrt(3600, 0.0).approxEqual(60));
    assert(manualSqrt(BigInt(9), BigInt(0)) == 3);
    assert(manualSqrt(BigInt(250000), BigInt(0)) == 500);
}