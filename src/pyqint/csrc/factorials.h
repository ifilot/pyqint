#pragma once

/**
 * @brief Compute n! (factorial) as a double.
 *
 * @param n  Non-negative integer
 *
 * @return factorial of n
 */
static double factorial(size_t n) {
    static const double ans[] = {
        1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000
    };
 
    if(n > 15) {
        return n * factorial(n-1);
    } else {
        return ans[n];
    }
}

/**
 * @brief Compute n!! (double factorial) as a double.
 *
 * @param n  Non-negative integer
 *
 * @return double factorial of n
 */
static double double_factorial(size_t n) {
    static const double ans[] = {
        1,1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,645120,2027025
    };
 
    if(n > 15) {
        return n * double_factorial(n-2);
    } else {
        return ans[n];
    }
}
