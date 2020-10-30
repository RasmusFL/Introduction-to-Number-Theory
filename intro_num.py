import math
import random
import time

# The Euclidean algorithm, the extended Euclidean algorithm with Bezout's lemma, Euler's phi function,
# a solution to E2.10
# --------------------------------------------------------------------------------------------------------------


# a simple recursive version, see CLRS chapter 31
def euclid(a: int, b: int):
    if a % b == 0:
        return abs(b)
    return euclid(b, a % b)


# see CLRS chapter 31, does not force the gcd to be positive
def extended_euclid(a: int, b: int):      # extended Euclidean Algorithm, see Bezout's lemma
    if b == 0:
        return [a, 1, 0]
    else:
        c = extended_euclid(b, a % b)
        return [c[0], c[2], c[1] - (a // b)*c[2]]


# iterative version
def i_euclid(a: int, b: int):     # an iterative version of the Euclidean Algorithm
    while True:
        q = a // b
        r = a % b
        if r == 0:
            return abs(b)
        a = b
        b = r


# iterative extended version
def i_extended_euclid(a: int, b: int):
    x, y, u, v = 1, 0, 0, 1
    while b != 0:
        q = a // b
        r = a - q*b
        x1 = x - q*u
        y1 = y - q*v
        a = b
        b = r
        x = u
        y = v
        u = x1
        v = y1
    return [a, x, y]


# Euclidean algorithm with instructions
def instructive_euclid(a: int, b: int):
    while True:
        q = a // b
        r = a - q*b
        if q < 0 and b < 0:
            print(a, " = ", "(", q, ")*(", b, ") + ", r)
        elif q < 0 and b >= 0:
            print(a, " = ", "(", q, ")*", b, " + ", r)
        elif q >= 0 and b < 0:
            print(a, " = ", q, "*(", b, ") + ", r)
        else:
            print(a, " = ", q, "*", b, " = ", r)
        if r == 0:
            return abs(b)
        a = b
        b = r


# Euler's phi-function by computing directly all integers less than or equal to n that are coprime to n
def euler_phi1(n: int):
    phi = 0
    for i in range(1, n + 1):
        if i_euclid(i, n) == 1:
            phi += 1
    return phi


def exercise2_10(lim: float):   # only works with a positive floating-point value
    file = open("exercise2_10.txt", 'w')    # print to a .txt-file
    n = 1
    while True:
        q = euler_phi1(n)/n
        print("(n, phi(n)/n) = (", n, ", ", q, ")")
        file.write("(n, phi(n)/n) = (" + str(n) + ", " + str(q) + ")\n")
        if q <= lim:
            return n
        n += 1


# Functions for primality testing, prime factorization, Carmichael numbers, the Sieve of Erastothenes
# --------------------------------------------------------------------------------------------------------------

# check that a is prime through trial division
def is_prime(a: int):
    if a == 1:
        return False
    sqa = math.sqrt(a)
    for i in range(2, int(sqa) + 1):    # sufficient to check all non-trivial divisors less than or equal to sqrt(a)
        if a % i == 0:
            return False
    return True


# returns a list of prime factors of a (positive integer)
def prime_divisors(a: int):
    factors = []
    for d in [2, 3, 5, 7]:
        while a % d == 0:
            factors.append(d)
            a //= d
    increments = [2, 4, 2, 4, 6, 2, 6, 4, 2]
    i = 0
    d = 11
    while d * d <= a:
        while a % d == 0:
            factors.append(d)
            a //= d
        if i == 8:
            i = 0
        d += increments[i]
        i += 1
    if a > 1:
        factors.append(int(a))
    return factors


# returns a list of (a, n), where n is the power of a
def prime_power_pairs(v: list):
    temp = v
    temp.sort()
    temp.append(-1)
    n = 1
    res = []
    for i in range(0, len(v) - 1):
        if temp[i] != temp[i+1]:
            res.append([temp[i], n])
            n = 1
        else:
            n += 1
    return res


# Euler's phi-function using the prime factorization of n, equation (2.4) and Corollary 2.3.5 in the course notes
def euler_phi2(n: int):
    phi = 1
    factors = prime_power_pairs(prime_divisors(n))
    for p in factors:
        phi *= (pow(p[0], p[1]) - pow(p[0], p[1] - 1))
    return phi


# returns true if n is a Carmichael number
def is_carmichael(n: int):
    factors = prime_divisors(n)
    s = len(factors)
    factors.append(0)
    if s == 1:          # n has to be composite by definition
        return False

    for i in range(0, s):
        if factors[i] == factors[i + 1] and i != s:    # n is not square-free
            return False
        if (n - 1) % (factors[i] - 1) != 0:     # see Korselt's criterion
            return False
    return True


# makes a file of Carmichael numbers less than or equal to b
def list_carmichael_numbers(b: int):
    file = open("Carmichael.txt", 'w')
    for n in range(2, b + 1):
        if is_carmichael(n):
            print("n = " + str(n))
            file.write("n = " + str(n) + "\n")


# makes a file of primes less than or equal to n
def erastothenes(n):
    file = open("Erastothenes.txt", 'w')    # write primes to the file Erastothenes.txt
    checklist = [False] * (n + 1)
    sqrn = int(math.sqrt(n))
    primes = []

    for i in range(2, sqrn + 1):
        if not checklist[i]:
            j = i * i
            while j <= n:
                checklist[j] = True
                j += i
    for i in range(2, n + 1):
        if not checklist[i]:
            file.write(str(i) + ", ")
            primes.append(i)
    return primes


# Functions for solving modular linear equations and the Chinese remainder theorem
# --------------------------------------------------------------------------------------------------------------


# solves a*x = b (mod n) if possible
def mod_linear_solve(a: int, b: int, n: int):
    s = extended_euclid(a, n)
    if b % s[0] == 0:
        x = (s[1] * b)/s[0]
        x %= n
        return int(x)
    print("No solution to equation: ", a, "x", " = ", b, " mod ", n)


# returns true if the integers in a are pairwise coprime
def pairwise_coprime(a: list):
    la = len(a)
    for i in range(0, la - 1):
        for j in range(i + 1, la):
            if i_euclid(a[i], a[j]) != 1:
                return False
    return True


# finds a solution x to the modular equations a[i] = x mod n[i] for 0<=i<n if the integers in n are pairwise coprime
def chinese_remainder_theorem(a: list, n: list):
    la = len(a)
    ln = len(n)
    if la != ln:
        raise Exception("Error, input lengths are not identical.")

    x = a[0]
    q = n[0]
    if pairwise_coprime(n):
        for i in range(1, la):
            x += q * mod_linear_solve(q, a[i] - x, n[i])
            q *= n[i]
        return x
    raise Exception("Error, moduli must be pairwise co-prime")


# converts a positive integer in base 10 to binary (as a list in reverse!)
def dec_to_bin(a: int):
    binary = []
    while a != 0:
        if a % 2 == 0:
            a //= 2
            binary.append(0)
        else:
            a = (a - 1) // 2
            binary.append(1)
    return binary


# reverse the previous function
def bin_to_dec(binary: list):
    a = 0
    for i in range(0, len(binary)):
        if binary[i] == 1:
            a += pow(2, i)
    return a


# calculates a^b mod n using modular exponentiation
def modular_exponentiation(a: int, b: int, n: int):
    if euclid(a, n) == 1:
        b %= euler_phi2(n)      # reduce the power by phi(n) by Euler's theorem
    binary = dec_to_bin(b)      # get the power in binary
    powers = [a % n]

    for i in range(1, len(binary)):
        powers.append(powers[i-1]*powers[i-1])
        powers[i] %= n

    res = 1
    for i in range(0, len(binary)):
        if binary[i] == 1:
            res *= powers[i]
            res %= n
    return res


# a very simplified demonstration of RSA
def rsa():
    p1 = 449
    p2 = 1033
    n = p1 * p2
    print("n = ", n)

    phi_n = (p1 - 1)*(p2 - 1)       # secret!
    e = random.randint(2, phi_n)    # pick a random e coprime to phi_n
    while(i_euclid(e, phi_n) != 1):
        e = random.randint(2, phi_n)
    d = mod_linear_solve(e, 1, phi_n)   # secret!
    print("Public key = (", e, ", ", n, ")")
    print("Private key = (", d, ", ", n, ")")

    time.sleep(5)

    a = 1300
    message = modular_exponentiation(a, e, n)
    print("I wish to send the message of my exam: ", a)
    print("The message becomes ", a, "^", e, " mod ", n, ", which is ", message)
    print("To decrypt, we compute ", message, "^", d, " mod ", n, ": ", modular_exponentiation(message, d, n))

# Quadratic residues and the Theorem of Quadratic Reciprocity
# --------------------------------------------------------------------------------------------------------------


# directly calculate the quadratic residues modulo a
def quadratic_residues(a):
    res = []
    for i in range(0, a // 2 + 1):
        res.append((i * i) % a)
    return res


# computes the Jacobi symbol (a/b), where b is an odd positive integer
def jacobi_symbol(a, b):
    if b % 2 == 0 or b < 1:     # check that the Jacobi symbol is well defined
        raise Exception("The second argument has to be an odd positive integer")
    a %= b
    res = 1
    while a != 0:
        while a % 2 == 0:   # use the second supplementary law
            a //= 2
            if b % 8 == 3 or b % 8 == 5:
                res = -res
        a, b = b, a
        if a % 4 == 3 and b % 4 == 3:   # use quadratic reciprocity
            res = -res
        a %= b

    if b == 1:
        return res
    else:
        return 0


# Continued fractions, rationals to simple continued fractions and vice verse, partial convergents
# --------------------------------------------------------------------------------------------------------------


# convert a rational number a/b to a simple continued fraction
def rational_to_con_frac(a: int, b: int):
    if b == 0:
        raise Exception("The denominator cannot be zero")

    d = a / b
    cf = []
    while b != 0:       # see the proof of lemma 7.14
        q = abs(a) // abs(b)
        if a / b < 0:
            q = -q
        if q <= 0 and d < 0:   # ensures that only a_0 can be negative
            q -= 1
        r = a - q * b
        cf.append(q)
        a = b
        b = r

    return cf


# calculates the partial convergent (p_m, q_m) of [a_0, a_1, ..., a_m]
def partial_convergents(cf: list, m: int):
    m += 2
    if m >= len(cf) + 2 or m < 0:
        raise Exception("Error, index m is out of range")

    p = [0, 1]
    q = [1, 0]
    for i in range(2, m + 1):
        p.append(cf[i - 2]*p[i - 1] + p[i - 2])
        q.append(cf[i - 2]*q[i - 1] + q[i - 2])
    return [p[m], q[m]]


# convert a simple continued fraction to a rational number (a list of numerator and denominator)
def con_frac_to_rational(cf: list):
    return partial_convergents(cf, len(cf) - 1)

# Testing
# --------------------------------------------------------------------------------------------------------------


# print(euler_phi1(2003))
# exercise2_10(0.25)
# print(is_prime(2003))
# print(prime_divisors(318871))
# print(is_carmichael(561))
# list_carmichael_numbers(100000)
# erastothenes(2000000)
# print(pairwise_coprime([16, 101, 14715, 15131]))
# print(mod_linear_solve(518517, 14185, 800981))
# print(chinese_remainder_theorem([417559561, -15156266, 1518890081], [546175, 471, 6763]))
# print(quadratic_residues(23))
# print(jacobi_symbol(681758782589, 18978917876184379))
# print(prime_divisors(66314079629875))
# print(prime_power_pairs(prime_divisors(1305355329341905924810625000)))
# print(euler_phi2(149787591))
# print(bin_to_dec(dec_to_bin(15195185901885817589)))
# print(modular_exponentiation(9174917984,4716,5761559))
# rsa()
# print(rational_to_con_frac(817491589891, 857857815816710099))
# print(partial_convergents(rational_to_con_frac(-7174817, 91575535357), 15))
# print(con_frac_to_rational(rational_to_con_frac(-74617, 19857611)))
