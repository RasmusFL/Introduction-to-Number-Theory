#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#define llint long long int

using namespace std;

/*
    Simple driver functions like print
*/

void print(const vector<llint>& v)
// prints the elements of a vector of llints
{
    cout << "{";
    for (int i = 0; i < v.size(); ++i) {
        if (i != v.size() - 1)
            cout << v[i] << ", ";
        else
            cout << v[i];
    }
    cout << "}";
}

//-----------------------------------------------------------------------------------------

/*
    The Euclidean algorithm, the extended Euclidean algorithm with Bezout's lemma, Euler's phi
    function, a solution to E2.10
*/

unsigned int Euclid(llint a, llint b)
// simple recursive version, see CLRS chapter 31
{
    if (a % b == 0)
        return abs(b);
    return Euclid(b, a % b);
}

vector<llint> ExtendedEuclid(llint a, llint b)
// see CLRS chapter 31, does not force the gcd to be positive
{
    if (b == 0) {
        vector<llint> s = {a, 1, 0};
        return s;
    }
    else {
        vector<llint> c = ExtendedEuclid(b, a % b);
        vector<llint> s = {c[0], c[2], c[1] - (a/b)*c[2]};
        return s;
    }
}

unsigned int IEuclid(llint a, llint b)
// iterative version
{
    llint q, r;

    while(true) {
        q = a / b;
        r = a - q*b;
        if (r == 0)
            return abs(b);
        a = b;
        b = r;
    }
}

vector<llint> IExtendedEuclid(llint a, llint b)
// iterative extended version
{
    llint x, y, u, v;
    x = 1;
    y = 0;
    u = 0;
    v = 1;
    llint q, r;

    while(b != 0) {
        q = a / b;
        r = a - q*b;
        llint x1, y1;
        x1 = x - q*u;       // calculate the coefficients from previous coefficients
        y1 = y - q*v;
        a = b;
        b = r;
        x = u;
        y = v;
        u = x1;
        v = y1;
    }
    vector<llint> s = {a, x, y};
    return s;
}

unsigned int InstructiveEuclid(llint a, llint b)
// Euclidean algorithm with instructions
{

    llint q;
    llint r;

    while(true) {
        q = a / b;
        r = a - q*b;    // note that r may be negative if a or b is negative
        cout << a << " = ";
        if (q < 0)
            cout << "(" << q << ")" << "*";
        else
            cout << q << "*";
        if (b < 0)
            cout << "(" << b << ")" << " + " << r << endl;
        else
            cout << b << " + " << r << endl;
        if (r == 0)
            return abs(b);
        a = b;
        b = r;
    }
}

unsigned llint Eulerphi1(unsigned int n)
// Euler's phi-function by computing directly all integers less than or equal to n
// that are coprime to n
{
    unsigned int phi = 0;
    for (int i = 1; i <= n; ++i)
        if (Euclid(i, n) == 1) ++phi;
    return phi;
}

unsigned int exercise2_10(double lim)   // use a positive double
{
    ofstream ofs{ "exercise2_10.txt" };     // print to a .txt-file
    double q;
    double n = 1;

    while(true) {
        q = double(Eulerphi1(n))/double(n);
        cout << "(n, phi(n)/n) = (" << n << ", " << q << ")" << endl;
        ofs << "(n, phi(n)/n) = (" << n << ", " << q << ")" << endl;
        if (q <= lim) {
            n = int(n);
            return n;
        }
        ++n;
    }
}

//-----------------------------------------------------------------------------------------

/*
    Functions for primality testing, prime factorization, Carmichael numbers, the Sieve of
    Erastothenes
*/

bool isPrime(unsigned llint a)
// check that a is prime through trial division
{
    if (a == 1)
        return false;

    int sqa = sqrt(a);      // it is sufficient to check all non-trivial divisors less than or equal to sqrt(a)
    for (int i = 2; i <= sqa; ++i) {
        if (a % i == 0)
            return false;
    }

    return true;
}

vector<llint> primeDivisors(unsigned llint a)
// returns a vector of prime factors of a (using wheel factorization)
{
    vector<llint> prime_divisors;
    for (llint d : {2, 3, 5, 7}) {
        while (a % d == 0) {
            prime_divisors.push_back(d);
            a /= d;
        }
    }
    vector<llint> increments = {2, 4, 2, 4, 6, 2, 6, 4, 2};
    int i = 0;
    for (llint d = 11; d*d <= a; d += increments[i++]) {
        while (a % d == 0) {
            prime_divisors.push_back(d);
            a /= d;
        }
        if (i == 8) i = 0;
    }
    if (a > 1)
        prime_divisors.push_back(a);
    return prime_divisors;
}

struct primepower {
// a simple struct for representing pairs of a prime and its power
    llint a;
    unsigned int n;

    primepower(llint aa, unsigned int nn) : a { aa }, n { nn } {}
};

vector<primepower> primepowerpairs(const vector<llint>& v)
// returns a vector of (a, n)s where n is the power of a
{
    vector<llint> temp = v;
    sort(temp.begin(), temp.end());

    vector<primepower> res;
    unsigned int n = 1;
    for (int i = 0; i < temp.size(); ++i) {
        if (temp[i] != temp[i + 1]) {
            primepower p = primepower(temp[i], n);
            res.push_back(p);
            n = 1;
        }
        else
            ++n;
    }
    return res;
}

unsigned llint Eulerphi2(unsigned llint n)
// Euler's phi-function using the prime factorization of n, equation (2.4)
// and Corollary 2.3.5 in the Course notes
{
    unsigned llint phi = 1;
    vector<primepower> factors = primepowerpairs(primeDivisors(n));
    for (primepower p : factors) {
        phi *= (pow(p.a, p.n) - pow(p.a, p.n - 1));
    }
    return phi;
}

bool isCarmichael(const int& n)
// returns true if n is a Carmichael number
{
    vector<llint> pd = primeDivisors(n);
    int s = pd.size();
    pd.push_back(0);
    if (s == 1) return false;   // n has to be composite by definition

    for (int i = 0; i < s; ++i) {
        if (pd[i] == pd[i+1] && i != s) return false;       // n is not square-free
        if ((n-1) % (pd[i]-1) != 0) return false;           // see Korselt's criterion
    }

    return true;
}

void listCarmichaelnumbers(const unsigned int& b)   // enter an upper bound b
// makes a file of Carmichael numbers less than or equal to b
{
    ofstream ofs {"Carmichael.txt"};

    for (int n = 2; n <= b; ++n) {
        if (isCarmichael(n)) {
            cout << "n = " << n << endl;
            ofs << "n = " << n << endl;
        }
    }
}

vector<int> Erastothenes(const unsigned int n)    // max n is approximately 2*10^6
// makes a file of primes less than or equal to n
{
    ofstream ofs {"Erastothenes.txt"};  // write primes to the file Erastothenes.txt
    bool checklist[n+1] {false};
    unsigned int sqrn = sqrt(n);
    vector<int> primes;

    for (int i = 2; i <= sqrn; ++i) {
        if (checklist[i] == false) {
            for (int j = i*i; j <= n; j += i) {     // we can start from i^2
                checklist[j] = true;
            }
        }
    }
    for (int i = 2; i < n+1; ++i) {
        if (checklist[i] == false) {
            ofs << i << ", ";
            primes.push_back(i);
        }
    }
    return primes;
}

//-----------------------------------------------------------------------------------------

/*
    Functions for solving modular linear equations and the Chinese remainder theorem,
    modular exponentiation, a simple demonstration of RSA
*/

llint modLinearSolve(llint a, llint b, unsigned int n)
// solves a*x = b (mod n) for x (if possible)
{
    vector<llint> s = ExtendedEuclid(a, n);
    if (b % s[0] == 0) {
        llint x = (s[1] * b)/s[0];
        x %= n;
        if (x < 0)
            return x + n;   // return x as smallest non-negative representative
        return x;
    }

    cout << "No solution to equation: " << a << "x" << " = " << b << " mod " << n << endl;
    return NULL;
}

bool pairwiseCoprime(const vector<llint>& a)
// returns true if the integers in a are pairwise coprime
{
    for (int i = 0; i < a.size() - 1; ++i) {
        for (int j = i + 1; j < a.size(); ++j)
            if (IEuclid(a[i], a[j]) != 1) return false;
    }
    return true;

}

llint ChRem(const vector<llint>& a, const vector<llint>& n)
// finds a solution x to the modular equations a[i] = x mod n[i] for 0<=i<n if
// the integers in n are pairwise coprime
{
    int la = a.size();
    int ln = n.size();
    if (la != ln)
        throw runtime_error("Error, input lengths are not identical.");
    if (la == 0)
        return NULL;

    llint x = a[0];
    llint q = n[0];
    if (pairwiseCoprime(n)) {
        for (int i = 1; i < la; ++i) {      // see the proof of the theorem
            x += q*modLinearSolve(q, a[i] - x, n[i]);
            q *= n[i];
        }
        return x;
    }
    throw runtime_error("Error, moduli must be pairwise co-prime.");

}

string decToBin(unsigned llint a)
// converts a positive integer in base 10 to binary (as a string in reverse!)
{
    string bin;
    while(a != 0) {
        if (a % 2 == 0) {
            a = a/2;
            bin += "0";
        }
        else {
            a = (a - 1)/2;
            bin += "1";
        }

    }
    return bin;
}

llint binToDec(string bin)
// reverse the previous function
{
    int a = 0;
    for (int i = 0; i < bin.size(); ++i) {
        if (bin[i] == '1')
            a += pow(2, i);
    }
    return a;
}

llint modularExponentiation(unsigned llint a, unsigned llint b, unsigned int n)
// calculates a^b mod n using modular exponentiation
{
    if (Euclid(a, n) == 1)
        b %= Eulerphi2(n);                  // reduce the power by phi(n) by Euler's theorem
    string bin = decToBin(b);               // get the power in binary
    vector<llint> powers;
    powers.push_back(a % n);                                // the initial power is just a itself
    for (int i = 1; i < bin.size(); ++i) {                  // compute 2^i powers of a (mod n)
        powers.push_back(powers[i - 1]*powers[i - 1]);      // do successive squarings using the previous result
        powers[i] %= n;
    }
    llint res = 1;
    for (int i = 0; i < bin.size(); ++i) {
        if (bin[i] == '1') {
            res *= powers[i];
            res %= n;
        }
    }
    return res;
}

void RSA()
// a very simplified demonstration of RSA
{
    // choose two (somewhat) large primes
    int p1 = 449;
    int p2 = 1033;
    llint n = p1*p2;
    cout << "n = " << n << endl;
    llint phi_n = (p1 - 1)*(p2 - 1);    // secret!
    llint e = phi_n / 3;                // pick "random" e
    while(IEuclid(e, phi_n) != 1)
        e += 1;
    llint d = modLinearSolve(e, 1, phi_n);  // secret!
    cout << "Public key = (" << e << ", " << n << ")" << endl;
    cout << "Private key = (" << d << ", " << n << ")" << endl << endl;

    cin.get();  // just a delay

    llint a = 1300;
    llint message = modularExponentiation(a, e, n);
    cout << "I wish to send the message of my exam: " << a << endl;
    cout << "The message becomes " << a << "^" << e << " mod " << n << ", which is: " << message << endl;
    cout << "To decrypt, we compute " << message << "^" << d << " mod " << n << ":" << endl;
    cout << modularExponentiation(message, d, n);

}

//-----------------------------------------------------------------------------------------

/*
    Quadratic residues and the Theorem of Quadratic Reciprocity
*/

vector<llint> quadraticResidues(unsigned llint a)
// directly calculate the quadratic residues modulo a
{
    vector<llint> res;
    for (llint i = 0; i < a/2 + 1; ++i) {
        res.push_back((i*i) % a);
    }
    return res;
}

int JacobiSymbol(llint a, llint b)
// computes the Jacobi symbol (a/b), where b is an odd positive integer
{
    if (b % 2 == 0 || b < 1)    // check that the Jacobi Symbol is well defined
        throw runtime_error("The second argument has to be an odd positive integer");
    a %= b;
    if (a < 0)
        a = a + b;  // ensure the first argument is positive (no need to use the first supplementary law)

    int result = 1;
    while (a != 0) {
        while (a % 2 == 0) {            // use the second supplementary law
            a /= 2;
            if (b % 8 == 3 || b % 8 == 5)
                result = -result;
        }
        swap(a,b);
        if (a % 4 == 3 && b % 4 == 3)   // use quadratic reciprocity
            result = -result;
        a %= b;
    }

    if (b == 1)
        return result;
    else
        return 0;
}

//-----------------------------------------------------------------------------------------

/*
    Continued fractions, rationals to simple continued fractions and vice verse, partial
    convergents
*/

void printConFrac(const vector<llint>& cf)
// simply prints a continued fraction with "[]"
{
    cout << "[";
    for (int i = 0; i < cf.size(); ++i) {
        if (i != cf.size() - 1) cout << cf[i] << ", ";
        else cout << cf[i] << "]";
    }
}

vector<llint> rationalToConFrac(llint a, llint b)
// convert a rational number a/b to a simple continued fraction
{
    if (b == 0)
        throw runtime_error("The denominator cannot be zero");
    double d = double(a) / double(b);   // for sign check

    vector<llint> cf;
    llint q, r;
    while(b != 0) {     // see the proof of lemma 7.14
        q = a / b;
        if (q <= 0 && d < 0)    // ensures that only a_0 can be negative
            --q;
        r = a - q*b;
        cf.push_back(q);
        a = b;
        b = r;
    }
    return cf;
}

vector<llint> partialConvergents(const vector<llint>& cf, int m)
// calculates the partial convergent (p_m, q_m) of [a_0, a_1, ..., a_m]
{
    m += 2;
    if (m >= cf.size() + 2 || m < 0)
        throw runtime_error("Error, index m is out of range");
    llint p[m + 2];
    llint q[m + 2];
    p[0] = 0; p[1] = 1;
    q[0] = 1; q[1] = 0;

    for (int i = 2; i < m + 1; ++i) {
        p[i] = cf[i - 2]*p[i - 1] + p[i - 2];
        q[i] = cf[i - 2]*q[i - 1] + q[i - 2];
    }
    vector<llint> pq = {p[m], q[m]};

    return pq;
}

vector<llint> conFracToRational(const vector<llint>& cf)
// convert a simple continued fraction to a rational number
// (a vector of numerator and denominator)
{
    return partialConvergents(cf, cf.size() - 1);
}

int main()
{
    //cout << Euclid(2567251, -432341) << endl;
    //cout << exercise2_10(0.25);
    //print(primeDivisors(2465));
    //listCarmichaelnumbers(100000);
    //Erastothenes(2000000);
    //cout << Euclid(110, 51);
    //cout << IEuclid(2567251, -432341) << endl;
    //cout << ChRem({417559561,-15156266,1518890081},{546175,471,6763});
    //cout << modLinearSolve(518517, 14185, 800981) << endl;
    //print(quadraticResidues(23));
    //cout << JacobiSymbol(681758782589, 18978917876184379);
    //cout << InstructiveEuclid(284718, -187814);
    //print(primeDivisors(8418741784788414244));
    //print(IExtendedEuclid(122, 42));
    //printConFrac(rationalToConFrac(-4174817, 9157377922));
    //cout << modularExponentiation(748714, 214819, 2003);
    //vector<llint> test = rationalToConFrac(-7174817, 91575535357);
    //printConFrac(test);
    //cout << partialConvergents(test, 15)[0] << "/" << partialConvergents(test, 15)[1];
    //print(conFracToRational(test));
    //cout << JacobiSymbol(3, 15);

    return 0;
}






