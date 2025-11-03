complex<double> operator+(const int &lhs, const complex<double> &rhs)
{
    // Convert int to complex<double>
    complex<double> lhs_complex = complex<double>(lhs, 0);

    return lhs_complex + rhs;
};
complex<double> operator+(const complex<double> &lhs, const int &rhs)
{
    // Convert int to complex<double>
    complex<double> rhs_complex = complex<double>(rhs, 0);

    return lhs + rhs_complex;
};
complex<double> operator-(const int &lhs, const complex<double> &rhs)
{
    // Convert int to complex<double>
    complex<double> lhs_complex = complex<double>(lhs, 0);

    return lhs_complex - rhs;
};
complex<double> operator-(const complex<double> &lhs, const int &rhs)
{
    // Convert int to complex<double>
    complex<double> rhs_complex = complex<double>(rhs, 0);

    return lhs - rhs_complex;
};
complex<double> operator*(const int &lhs, const complex<double> &rhs)
{
    // Convert int to complex<double>
    complex<double> lhs_complex = complex<double>(lhs, 0);

    return lhs_complex * rhs;
};
complex<double> operator*(const complex<double> &lhs, const int &rhs)
{
    // Convert int to complex<double>
    complex<double> rhs_complex = complex<double>(rhs, 0);

    return lhs * rhs_complex;
};
complex<double> operator/(const complex<double> &lhs, const int &rhs)
{
    // Convert int to complex<double>
    complex<double> rhs_complex = complex<double>(rhs, 0);

    return lhs / rhs_complex;
};
complex<double> operator/(const int &lhs, const complex<double> &rhs)
{
    // Convert int to complex<double>
    complex<double> lhs_complex = complex<double>(lhs, 0);

    return lhs_complex / rhs;
};