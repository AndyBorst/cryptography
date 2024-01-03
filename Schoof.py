from sage.all import *

# Given parameters
q = 257743850762632419871495
p = 11 * q * (q + 1) + 3
a = 425706413842211054102700238164133538302169176474
b = 203362936548826936673264444982866339953265530166

# Define the finite field
Fp = GF(p)

# Define the elliptic curve
E = EllipticCurve(Fp, [a, b])

# Function to compute the number of points using Schoof's algorithm
def schoof_algorithm(E):
    t = E.trace_of_frobenius()
    n = E.order()
    
    # Use the Chinese Remainder Theorem (CRT)
    # to solve the congruence equation for each prime divisor of n
    factors = n.factor()
    points = 1
    for prime, power in factors:
        q = prime ** power
        f = (p + 1 - t) % q
        g = 1
        for i in range(2 * int(sqrt(q))):  # Square root of q
            g = g * f % q
        points *= q - g
    return points

# Calculate the number of points on the curve using Schoof's algorithm
number_of_points = schoof_algorithm(E)
print("Number of points on the curve:", number_of_points)
