from Cryptodome.Util.number import getPrime
import sys
import subprocess

nbits = sys.argv[1]
nbits = int(nbits, 10)

p = getPrime(nbits)
q = p

while q==p:
    q = getPrime(nbits)

n = p*q

print("n: ", n)
print("(p,q): ", (p,q))

print("Result:")
subprocess.run(["./factoring_algorithms", "2", str(n)])