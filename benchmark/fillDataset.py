from Cryptodome.Util.number import getPrime
import json
import math

x = math.pow(256, 1/100)
numbits = [int(pow(x, i)) for i in range(116)]
numbits = [num for num in numbits if num > 8]


f = open("dataset.json", "w")
f.write('[\n')

last = 1
for i, numbit in enumerate(numbits):
    p = getPrime(numbit//2)
    q = p
    n = p*q

    while q==p or n.bit_length() < numbit or n == last:
        q = getPrime((numbit+1)//2)
        n = p*q

    #print(n)
    last = n
    line = ({'index': i, 'numbits (of n)': numbit, 'n': n, '(p, q)': (p, q)})
    json.dump(line, f)
    if i != len(numbits)-1:
        f.write(',')
    f.write('\n')

f.write(']')
f.close()