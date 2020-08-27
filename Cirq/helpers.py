import numpy as np
from typing import *

# converting between binary and int
def toBinary(num: int, numBits: int) -> List[int]:
    return list(reversed([int((num & (2**i)) > 0) for i in range(numBits)]))

# converts entries in vector to binary and concatenates results
def convertVectToBinary (x : List[int], modulus : int) -> List[int]:
    bits = int(np.log2(modulus))
    lst = []
    n = len(x)
    for i in range(n):
        lst += reversed(toBinary(x[i], bits)) # reversed because of little endian convention
    return lst

# checks if d . (x1 xor x2) = b
def checkEquation (x1 : List[int], x2: List[int], d : List[int], b : int, modulus : int) -> bool:
    x1bin = convertVectToBinary(x1, modulus)
    x2bin = convertVectToBinary(x2, modulus)
    x1xorx2 = [(a ^ b) for (a, b) in zip(x1bin, x2bin)]
    ds = np.dot(d, x1xorx2) % 2
    return (ds == b)

# converts from num in base 10 to base given
def toBase(num : int, numBits : int, base : int) -> List[int]:
    if (num == 0):
        return [0] * numBits
    digits = []
    while (num > 0):
        digits.append(num % base)
        num //= base
    if (len(digits) < numBits):
        digits += [0] * (numBits - len(digits))
    return digits[::-1]

# finds preimage x such that [Ax mod modulus] = y
def findPreimage (y: List[int], A: List[List[int]], modulus: int) -> List[int]:
    n = len(A[0])
    for i in range(0, modulus ** n):
        x = np.array(toBase(i, n, modulus))
        res = (np.dot(A, x) % modulus) >= modulus / 2
        if ((res == y).all()):
            return np.array(list(map(lambda e: [e], x)))
    return []

def roundingBinary(x : List[int], q : int) -> List[int]:
    roundX = list(map(lambda e: 0 if e < q/2 else 1, x))
    return roundX
