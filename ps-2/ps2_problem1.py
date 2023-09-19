# Graduate Computational Physics
# Problem Set 2
# problem 1
# By Hayden Orth

import numpy as np

# function to get the binary representation of a number
def get_bits(number):
    # Takes a number and returns a list of its representation in bits from
    # highest to lowest significance
    
    bytes = number.tobytes()
    bits = []
    for byte in bytes:
        bits = bits + np.flip(np.unpackbits(np.uint8(byte)), np.uint8(0)).tolist()
    return list(reversed(bits))

# get 32 bit floating point binary representation for desired value
value = 100.98763
bitlist=get_bits(np.float32(value))

# pull desired values from binary string. sign, exponent, and mantissa
print('32 bit floating point:')
print(bitlist)
print('sign bit:')
print(bitlist[0])
print('exponent:')
print(bitlist[1:9])
print('mantissa:')
print(bitlist[9:32])

# find the difference between the way Python natively represents the value versus numpy's
# 32-bit floating point
print('difference:')
print(value - np.float32(value))
