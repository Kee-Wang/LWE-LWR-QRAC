import cirq
from cirq.contrib.svg import circuit_to_svg
from cirq.contrib.qcircuit import circuit_to_latex_using_qcircuit
from helpers import *

# computes y = As + e
# A is m x n, s is an n-dimensional vectors (n x 1), e is an m-dimensional vector (m x 1)
def LWEfunEval (A : List[List[int]], s : List[int], e : List[int], modulus : int) -> List[int]:
    y = (np.dot(A, s) + e) % modulus
    return y

m = 6
n = 2
modulus = 8
bits = int(np.ceil(np.log2(modulus)))

# n = 2, m = 6
A = np.array([[4, 4], [0, 6], [6, 0], [4, 7], [0, 4], [7, 5]])
s = np.array([[1], [0]])
e = np.array([[1], [1], [0], [0], [0], [0]])

y = LWEfunEval(A, s, e, modulus)
y = (list(map(lambda w: w[0], y.tolist())))

nQubits = 1 + bits * n + m

# initialize registers
b = cirq.LineQubit(0)
x = []
for i in range(n):
    x.append([cirq.LineQubit(i) for i in range(1 + bits * i, bits * (i + 1) + 1)])

r = [cirq.LineQubit(i) for i in range(1 + bits * n, 1 + bits * n + m)]

# initialize circuit
circuit = cirq.Circuit()

# prepare superpositions over b and x
# have \sum_{b, x} |b>|x>|0>
circuit.append(cirq.H(b))
for i in range(n):
    for j in range(bits):
        circuit.append(cirq.H(x[i][j]))
        
# initialize phase qubis
for i in range(m):
    circuit.append(cirq.H(r[i]))

# prepare phase qubits that recover MSB with high probability when measured
# after each iteration of the loop, the i-th phase qubit should be in the state:
# |+_theta>, where theta = pi/4 * (<a, x> + b . y_i) and a is the i-th row of A

for i in range(m):
    for j in range(n):
        for k in range(bits):
            power = 2 ** k
            tot = A[i][j] * power
            if ((tot % 4) != 0 or ((tot % 4) == 0 and ((tot / 4) % 2) != 0)):
                crz = cirq.rz(tot * np.pi / 4.0).controlled()
                circuit.append(crz(x[j][k], r[i]))
    if ((y[i] % 4) != 0 or ((y[i] % 4) == 0 and ((y[i] / 4) % 2) != 0)):
        crz2 = cirq.rz(y[i] * np.pi / 4.0).controlled()
        circuit.append(crz2(b, r[i]))

# get an equation in s
circuit.append(cirq.H(b))
for i in range(n):
    for j in range(bits):
        circuit.append(cirq.H(x[i][j]))

# need to measure in basis {|+_{3pi/8}>, |-_{3pi/8}>}
for i in range(m):
    circuit.append(cirq.rz(-3.0 * np.pi / 8.0)(r[i]))
    circuit.append(cirq.H(r[i]))

circuit.append(cirq.measure(b, key='bMeas'))
for i in range(n):
    circuit.append(cirq.measure(*x[i], key='d{}'.format(i)))
circuit.append(cirq.measure(*r, key='roundRes'))

numRuns = 1000
simulator = cirq.Simulator()
result = simulator.run(circuit, repetitions=numRuns)
meas = result.measurements

# depth
print(len(cirq.Circuit(circuit.all_operations())))

# total gate count
print(len(list(circuit.all_operations())))

# latex file
with open('circuits/circuit_cirq_n=2.tex', 'w') as file:
    file.write('\\documentclass[draft]{beamer}\n')
    file.write('\\usepackage[size=custom,height=39,width=120,scale=0.7]{beamerposter}\n')
    file.write('\\usepackage[braket]{qcircuit}\n')
    file.write('\\begin{document}\n')
    file.write(circuit_to_latex_using_qcircuit(circuit))
    file.write('\\end{document}')

# circuit svg file
with open('circuits/circuit_cirq_n=2.svg', 'w') as file:
   file.write(circuit_to_svg(circuit))

# convert to OpenQASM
# with open('circuitQASM.qasm', 'w') as file:
#     file.write(circuit.to_qasm())

# ------------------------------------------------------------------------------
# Checking Results

successTot = 0
validRuns = 0
bCount = 0

for i in range(numRuns):
    d = []
    for j in range(n):
        key = 'd{}'.format(j)
        d += list(meas[key][i])
    b = int(meas['bMeas'][i])
    roundRes = list(meas['roundRes'][i])
    if (sum(d) > 0):
        bCount += b
        x = findPreimage(roundRes, A, modulus)
        if (len(x) > 0):
            z = (x - s) % modulus
            if (checkEquation(x, z, d, b, modulus)):
                successTot += 1
        validRuns += 1

print(successTot / validRuns)
print(bCount / validRuns)
