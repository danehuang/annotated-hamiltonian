from collections import OrderedDict
import numpy as np

from qutip import Qobj
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector, DensityMatrix
import qiskit.visualization 


# -----------------------------------------------------------------------------
# Qubits
# -----------------------------------------------------------------------------

zero = Statevector(np.array([1.0, 0.0]))
one = Statevector(np.array([0.0, 1.0]))


# -----------------------------------------------------------------------------
# Utility
# -----------------------------------------------------------------------------

def mask_small_values(qobj, threshold=1e-5):
    # Convert the Qobj to a numpy array
    qarray = qobj.full()
    
    # Apply the mask
    qarray[np.abs(qarray) < threshold] = 0
    
    # Return a new Qobj with small values masked
    return Qobj(qarray)


def invert_keys(D):  # maps g-->0, r -->1  to stick with computer science convention
    # Lambda function to invert the '0's and '1's
    invert_binary = lambda s: s.replace('0', '2').replace('1', '0').replace('2', '1')
    D_inverted = [OrderedDict((invert_binary(k), v) for k, v in d.items()) for d in D]
    return D_inverted


# -----------------------------------------------------------------------------
# Display
# -----------------------------------------------------------------------------

def pretty(qobj):
    if isinstance(qobj, Qobj):
        return Statevector(qobj.full()).draw("latex")
    elif isinstance(qobj, np.ndarray):
        if len(qobj.shape) == 1:
            return Statevector(qobj).draw("latex")
        elif len(qobj.shape) == 2:
            return DensityMatrix(qobj).draw("latex")
    elif isinstance(qobj, Statevector):
        return qobj.draw("latex")
    else:
        raise ValueError(f"Type {type(qobj)} not expected ...")


def histogram_final_state(final_state):
    if isinstance(final_state, Statevector):
        pass
    elif isinstance(final_state, Qobj):
        final_state = Statevector(final_state.full())
    elif isinstance(final_state, np.ndarray):
        final_state = Statevector(final_state)
    else:
        raise ValueError(f"Unexpected type {type(final_state)}")
    n = int(np.log2(len(final_state)))
    qc_extract = QuantumCircuit(n, n)
    qc_extract.initialize(final_state)
    qc_extract.measure(range(n), range(n))
    results = AerSimulator().run(qc_extract, shots=2048).result()
    answer = results.get_counts()
    return answer


def simulate(qc, shots=2048):
    results = AerSimulator().run(qc, shots=shots).result()
    return results.get_counts()


def plot_histogram(counts):
    return qiskit.visualization.plot_histogram(counts)


# -----------------------------------------------------------------------------
# Hamiltonians
# -----------------------------------------------------------------------------

def is_hamiltonian(H: np.ndarray) -> bool:
    return np.allclose(H, np.conjugate(H).T)


def linear_interpolation(H_0, H_T, T):
    def s(t):
        return t / T
    H_t = [[H_0, lambda t: 1 - s(t)], [H_T, s]]
    return H_t


# -----------------------------------------------
# Rydberg Hamiltonian
# -----------------------------------------------

def mk_H_couple(n, Omega, phi):
    sigma_x = Omega / 2 * (np.exp(1j * phi) * np.outer(zero, one) + np.exp(-1j * phi) * np.outer(one, zero))
    H = np.zeros((2**n, 2**n), dtype=np.complex128)
    for j in range(n):
        if j == 0:
            H += np.kron(sigma_x, np.eye(2**(n-1)))
        elif j == n - 1:
            H += np.kron(np.eye(2**(n-1)), sigma_x)
        else:
            H += np.kron(np.eye(2**j), np.kron(sigma_x, np.eye(2**(n-j-1))))
    return H


def nhats(n: int, j: int) -> np.ndarray:
    oo_op = np.outer(one, one)
    if j == 0:
        return np.kron(oo_op, np.eye(2**(n-1)))
    elif j == n - 1:
        return np.kron(np.eye(2**(n-1)), oo_op)
    else:
        return np.kron(np.eye(2**j), np.kron(oo_op, np.eye(2**(n-j-1))))


def mk_H_detune(n: int, Delta: float) -> np.ndarray:
    H = np.zeros((2**n, 2**n), dtype=np.complex128)
    for j in range(n):
        H += nhats(n, j)
    return Delta * H


def mk_H_interaction(n, register):
    V_jk = register.rydberg_interaction()
    H_interaction = np.zeros((2**n, 2**n), dtype=np.complex128)
    for row in range(0, n):
        for col in range(0, row+1):
            if row != col:
                print(row, col)
                H_interaction += V_jk[row, col] * nhats(n, col) @ nhats(n, row)
    return H_interaction


def mk_H_Rydberg(n, Omega, phi, Delta, register):
    H_couple = mk_H_couple(n, Omega, phi)
    H_detune = mk_H_detune(n, Delta)
    H_interaction = mk_H_interaction(n, register)
    return H_couple - H_detune + H_interaction
