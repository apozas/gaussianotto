from scipy.linalg import block_diag, expm, eigvals, inv, sqrtm, svd
from numpy import array, ceil, diag, dot, eye, exp, floor, hstack, kron, log, pi, shape, sort, tan, tanh, trace, transpose, vstack, zeros

def Energy(sig, Ffree, Om):
    """Computes the free energy of the bath(s)-WM system.

    :param sig: WM-bath(s) covariance matrix
    :type sig: numpy.ndarray
    :param Ffree: Bath(s) free Hamiltonian matrix
    :type Ffree: numpy.ndarray
    :param Om: WM's frequency
    :type Om: float

    :returns: numpy.float64
    """
    Ffree = block_diag(0.5 * Om * eye(2), Ffree)
    
    return 0.5 * trace(dot(Ffree, sig))


def Entropy(sig):
    """ Computes the entropy of a given covariance matrix
    
    :param sig: Covariance matrix
    :type sig: numpy.ndarray
    
    :returns: numpy.float64
    """
    N      = int(shape(sig)[0] / 2)
    Form   = kron(eye(N), [[0, 1], [-1, 0]])
    nufull = sort(abs(eigvals(1j * dot(Form, sig)))) # Symplectic eigenvalues
    nu     = zeros((N, 1))
    for m in range(N):
        nu[m] = nufull[2 * m]
    H = sum((nu + 1) * log((nu + 1) / 2) / 2 - (nu - 1) * log((nu - 1) / 2) / 2)
    return H


def FreeRing(N, freqs, alpha):
    """ Constructs the free Hamiltonian matrix for an oscillator chain of the
    form H = 1/2 sum_n (freqs(n)*(q_n^2+p_n^2)+alpha(n)*q_n*q_(n+1), in the
    mode-mode basis
    
    :param N: Number of oscillators in the chain
    :type N: int
    :param freqs: Frequency of each oscillator
    :type freqs: list of float of length N
    :param alpha: Coupling strength between each pair of neighbor oscillators
    :type alpha: list of float of length N-1

    :returns: numpy.ndarray
    """
    Fq         = 0.5 * (diag(freqs) + diag(alpha, 1) + diag(alpha, -1))
    Fq[0, N-1] = 0.5 * alpha[0]
    Fq[N-1, 0] = 0.5 * alpha[0]
    Fp         = 0.5 * diag(freqs)
    F          = block_diag(Fq, Fp)
    
    # For consistency, transform to mode-mode ordering
    Z = Ordering(N)
    return dot(Z, dot(F, transpose(Z)))


def Initialize(N, T, Ffree):
    """ Computes covariance matrix of bath at temperature T, given the local
    free Hamiltonian matrix Ffree

    NOTE: This algorithm only works assuming that the free Hamiltonian
    contains no q-p interaction terms. q-q and p-p terms are valid.

    :param N: Number of oscillators in the bath
    :type N: int
    :param T: Bath temperature
    :type T: float
    :param Ffree: Free bath Hamiltonian matrix
    :type Ffree: numpy.ndarray

    :returns: numpy.ndarray
    """
    # First transform Ffree to q+p ordering
    Z     = Ordering(N)
    Ffree = dot(dot(transpose(Z), Ffree), Z)
    
    # Computation of transformation matrix to normal-mode basis
    Fq    = Ffree[:N, :N]
    Fp    = Ffree[N:, N:]
    A     = dot(sqrtm(Fq), sqrtm(Fp))
    
    [O1, d, O2] = svd(A)
        
    O1 = transpose(O1)                                              
    O  = block_diag(O1, O2)

    # Now ready to compute symplectic transformation from global 
    # (normal) basis to local
    S = dot(dot(sqrtm(block_diag(diag(d), diag(d))), O), inv(sqrtm(Ffree)))
    S = transpose(S)    # For consistency with the manuscript's notation

    # Now construct state in normal basis
    
    normals = array(2 * d)    # Frequencies of normal modes
    
    # Thermal symplectic eigenvalues
    nu = (exp(normals / T) + 1) / (exp(normals / T) - 1)

    sigNormal = block_diag(diag(nu), diag(nu))

    # Now transform to local basis
    sigLocal  = dot(dot(S, sigNormal), transpose(S))

    # Lastly, transform back to mode-mode ordering
    return dot(dot(Z, sigLocal), transpose(Z))


def MakeInt(N, interact):
    """ Constructs the interaction part of the Hamiltonian matrix.
    Called by MakeStimeIndep
    
    :param N: Number of modes in the bath
    :type N: int
    :param interact: Set of modes the machine interacts with
    :type interact: list of int

    :returns: numpy.ndarray
    """
    # First make submatrix X
    # Row of interactions with the machines' q
    X                   = zeros(2 * N)
    X[2 * interact - 2] = 0.5
    # Add row of interactions with the machines' p
    X                   = vstack((X, zeros(2 * N)))
    
    # X is in the off-diagonal blocks of Fint    
    Fint = hstack((zeros((2, 2)), X))
    Fint = vstack((Fint, hstack((transpose(X), zeros((2 * N, 2 * N))))))
    
    return Fint


def MakeS(N, Om, strength, interact, t, delta, lambd, Ffree):
    """ Constructs the symplectic evolution matrix of an interaction.
    The output is given in an array where the matrix in position [i] is the
    time evolution matrix from the beginning until t[i]
    :param N: Number of modes in the bath
    :type N: int
    :param Om: WM's frequency
    :type Om: float
    :param strength: Bath-WM coupling strength
    :type strength: float
    :param interact: Set of modes the WM interacts with
    :type interact: list of int
    :param t: Array of time points used for numerical integration
    :type t: list
    :param delta: Ramp-up time
    :type delta: float
    :param lambd: Values of the switching function at the steps in t
    :type lambd: list of float
    :param Ffree: Bath free Hamiltonian matrix
    :type Ffree: numpy.ndarray

    :returns: numpy.ndarray
    """
    Form  = kron(eye(1 + N), [[0, 1], [-1, 0]])    # Symplectic form in q-p basis
    Ffree = block_diag(Om * eye(2) / 2, Ffree)     # Complete free Hamiltonian matrix
    Fint  = MakeInt(N, interact)                   # Interaction Hamiltonian matrix

    dt        = t[1] - t[0]
    ndeltaON  = int(ceil(delta / dt))
    ndeltaOFF = int(shape(t)[0]) - ndeltaON
    
    # Minimal array of time evolution matrices
    Sdt       = zeros(shape=(2 * (N + 1), 2 * (N + 1), ndeltaON + 1))
    
    # Note that the Hamiltonians involved are symmetric, and thus we compute
    # 2 * F instead of F+F^T
    for i in list(range(1, ndeltaON + 2)):
        Sdt[:, :, i - 1] = expm(2 * dot(Form, Ffree + strength * lambd[i - 1] * Fint) * dt)
    
    # Complete time evolution array
    S          = zeros(shape=(2 * (N + 1), 2 * (N + 1), shape(t)[0]))
    S[:, :, 0] = Sdt[:, :, 0]
    # Ramp-up
    for i in list(range(2, ndeltaON + 1)):
        S[:, :, i - 1] = dot(Sdt[:, :, i - 1], S[:, :, i - 2])
    # Plateau
    Splat = expm(2 * dot(Form, Ffree + strength * Fint) * dt)
    for i in list(range(ndeltaON + 1, ndeltaOFF + 1)):
        S[:, :, i - 1] = dot(Splat, S[:, :, i - 2])
    # Ramp-down
    for i in list(range(ndeltaOFF + 1, shape(t)[0] + 1)):
        S[:, :, i - 1] = dot(Sdt[:, :, shape(t)[0] - (i - 1) - 1], S[:, :, i - 2])
    
    return S


def MakeStimeIndep(N, Om, strength, interact, t, delta, lambd, Ffree):
    """ Constructs the symplectic evolution matrix of a bath-WM system
    during an interaction.
    :param N: Number of modes in the bath
    :type N: int
    :param Om: WM's frequency
    :type Om: float
    :param strength: Bath-WM coupling strength
    :type strength: float
    :param interact: Set of modes the WM interacts with
    :type interact: list of int
    :param t: Array of time points used for numerical integration
    :type t: list
    :param delta: Ramp-up time
    :type delta: float
    :param lambd: Values of the switching function at the steps in t
    :type lambd: list of float
    :param Ffree: Bath free Hamiltonian matrix
    :type Ffree: numpy.ndarray

    :returns: numpy.ndarray
    """
    Form  = kron(eye(1 + N), [[0, 1], [-1, 0]])     # Symplectic form in q-p basis
    Ffree = block_diag(Om * 0.5 * eye(2), Ffree)    # WM-bath free Hamiltonian matrix
    Fint  = MakeInt(N, interact)                    # Interaction Hamiltonian matrix
    
    dt        = t[1] - t[0]
    ndeltaON  = int(ceil(delta / dt))
    ndeltaOFF = shape(t)[0] - ndeltaON
    
    # Computation of S
    # Evolution during the plateau. Since F is constant and lambda as well,
    # integration is trivial.
    # Note that the Hamiltonians involved are symmetric, and thus we compute
    # 2*F instead of F+F^T
    
    Sev = expm(2 * dot(Form, (Ffree + strength * Fint)) * dt * (len(t) - 2 * ndeltaON))
    
    # Evolution during ramp-up and down. It is symmetric going outwards from the plateau
    for i in range(int(ndeltaOFF), int(shape(t)[0])):
        tempS = expm(2 * dot(Form, (Ffree + strength * lambd[i] * Fint)) * dt)
        Sev   = dot(tempS, dot(Sev, tempS))
            
    return Sev


def Ordering(N):
    """ Produces the 2N*2N transformation between the mode-mode ordering
    and the q-p ordering
    :param N: Number of modes in the system
    :type N: int

    :returns: numpy.ndarray
    """
    Z = zeros((2 * N, 2 * N))
    
    for j in range(N):
        Z[2 * j, j]         = 1
        Z[2 * j + 1, j + N] = 1
    
    return Z


def Switching(t, delta):
    """ Produces a smooth switching function that switches on over
    time delta, remains constant, and then switches off over time delta.
    
    :param t: Array of time points
    :type t: list of float
    :param delta: Ramp-up time
    :type delta: float
    
    :returns: numpy.ndarray
    """
    dt = t[1] - t[0]
    tf = t[len(t) - 1]
    
    ndeltaON  = int(ceil(delta / dt))
    ndeltaOFF = len(t) - ndeltaON
    
    xi           = [1 for i in range(0,len(t))]
    xi[0]        = 0
    xi[len(t)-1] = 0
    
    for i in range(1, ndeltaON):
        xi[i] = (1 - tanh(1 / tan(pi * t[i] / delta))) / 2
        
    for i in range(ndeltaOFF, len(t)-1):
        xi[i] = (1 - tanh(1 / tan(pi * (tf - t[i]) / delta))) / 2
        
    return xi