"""Geometry of N points on the unit sphere S^{d-1} in R^d.

Status: untrusted (exploration or tooling; no part of the proof).

Conventions
-----------
* A configuration is an (N, d) array X with unit rows.
* c_ij = <x_i, x_j>.  The objective is the *max cosine* c(X) = max_{i<j} c_ij; the kissing
  radius r = s / (1 - s) with s = sin(theta/2) = sqrt((1 - c)/2) increases exactly when c decreases.
* Contacts are the pairs attaining the max cosine (within a tolerance).
* Tangent coordinates: at x_i an orthonormal basis B_i ((d-1) x d) of the tangent space x_i^perp.
  A tangent vector field is stored as a flat array of length N (d-1).
"""
import numpy as np


def normalize(X):
    X = np.asarray(X, float)
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def r_of_c(c):
    """Kissing radius of spheres touching a unit sphere whose centers have angular distance arccos(c)."""
    s = np.sqrt((1 - np.asarray(c, float)) / 2)
    return s / (1 - s)


def c_of_r(r):
    s = np.asarray(r, float) / (1 + np.asarray(r, float))
    return 1 - 2 * s * s


def theta_of_c(c):
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def cosines(X):
    """Gram matrix with the diagonal set to -2 (so max() ignores it)."""
    C = X @ X.T
    np.fill_diagonal(C, -2.0)
    return C


def max_cos(X):
    return cosines(X).max()


def contacts(X, tol=1e-9):
    """Pairs (i<j) within tol of the max cosine; returns (pairs, cmax)."""
    C = cosines(X)
    cm = C.max()
    I, J = np.where(np.triu(C >= cm - tol, 1))
    return list(zip(I.tolist(), J.tolist())), float(cm)


def tangent_bases(X):
    """B[i] is a (d-1) x d matrix whose rows are an orthonormal basis of x_i^perp."""
    N, d = X.shape
    B = np.zeros((N, d - 1, d))
    for i, x in enumerate(X):
        # complete x to an orthonormal basis by QR of [x | I]
        Q, _ = np.linalg.qr(np.column_stack([x, np.eye(d)]))
        B[i] = Q[:, 1:d].T
    return B


def to_ambient(B, v):
    """Tangent coordinates (flat, length N(d-1)) -> ambient vectors (N, d)."""
    N, dm1, d = B.shape
    return np.einsum('ia,iad->id', v.reshape(N, dm1), B)


def to_tangent(B, W):
    return np.einsum('iad,id->ia', B, W).ravel()


def exp_move(X, B, v, h):
    """Move each point along the great circle with initial velocity h * v_i (tangent coordinates)."""
    W = h * to_ambient(B, v)
    n = np.linalg.norm(W, axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        Y = np.where(n > 1e-300, np.cos(n) * X + np.sin(n) * W / np.where(n > 0, n, 1), X)
    return normalize(Y)


def rotation_generators(d):
    """Basis of so(d): the matrices E_ab - E_ba, a < b."""
    gens = []
    for a in range(d):
        for b in range(a + 1, d):
            M = np.zeros((d, d)); M[a, b] = 1; M[b, a] = -1
            gens.append(M)
    return gens


def rotation_motions(X, B):
    """Infinitesimal rotations as tangent vector fields (rows), d(d-1)/2 x N(d-1)."""
    d = X.shape[1]
    return np.array([to_tangent(B, X @ M.T) for M in rotation_generators(d)])


def random_config(N, d, rng):
    return normalize(rng.normal(size=(N, d)))
