from setup import *
from typing import Optional
from dataclasses import dataclass
from transcript import Transcript, Message1, Message2, Message3
from poly import Polynomial, Basis
from collections import Counter
from curve import Scalar

@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3

    def flatten(self):
        proof = {}
        proof["a_eval"] = self.msg_1.a_eval
        proof["W_a"] = self.msg_2.W_a
        proof["W_a_quot"] = self.msg_2.W_a_quot

        return proof


@dataclass
class Prover:
    group_order: int
    setup: Setup
    table: list

    def __init__(self, setup: Setup, table: list, group_order: int):
        self.setup = setup
        self.table = table
        self.group_order = group_order
        self.powers_of_x = setup.powers_of_x
        self.roots_of_unity = Scalar.roots_of_unity(group_order)
    def prove(self, witness) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")
        self.witness = witness

        # Round 1
        msg_1 = self.round_1(witness)
        self.beta, self.gamma = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.alpha, self.fft_cofactor = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()
        self.zeta = transcript.round_3(msg_3)

        return Proof(msg_1, msg_2, msg_3)

    # Prover sends commitment of m(X)
    def round_1(self, witness) -> Message1:
        setup = self.setup
        duplicates = dict(Counter(witness))
        self.m_values = [Scalar(duplicates.get(val, 0)) for val in self.table]
        self.t_values = [Scalar(val) for val in self.table]
        print("m_values: ", self.m_values)
        m_poly = Polynomial(self.m_values, Basis.LAGRANGE)
        self.m_poly = m_poly.ifft()
        self.m_comm_1 = setup.commit(self.m_poly)
        print("Commitment of m(X): ", self.m_comm_1)

        return Message1(self.m_comm_1)

    # Prover sends commitment of A(X), Q_A(X), B_0(X), Q_B(X), P(X)
    def round_2(self) -> Message2:
        setup = self.setup
        group_order = self.group_order
        beta = self.beta
        m_values = self.m_values
        t_values = self.t_values
        # 1. commit A(X)
        # 1.a. calculate A_i values
        self.A_values = []
        for i, t_i in enumerate(t_values):
            A_i = m_values[i]/(beta + t_i)
            self.A_values.append(A_i)
            # sanity check
            assert A_i == m_values[i]/(beta + t_i), "A: not equal"
        print("A_values: ", self.A_values)

        # 1.b. calculate A(X) from A_i values
        A_poly = Polynomial(self.A_values, Basis.LAGRANGE)
        # in coefficient form
        self.A_poly = A_poly.ifft()
        # 1.c. commit A(X)
        self.A_comm_1 = setup.commit(self.A_poly)
        print("Commitment of A(X): ", self.A_comm_1)

        # 2. commit Q_A(X)
        # 2.a. T(X) in coefficient form
        T_poly = Polynomial(t_values, Basis.LAGRANGE)
        # in coefficient form
        self.T_poly = T_poly.ifft()
        # 2.b. vanishing polynomial: X^N - 1, N = group_order - 1
        ZH_array = [Scalar(-1)] + [Scalar(0)] * (group_order - 1) + [Scalar(1)]
        # in coefficient form
        ZH_poly = Polynomial(ZH_array, Basis.MONOMIAL)
        print("self.roots_of_unity: ", self.roots_of_unity)

        # sanity check
        for i, A_i in enumerate(self.A_values):
            point = self.roots_of_unity[i]
            a_value = self.A_poly.coeff_eval(point)
            m_value = self.m_poly.coeff_eval(point)
            t_value = self.T_poly.coeff_eval(point)
            assert a_value == m_value / (beta + t_value) , "Not equal"
        # 2.c. Q_A(X) in coefficient form
        self.Q_A_poly = (self.A_poly * (self.T_poly + beta) - self.m_poly) / ZH_poly
        print("Q_A_poly value: ", self.Q_A_poly.values)
        # 2.d. commit Q_A(X)
        self.Q_A_comm_1 = setup.commit(self.Q_A_poly)
        print("Commitment of Q_A(X): ", self.Q_A_comm_1)

        # 3. commit B_0(X)
        # 3.a. calculate B_0_i values
        self.B_values = []
        f_values = self.witness
        for i, f_i in enumerate(f_values):
            B_i = 1 / (beta + f_i)
            self.B_values.append(B_i)
            # sanity check
            assert B_i == 1 / (beta + f_i), "B: not equal"
        print("B_values: ", self.B_values)
        # 3.b. calculate B_0(X) from B_0_i values, B_0(X) = (B(X) - B(0)) / X
        B_poly = Polynomial(self.B_values, Basis.LAGRANGE)
        # in coefficient form
        self.B_poly = B_poly.ifft()
        # f(X) = X, coefficient form: [0, 1]
        x_poly = Polynomial([Scalar(0), Scalar(1)], Basis.MONOMIAL)
        B_0_eval = self.B_poly.coeff_eval(Scalar(0))
        print("B_0_eval: ", B_0_eval)
        self.B_0_poly = (self.B_poly - B_0_eval) / x_poly
        # sanity check
        for i in range(group_order):
            point = self.roots_of_unity[i]
            b_value = self.B_poly.coeff_eval(point)
            b_0_value = self.B_0_poly.coeff_eval(point)
            assert b_value == self.B_values[i], "B_value and self.B_values[i]: Not equal"
            assert b_0_value == (b_value - B_0_eval) / point, "B_0: Not equal"
        # 3.c. commit B_0(X)
        self.B_0_comm_1 = setup.commit(self.B_0_poly)
        print("Commitment of B_0(X): ", self.B_0_comm_1)

        # 4. commit Q_B(X)
        # 4.a. calculate Q_B_i values
        # 4.b. calculate Q_B(X) from Q_B_i values
        # 4.c. commit Q_B(X)

        # 5. commit P(X)
        # 5.a. calculate P_i values
        # 5.b. calculate P(X) from P_i values
        # 5.c. commit P(X)

        print("self.beta, self.gamma:", self.beta, self.gamma)
        return Message2(self.powers_of_x[3])

    def round_3(self) -> Message3:
        setup = self.setup

        return Message3(self.powers_of_x[4])

    def rlc(self, term_1, term_2):
        return term_1 + term_2 * self.beta + self.gamma

    def generate_commitment(self, coeff: Polynomial, eval: Scalar):
        setup = self.setup
        zeta = self.zeta
        # Polynomial for (X - zeta)
        ZH_zeta_coeff = Polynomial([-zeta, Scalar(1)], Basis.MONOMIAL)
        quot_coeff = (coeff - eval) / ZH_zeta_coeff
        # witness for polynomial itself
        w = setup.commit(coeff)
        # witness for quotient polynomial
        w_quot = setup.commit(quot_coeff)
        return w, w_quot
