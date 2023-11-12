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
        # msg_1
        proof["m_comm_1"] = self.msg_1.m_comm_1
        # msg_2
        proof["A_comm_1"] = self.msg_2.A_comm_1
        proof["Q_A_comm_1"] = self.msg_2.Q_A_comm_1
        proof["B_0_comm_1"] = self.msg_2.B_0_comm_1
        proof["Q_B_comm_1"] = self.msg_2.Q_B_comm_1
        proof["P_comm_1"] = self.msg_2.P_comm_1
        # msg_3
        proof["b_0_gamma"] = self.msg_3.b_0_gamma
        proof["f_gamma"] = self.msg_3.f_gamma
        proof["a_0"] = self.msg_3.a_0
        proof["h_comm_1"] = self.msg_3.h_comm_1
        proof["a_0_comm_1"] = self.msg_3.a_0_comm_1

        return proof

@dataclass
class Prover:
    group_order_N: int
    group_order_n: int
    setup: Setup
    table: list

    def __init__(self, setup: Setup, table: list, group_order_N: int, group_order_n: int):
        self.setup = setup
        self.table = table
        self.group_order_N = group_order_N
        self.group_order_n = group_order_n
        self.roots_of_unity_N = Scalar.roots_of_unity(group_order_N)
        self.roots_of_unity_n = Scalar.roots_of_unity(group_order_n)
        self.powers_of_x = setup.powers_of_x
    def prove(self, witness) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")
        self.lookup_table = [Scalar(val) for val in witness]

        # Round 1
        msg_1 = self.round_1(witness)
        self.beta = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.gamma, self.eta = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()

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
        group_order_N = self.group_order_N
        group_order_n = self.group_order_n
        beta = self.beta
        m_values = self.m_values
        t_values = self.t_values
        # 1. commit A(X)
        # 1.a. compute A_i values
        self.A_values = []
        for i, t_i in enumerate(t_values):
            A_i = m_values[i]/(beta + t_i)
            self.A_values.append(A_i)
            # sanity check
            assert A_i == m_values[i]/(beta + t_i), "A: not equal"
        print("A_values: ", self.A_values)

        # 1.b. compute A(X) from A_i values
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
        # 2.b. vanishing polynomial: X^N - 1, N = group_order_N - 1
        ZH_array = [Scalar(-1)] + [Scalar(0)] * (group_order_N - 1) + [Scalar(1)]
        # in coefficient form
        ZH_poly = Polynomial(ZH_array, Basis.MONOMIAL)
        print("self.roots_of_unity_N: ", self.roots_of_unity_N)

        # sanity check
        for i, A_i in enumerate(self.A_values):
            point = self.roots_of_unity_N[i]
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
        # 3.a. compute B_0_i values
        self.B_values = []
        f_values = self.lookup_table
        for i, f_i in enumerate(f_values):
            B_i = 1 / (beta + f_i)
            self.B_values.append(B_i)
            # sanity check
            assert B_i == 1 / (beta + f_i), "B: not equal"
        print("B_values: ", self.B_values)
        # 3.b. compute B_0(X) from B_0_i values, B_0(X) = (B(X) - B(0)) / X
        B_poly = Polynomial(self.B_values, Basis.LAGRANGE)
        # in coefficient form
        self.B_poly = B_poly.ifft()
        # f(X) = X, coefficient form: [0, 1]
        self.x_poly = Polynomial([Scalar(0), Scalar(1)], Basis.MONOMIAL)
        B_0_eval = self.B_poly.coeff_eval(Scalar(0))
        print("B_0_eval: ", B_0_eval)
        self.B_0_poly = (self.B_poly - B_0_eval) / self.x_poly
        # sanity check
        for i in range(group_order_n):
            point = self.roots_of_unity_N[i]
            b_value = self.B_poly.coeff_eval(point)
            b_0_value = self.B_0_poly.coeff_eval(point)
            assert b_value == self.B_values[i], "B_value and self.B_values[i]: Not equal"
            assert b_0_value == (b_value - B_0_eval) / point, "B_0: Not equal"
        # 3.c. commit B_0(X)
        self.B_0_comm_1 = setup.commit(self.B_0_poly)
        print("Commitment of B_0(X): ", self.B_0_comm_1)

        # 4. commit Q_B(X)
        # 4.a. f(X) in coefficient form
        f_poly = Polynomial(f_values, Basis.LAGRANGE)
        # in coefficient form
        self.f_poly = f_poly.ifft()

        # sanity check
        for i, B_i in enumerate(self.B_values):
            point = self.roots_of_unity_N[i]
            b_value = self.B_poly.coeff_eval(point)
            f_value = self.f_poly.coeff_eval(point)
            assert b_value == 1 / (beta + f_value) , "B quotient: Not equal"
        # 4.b. Q_B(X) in coefficient form
        self.Q_B_poly = (self.B_poly * (self.f_poly + beta) - Scalar(1)) / ZH_poly
        print("Q_A_poly value: ", self.Q_A_poly.values)
        # 4.c. commit Q_A(X)
        self.Q_B_comm_1 = setup.commit(self.Q_B_poly)
        print("Commitment of Q_B(X): ", self.Q_B_comm_1)

        # 5. commit P(X)
        # N - 1 - (n - 2)
        # TODO: N should >> n
        x_exponent_order = group_order_N - 1 - (group_order_n - 2)
        x_exponent_values_in_coeff = [Scalar(0)] * (x_exponent_order) + [Scalar(1)]
        x_exponent_poly = Polynomial(x_exponent_values_in_coeff, Basis.MONOMIAL)
        self.P_poly = self.B_0_poly * x_exponent_poly
        # 5.c. commit P(X)
        self.P_comm_1 = setup.commit(self.P_poly)
        print("Commitment of P(X): ", self.P_comm_1)

        return Message2(
            self.A_comm_1,
            self.Q_A_comm_1,
            self.B_0_comm_1,
            self.Q_B_comm_1,
            self.P_comm_1
        )

    def round_3(self) -> Message3:
        # 1. V sends random γ,η ∈ F.
        setup = self.setup
        beta = self.beta
        gamma = self.gamma
        eta = self.eta
        group_order_N = self.group_order_N
        group_order_n = self.group_order_n

        # 2. compute b_0_gamma
        b_0_gamma = self.B_0_poly.coeff_eval(gamma)
        # compute f_gamma
        f_gamma = self.f_poly.coeff_eval(gamma)
        # 3. compute a_0
        a_0 = self.A_poly.coeff_eval(Scalar(0))
        # 4. compute b_0
        b_0 = group_order_N * a_0 / group_order_n
        # 5. compute b_gamma, and Q_b_gamma
        Z_H_gamma = gamma ** group_order_n - 1
        b_gamma = b_0_gamma * gamma + b_0
        Q_b_gamma = (b_gamma * (f_gamma + beta) - Scalar(1)) / Z_H_gamma

        # 6. batch KZG check
        # (a) both P and V compute v
        v = self.rlc(b_0_gamma, f_gamma, Q_b_gamma)
        # (b) compute commitment: pi_gamma = [h(X)]_1
        h_poly = (self.rlc(self.B_0_poly, self.f_poly, self.Q_B_poly) - v) / (self.x_poly - gamma)
        h_comm_1 = setup.commit(h_poly)

        # 3.7
        # (a) compute a_0_comm_1
        a_0_poly = (self.A_poly - a_0) / self.x_poly
        a_0_comm_1 = setup.commit(a_0_poly)

        return Message3(b_0_gamma, f_gamma, a_0, h_comm_1, a_0_comm_1)

    # random linear combination
    def rlc(self, term_1, term_2, term_3):
        return term_1 + term_2 * self.eta + term_2 * self.eta * self.eta
