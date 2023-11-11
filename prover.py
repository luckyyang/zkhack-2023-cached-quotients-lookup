from setup import *
from typing import Optional
from dataclasses import dataclass
from transcript import Transcript, Message1, Message2, Message3
from poly import Polynomial, Basis
from sortedcontainers import SortedSet
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

    def __init__(self, setup: Setup, table: list):
        self.setup = setup
        self.table = table
        self.powers_of_x = setup.powers_of_x
    def prove(self, witness) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")

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

    def round_1(self, witness) -> Message1:
        return Message1(self.powers_of_x[2])

    def round_2(self) -> Message2:
        setup = self.setup
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
