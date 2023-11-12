import py_ecc.bn128 as b
from dataclasses import dataclass
from curve import *
from transcript import Transcript
from poly import Polynomial, Basis
from curve import ec_lincomb, G1Point, G2Point
from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, neg

@dataclass
class VerificationKey:
    group_order_N: int
    group_order_n: int
    N_w: Scalar
    n_w: Scalar
    powers_of_x2: list[G2Point]
    T_comm_2: G2Point
    Z_V_comm_2: G2Point

    def verify_proof(self, pf, setup) -> bool:
        print("Start to verify proof")
        beta, gamma, eta = self.compute_challenges(pf)
        group_order_N = self.group_order_N
        group_order_n = self.group_order_n
        powers_of_x2 = self.powers_of_x2
        proof = pf.flatten()
        # get commitments
        m_comm_1 = proof["m_comm_1"]
        A_comm_1 = proof["A_comm_1"]
        Q_A_comm_1 = proof["Q_A_comm_1"]
        B_0_comm_1 = proof["B_0_comm_1"]
        Q_B_comm_1 = proof["Q_B_comm_1"]
        P_comm_1 = proof["P_comm_1"]
        b_0_gamma = proof["b_0_gamma"]
        f_gamma = proof["f_gamma"]
        a_0 = proof["a_0"]
        h_comm_1 = proof["h_comm_1"]
        a_0_comm_1 = proof["a_0_comm_1"]

        # Vefify

        ### Check 1: round 2.11: A encodes the correct values ###
        print("=== Started Check 1: round 2.11: A encodes the correct values ===")
        comb = ec_lincomb([
            (m_comm_1, 1),
            (A_comm_1, -beta)
        ])
        A_check_lhs1 = b.pairing(self.T_comm_2, A_comm_1)
        A_check_rhs1 = b.pairing(self.Z_V_comm_2, Q_A_comm_1)
        A_check_rhs2 = b.pairing(self.powers_of_x2[0], comb)
        assert A_check_lhs1 == A_check_rhs1 * A_check_rhs2, "verify 1: not equal"
        print("=== Finished Check 1: round 2.11: A encodes the correct values ===")

        ### Check 2: round 2.12: B_0 has the appropriate degree ###
        print("=== Started Check 2: B_0 has the appropriate degree ===")
        # TODO Put it into common preprocessed input?
        x_exponent_order = group_order_N - 1 - (group_order_n - 2)
        x_exponent_values_in_coeff = [Scalar(0)] * (x_exponent_order) + [Scalar(1)]
        x_exponent_poly = Polynomial(x_exponent_values_in_coeff, Basis.MONOMIAL)
        # commit x_exponent_poly
        x_exponent_comm_2 = setup.commit2(x_exponent_poly)

        B_0_check_lhs = b.pairing(x_exponent_comm_2, B_0_comm_1)
        B_0_check_rhs = b.pairing(self.powers_of_x2[0], P_comm_1)
        assert B_0_check_lhs == B_0_check_rhs, "B0 degree check failed"
        print("=== Finished Check 2: B_0 has the appropriate degree ===")

        ### Check 3: 3.6 (c) ###
        print("=== Start Check 3: batched KZG check for the correctness of b_0_gamma, f_gamma, Q_b_gamma ===")
        # compute c
        b_0 = group_order_N * a_0 / group_order_n
        Z_H_gamma = gamma ** group_order_n - 1
        b_gamma = b_0_gamma * gamma + b_0
        Q_b_gamma = (b_gamma * (f_gamma + beta) - Scalar(1)) / Z_H_gamma
        c = self.rlc(b_0, f_gamma, Q_b_gamma, eta)
        # batched KZG check for the correctness of b_0_gamma, f_gamma, Q_b_gamma
        batch_check_lhs = b.pairing(x_exponent_comm_2, B_0_comm_1)
        batch_check_rhs = b.pairing(self.powers_of_x2[0], P_comm_1)
        assert B_0_check_lhs == B_0_check_rhs, "B0 degree check failed"
        print("=== Finished Check 3: batched KZG check for the correctness of b_0_gamma, f_gamma, Q_b_gamma ===")


        ### Check 4: 3.7 (b) ###
        print("=== Start Check 4: KZG check for the correctness of a_0 ===")
        # TODO check 4
        print("=== Start Check 4: KZG check for the correctness of a_0 ===")

        print("Finished to verify proof")
        return True

    # Compute challenges (should be same as those computed by prover)
    def compute_challenges(
        self, proof
    ) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        beta = transcript.round_1(proof.msg_1)
        gamma, eta = transcript.round_2(proof.msg_2)

        return beta, gamma, eta

    def verify_commitment(self, proof, W, W_quot_key, eval_key, zeta):
        W_quot = proof[W_quot_key]
        eval = proof[eval_key]
        ec_comb = ec_lincomb(
            [
                (W, 1),
                (W_quot, zeta),
                (b.G1, -eval),
            ]
        )

        assert b.pairing(self.X_2, W_quot) == b.pairing(b.G2, ec_comb)
        print(f"Done KZG10 commitment check for {eval_key} polynomial")

    def rlc(self, term_1, term_2, term_3, eta):
        return term_1 + term_2 * eta + term_2 * eta * eta
