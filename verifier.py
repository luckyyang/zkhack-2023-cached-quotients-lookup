import py_ecc.bn128 as b
from dataclasses import dataclass
from curve import *
from transcript import Transcript
from poly import Polynomial, Basis
from curve import ec_lincomb, G1Point, G2Point


@dataclass
class VerificationKey:
    group_order_N: int
    group_order_n: int
    N_w: Scalar
    n_w: Scalar
    powers_of_x2: list[G2Point]

    def verify_proof(self, pf) -> bool:
        print("Start to verify proof")
        beta, gamma, eta = self.compute_challenges(pf)
        proof = pf.flatten()
        print("proof: ", proof)
        print("beta, gamma, eta: ", beta, gamma, eta)
        # round 2
        # 2.11 verification
        # 2.12 verification
        # round 3
        # 3.6 (c)
        # 3.7 (b)
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

    def rlc(self, term_1, term_2, term_3):
        return term_1 + term_2 * self.eta + term_2 * self.eta * self.eta
