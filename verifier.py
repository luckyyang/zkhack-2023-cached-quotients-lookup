import py_ecc.bn128 as b
from dataclasses import dataclass
from curve import *
from transcript import Transcript
from poly import Polynomial, Basis


@dataclass
class VerificationKey:
    def verify_proof(self, group_order: int, pf, public=[]) -> bool:
        print("Start to verify proof")
        # round 2
        # 2.11 verification
        # 2.12 verification

        print("Finished to verify proof")
        return True

    # Compute challenges (should be same as those computed by prover)
    def compute_challenges(
        self, proof
    ) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        beta, gamma = transcript.round_1(proof.msg_1)
        alpha, _fft_cofactor = transcript.round_2(proof.msg_2)
        zeta = transcript.round_3(proof.msg_3)

        return beta, gamma, alpha, zeta

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
