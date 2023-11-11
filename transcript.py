from curve import Scalar, G1Point
from merlin.merlin_transcript import MerlinTranscript
from py_ecc.secp256k1.secp256k1 import bytes_to_int
from dataclasses import dataclass


@dataclass
class Message1:
    # [m(x)]₁
    m_comm_1: G1Point


@dataclass
class Message2:
    # [z(x)]₁ (commitment to permutation polynomial)
    z_1: G1Point


@dataclass
class Message3:
    # [quot(x)]₁ (commitment to the quotient polynomial t(X))
    W_t: G1Point


# https://merlin.cool/
class Transcript(MerlinTranscript):
    def append(self, label: bytes, item: bytes) -> None:
        self.append_message(label, item)

    def append_scalar(self, label: bytes, item: Scalar):
        self.append_message(label, item.n.to_bytes(32, "big"))

    def append_point(self, label: bytes, item: G1Point):
        self.append_message(label, item[0].n.to_bytes(32, "big"))
        self.append_message(label, item[1].n.to_bytes(32, "big"))

    def get_and_append_challenge(self, label: bytes) -> Scalar:
        while True:
            challenge_bytes = self.challenge_bytes(label, 255)
            f = Scalar(bytes_to_int(challenge_bytes))
            if f != Scalar.zero():  # Enforce challenge != 0
                self.append(label, challenge_bytes)
                return f

    def round_1(self, message: Message1) -> tuple[Scalar, Scalar]:
        self.append_point(b"m_comm_1", message.m_comm_1)

        # The first two Fiat-Shamir challenges
        beta = self.get_and_append_challenge(b"beta")
        gamma = self.get_and_append_challenge(b"gamma")

        return beta, gamma

    def round_2(self, message: Message2) -> tuple[Scalar, Scalar]:
        self.append_point(b"z_1", message.z_1)

        alpha = self.get_and_append_challenge(b"alpha")
        # This value could be anything, it just needs to be unpredictable. Lets us
        # have evaluation forms at cosets to avoid zero evaluations, so we can
        # divide polys without the 0/0 issue
        fft_cofactor = self.get_and_append_challenge(b"fft_cofactor")

        return alpha, fft_cofactor

    def round_3(self, message: Message3) -> Scalar:
        self.append_point(b"W_t", message.W_t)

        zeta = self.get_and_append_challenge(b"zeta")
        return zeta
