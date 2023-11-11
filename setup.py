import py_ecc.bn128 as b
from curve import ec_lincomb, G1Point, G2Point, Scalar
from verifier import VerificationKey
from dataclasses import dataclass
from poly import Polynomial, Basis
from program import CommonPreprocessedInput

@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    powers_of_x2: list[G2Point]

    @classmethod
    # tau: a random number whatever you choose
    def generate_srs(cls, powers: int, tau: int):
        print("Start to generate structured reference string")

        # Initialize powers_of_x with 0 values
        powers_of_x = [0] * powers
        # powers_of_x[0] =  b.G1 * tau**0 = b.G1
        # powers_of_x[1] =  b.G1 * tau**1 = powers_of_x[0] * tau
        # powers_of_x[2] =  b.G1 * tau**2 = powers_of_x[1] * tau
        # ...
        # powers_of_x[i] =  b.G1 * tau**i = powers_of_x[i - 1] * tau
        powers_of_x[0] = b.G1

        for i in range(powers):
            if i > 0:
                powers_of_x[i] = b.multiply(powers_of_x[i - 1], tau)

        assert b.is_on_curve(powers_of_x[1], b.b)
        print("Generated G1 side, X^1 point: {}".format(powers_of_x[1]))

        powers_of_x2 = [0] * (powers + 1)
        powers_of_x2[0] = b.G2

        for i in range(powers + 1):
            if i > 0:
                powers_of_x2[i] = b.multiply(powers_of_x2[i - 1], tau)

        assert b.is_on_curve(powers_of_x2[1], b.b2)
        print("Generated G2 side, X^1 point: {}".format(powers_of_x2[1]))

        # assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(powers_of_x2[1], b.G1)
        print("X^1 points checked consistent")
        print("Finished to generate structured reference string")

        return cls(powers_of_x, powers_of_x2)

    # Encodes the KZG commitment that evaluates to the given values in the group
    def commit(self, values: Polynomial) -> G1Point:
        if (values.basis == Basis.LAGRANGE):
            # inverse FFT from Lagrange basis to monomial basis
            coeffs = values.ifft().values
        elif (values.basis == Basis.MONOMIAL):
            coeffs = values.values
        if len(coeffs) > len(self.powers_of_x):
            raise Exception("Not enough powers in setup")
        return ec_lincomb([(s, x) for s, x in zip(self.powers_of_x, coeffs)])

    def setup():
        # 1. generate_srs: will do in the runtime
        # 2. Compute and output [ZV(x)] * G2
        # 3. Compute and output [T(x)] * G2
        # 4. (a): qi = [Qi(x)] * G1
        # 4. (b): [Li(x)] * G1
        # 4. (c): [Li(x)−Li(0) / x] * G1
        print("setup complete")

    def verification_key(self, pk: CommonPreprocessedInput) -> VerificationKey:
        return VerificationKey(
            pk.group_order_N,
            pk.group_order_n,
            Scalar.root_of_unity(pk.group_order_N),
            Scalar.root_of_unity(pk.group_order_n),
            self.powers_of_x2,
        )
