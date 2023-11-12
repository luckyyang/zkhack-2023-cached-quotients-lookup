from setup import Setup
from program import CommonPreprocessedInput
from curve import Scalar
from prover import Prover

def prover(setup: Setup, group_order_N, group_order_n):
    print("Beginning prover test")

    # public table
    public_table = [1, 2, 3, 4, 5, 6, 7, 8]
    print("table: ", public_table)
    # values to lookup
    witness = [1, 2, 1, 7, 3, 6, 3, 2]

    prover = Prover(setup, public_table, group_order_N, group_order_n)
    proof = prover.prove(witness)
    print("Prover test success")

    return proof

def verifier(setup, proof, group_order_N, group_order_n):
    print("Beginning verifier test")
    common_preprocessed_input = CommonPreprocessedInput(group_order_N, group_order_n)
    vk = setup.verification_key(common_preprocessed_input)
    assert vk.verify_proof(proof)
    print("Verifier test success")

if __name__ == "__main__":
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100
    # group order for witness
    group_order_n = 8
    # group order for public table
    group_order_N = 8
    # number of powers of tau
    powers = group_order_N * 2
    # public table
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, group_order_N, group_order_n)
    # run verifier
    verifier(setup, proof, group_order_N, group_order_n)

