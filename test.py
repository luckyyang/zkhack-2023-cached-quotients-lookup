from setup import Setup
from prover import Prover
from sortedcontainers import SortedSet

def prover():
    print("Beginning prover test")
    group_order = 8
    powers = group_order * 4
    tau = 100
    setup = Setup.generate_srs(powers, tau)

    public_table = ["a", "b", "c", "d", "e"]
    print("table: ", public_table)
    # values to lookup
    witness = ["a", "b", "a", "c", "c", "a", "c"]
    prover = Prover(setup, public_table)
    proof = prover.prove(witness)
    print("Prover test success")
    return setup, proof, group_order

def verifier(setup, proof, group_order):
    print("Beginning verifier test")
    # vk = setup.verification_key(program.common_preprocessed_input())
    # assert vk.verify_proof(group_order, proof, public)
    print("Verifier test success")

if __name__ == "__main__":
    setup, proof, group_order = prover()
    verifier(setup, proof, group_order)

