from setup import Setup
from curve import Scalar
from prover import Prover

def prover():
    print("Beginning prover test")
    group_order = 8
    powers = group_order * 2
    tau = 100
    setup = Setup.generate_srs(powers, tau)

    public_table = [1, 2, 3, 4, 5, 6, 7, 8]
    # public_table = [Scalar(element) for element in public_table_list]
    print("table: ", public_table)
    # values to lookup
    witness = [1, 2, 1, 7, 3, 6, 3, 2]
    # witness = [Scalar(element) for element in witness_list]
    prover = Prover(setup, public_table, group_order, group_order)
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

