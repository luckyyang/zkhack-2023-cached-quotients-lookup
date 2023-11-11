from setup import Setup
from program import CommonPreprocessedInput
from curve import Scalar
from prover import Prover

def prover():
    print("Beginning prover test")
    group_order_n = 8
    group_order_N = 8
    powers = group_order_N * 2
    tau = 100
    setup = Setup.generate_srs(powers, tau)

    public_table = [1, 2, 3, 4, 5, 6, 7, 8]
    # public_table = [Scalar(element) for element in public_table_list]
    print("table: ", public_table)
    # values to lookup
    witness = [1, 2, 1, 7, 3, 6, 3, 2]
    # witness = [Scalar(element) for element in witness_list]
    prover = Prover(setup, public_table, group_order_N, group_order_n)
    proof = prover.prove(witness)
    print("Prover test success")
    return setup, proof, group_order_N, group_order_n

def verifier(setup, proof, group_order_N, group_order_n):
    print("Beginning verifier test")
    common_preprocessed_input = CommonPreprocessedInput(group_order_N, group_order_n)
    vk = setup.verification_key(common_preprocessed_input)
    assert vk.verify_proof(proof)
    print("Verifier test success")

if __name__ == "__main__":
    setup, proof, group_order_N, group_order_n = prover()
    verifier(setup, proof, group_order_N, group_order_n)

