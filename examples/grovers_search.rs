// This is an adapted example from: grovers_search.c in
// QuEST repository.  QuEST is distributed under MIT License

use std::f64::consts::PI;

use quest_bind::{
    get_prob_amp,
    hadamard,
    init_plus_state,
    multi_controlled_phase_flip,
    pauli_x,
    QuestEnv,
    QuestError,
    Qureg,
};

fn apply_oracle(
    qureg: &mut Qureg,
    num_qubits: i32,
    sol_elem: i64,
) {
    // apply X to transform |111> into |solElem>
    for q in 0..num_qubits {
        if ((sol_elem >> q) & 1) == 0 {
            pauli_x(qureg, q).unwrap();
        }
    }

    // effect |111> -> -|111>
    let ctrls = (0..num_qubits).collect::<Vec<_>>();
    multi_controlled_phase_flip(qureg, &ctrls, num_qubits).unwrap();

    // apply X to transform |solElem> into |111>
    for q in 0..num_qubits {
        if ((sol_elem >> q) & 1) == 0 {
            pauli_x(qureg, q).unwrap();
        }
    }
}

fn apply_diffuser(
    qureg: &mut Qureg,
    num_qubits: i32,
) {
    // apply H to transform |+> into |0>
    for q in 0..num_qubits {
        hadamard(qureg, q).unwrap();
    }

    // apply X to transform |11..1> into |00..0>
    for q in 0..num_qubits {
        pauli_x(qureg, q).unwrap();
    }

    //      // effect |11..1> -> -|11..1>
    let ctrls = (0..num_qubits).collect::<Vec<_>>();
    multi_controlled_phase_flip(qureg, &ctrls, num_qubits).unwrap();

    // apply X to transform |11..1> into |00..0>
    for q in 0..num_qubits {
        pauli_x(qureg, q).unwrap();
    }

    // apply H to transform |+> into |0>
    for q in 0..num_qubits {
        hadamard(qureg, q).unwrap();
    }
}

fn main() -> Result<(), QuestError> {
    //      // prepare the hardware-agnostic QuEST environment
    let env = &QuestEnv::new();

    // choose the system size
    let num_qubits = 15_i32;
    let num_elems = 2_i64.pow(num_qubits as u32);
    let num_reps = (f64::ceil(PI / 4.0 * (num_elems as f64).sqrt())) as usize;

    println!(
        "num_qubits: {num_qubits}, num_elems: {num_elems}, num_reps: \
         {num_reps}"
    );

    // randomly choose the element for which to search
    let sol_elem = 34325 % num_elems;

    // prepare |+>
    let qureg = &mut Qureg::try_new(num_qubits, env).unwrap();
    init_plus_state(qureg);

    // apply Grover's algorithm
    for _ in 0..num_reps {
        apply_oracle(qureg, num_qubits, sol_elem);
        apply_diffuser(qureg, num_qubits);

        // monitor the probability of the solution state
        println!(
            "prob of solution |{sol_elem}> = {}",
            get_prob_amp(qureg, sol_elem).unwrap()
        );
    }

    Ok(())
}
