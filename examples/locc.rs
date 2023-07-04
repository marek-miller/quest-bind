use std::{
    sync::mpsc::channel,
    thread,
};

use quest_bind::*;

fn main() -> Result<(), QuestError> {
    // Initialize QuEST environment
    let env = &QuestEnv::new();
    report_quest_env(env);

    // Initialize 2-qubit register in |00> state
    let qureg = &mut Qureg::try_new(2, env)?;
    report_qureg_params(qureg);
    init_zero_state(qureg);

    println!("---\nPrepare Bell state: |00> + |11>");
    hadamard(qureg, 0).and(controlled_not(qureg, 0, 1))?;

    let mut qb0 = Qubit::new(qureg, 0).unwrap();
    let mut qb1 = Qubit::new(qureg, 1).unwrap();

    let (tx, rx) = channel();
    thread::scope(|s| {
        s.spawn(|| {
            println!("[A] performs local measurement...");
            if let Ok(outcome) = qb0.measure() {
                println!("[A] Qubit measured in state: |{outcome}>");

                println!("[A] Send result via classical channel to B");
                tx.send(outcome).unwrap();
            } else {
                eprint!("[A] Error while measuring qubit");
            }
        });

        s.spawn(move || {
            println!("[B] performs local measurement...");
            if let Ok(outcome) = qb1.measure() {
                println!("[B] Qubit measured in state: |{outcome}>");

                let msg = rx.recv().unwrap();
                println!("[B] Receive classical message from A: \"{}\"", msg);

                assert!(msg == outcome);
                println!("[B] They match!");
            } else {
                eprint!("[B] Error while measuring qubit");
            }
        });
    });

    // Both `env` and `qureg` are safely discarded here
    Ok(())
}
