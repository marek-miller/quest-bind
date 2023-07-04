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

    let qb0 = &mut Qubit::new(qureg, 0).unwrap();
    let qb1 = &mut Qubit::new(qureg, 1).unwrap();
    
    println!("---\nPrepare Bell state: |00> + |11>");
    hadamard(qb0).and(controlled_not(qb0, qb1))?;

    let (tx, rx) = channel();

    thread::scope(|s| {
        s.spawn(move || {
            println!("[A] Perform local measurement...");
            if let Ok(outcome) = qb0.measure() {
                println!("[A] Qubit measured in state: |{outcome}>");

                println!("[A] Send result via classical channel to B");
                tx.send(outcome).unwrap();
            } else {
                eprintln!("[A] Error while measuring qubit");
            }
        });

        s.spawn(move || {
            println!("[B] Perform local measurement...");
            if let Ok(outcome) = qb1.measure() {
                println!("[B] Qubit measured in state: |{outcome}>");

                let msg = rx.recv().unwrap();
                println!("[B] Receive classical message from A: \"{}\"", msg);

                assert!(msg == outcome);
                println!("[B] They match!");
            } else {
                eprintln!("[B] Error while measuring qubit");
            }
        });
    });

    // Both `env` and `qureg` are safely discarded here
    Ok(())
}
