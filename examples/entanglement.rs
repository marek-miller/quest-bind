use quest_bind::{
    controlled_not,
    hadamard,
    init_zero_state,
    measure,
    report_quest_env,
    report_qureg_params,
    QuestEnv,
    QuestError,
    Qureg,
};

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

    // Measure each qubit
    let outcome0 = measure(qureg, 0)?;
    let outcome1 = measure(qureg, 1)?;

    println!("Qubit \"0\" measured in state: |{outcome0}>");
    println!("Qubit \"1\" measured in state: |{outcome1}>");
    assert_eq!(outcome0, outcome1);
    println!("They match!");

    // Both `env` and `qureg` are safely discarded here
    Ok(())
}
