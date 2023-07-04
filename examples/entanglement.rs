use quest_bind::{
    controlled_not,
    hadamard,
    init_zero_state,
    report_quest_env,
    report_qureg_params,
    Qubit,
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

    let qb0 = &mut Qubit::new(qureg, 0).unwrap();
    let qb1 = &mut Qubit::new(qureg, 1).unwrap();

    println!("---\nPrepare Bell state: |00> + |11>");
    hadamard(qb0).and(controlled_not(qb0, qb1))?;

    // Measure each qubit
    let outcome0 = qureg.qubit(0).unwrap().measure().unwrap();
    let outcome1 = qureg.qubit(1).unwrap().measure().unwrap();

    println!("Qubit \"0\" measured in state: |{outcome0}>");
    println!("Qubit \"1\" measured in state: |{outcome1}>");
    assert_eq!(outcome0, outcome1);
    println!("They match!");

    // Both `env` and `qureg` are safely discarded here
    Ok(())
}
