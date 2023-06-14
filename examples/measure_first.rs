use quest_bind::{
    init_plus_state,
    measure_with_stats,
    report_quest_env,
    report_qureg_params,
    QuestEnv,
    QuestError,
    Qureg,
};

fn main() -> Result<(), QuestError> {
    let env = QuestEnv::new();
    report_quest_env(&env);

    let qureg = &mut Qureg::try_new(0x10, &env)?;
    report_qureg_params(qureg);
    init_plus_state(qureg);

    let outcome_prob = &mut -1.;
    let outcome = measure_with_stats(qureg, 1, outcome_prob)?;
    println!("Measure first qubit.");
    println!("Outcome: {outcome} with prob.: {outcome_prob:.2}");

    Ok(())
}
