use quest_bind::{
    create_quest_env,
    create_qureg,
    destroy_quest_env,
    destroy_qureg,
    init_plus_state,
    measure_with_stats,
    report_quest_env,
    report_qureg_params,
    Error as QuestError,
};

fn main() -> Result<(), QuestError> {
    let env = create_quest_env();
    report_quest_env(&env);

    let mut qureg = create_qureg(0x10, &env)?;
    {
        let qureg = &mut qureg;
        init_plus_state(qureg);
        report_qureg_params(qureg);

        let mut outcome_prob = 0.;
        let outcome = measure_with_stats(qureg, 1, &mut outcome_prob);

        println!("Measure first qubit.");
        println!("Outcome: {outcome} with prob: {outcome_prob:.2}");
    }

    destroy_qureg(qureg, &env);
    destroy_quest_env(env);
    Ok(())
}
