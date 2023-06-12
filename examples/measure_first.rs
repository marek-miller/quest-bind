use quest_bind::{
    init_plus_state,
    measure_with_stats,
    report_quest_env,
    report_qureg_params,
    QuESTEnv,
    QuestError,
    Qureg,
};

fn main() -> Result<(), QuestError> {
    let env = QuESTEnv::new();
    report_quest_env(&env);

    let mut qureg = Qureg::try_new(0x10, &env)?;
    {
        let qureg = &mut qureg;
        init_plus_state(qureg);
        report_qureg_params(qureg);

        let mut outcome_prob = 0.;
        let outcome = measure_with_stats(qureg, 1, &mut outcome_prob)?;

        println!("Measure first qubit.");
        println!("Outcome: {outcome} with prob: {outcome_prob:.2}");
    }
    Ok(())
}
