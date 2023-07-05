use quest_bind::*;

#[test]
fn sync_qureg_01() {
    // Qureg must be Sync for this to compile
    let env = &QuestEnv::new();
    let qureg = &Qureg::try_new(2, env).unwrap();

    std::thread::scope(|s| {
        s.spawn(|| {
            report_qureg_params(qureg);
            assert!(!qureg.is_density_matrix());
        });

        report_qureg_params(qureg);
        assert!(!qureg.is_density_matrix());
    });

    report_qureg_params(qureg);
    assert!(!qureg.is_density_matrix());
}
