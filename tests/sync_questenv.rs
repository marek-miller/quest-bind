use std::thread;

use quest_bind::{
    get_quest_seeds,
    report_quest_env,
    seed_quest,
    seed_quest_default,
    QuestEnv,
};

// This won't compile if QuestEnv is not Sync
#[test]
fn questenv_is_sync() {
    let env = QuestEnv::new();

    thread::scope(|s| {
        s.spawn(|| report_quest_env(&env));
    });
}

// Assure no data race happens here
#[test]
fn questenv_sync() {
    let env = QuestEnv::new();

    thread::scope(|s| {
        s.spawn(|| env.sync());
        s.spawn(|| env.sync());
        env.sync();
    });
}

#[test]
fn questenv_seed_01() {
    let mut env = QuestEnv::new();
    seed_quest(&mut env, &[1, 2, 3]);

    thread::scope(|s| {
        s.spawn(|| {
            let seeds = get_quest_seeds(&env);
            assert_eq!(seeds, &[1, 2, 3])
        });
    });

    seed_quest(&mut env, &[3, 2, 1]);
    thread::scope(|s| {
        s.spawn(|| {
            let seeds = get_quest_seeds(&env);
            assert_eq!(seeds, &[3, 2, 1])
        });
    });
}

#[test]
fn questenv_seed_quest_default_01() {
    let mut env = QuestEnv::new();
    seed_quest_default(&mut env);

    let seeds = get_quest_seeds(&env);
    let mut seeds_check = get_quest_seeds(&env);
    assert_eq!(seeds, seeds_check);

    thread::scope(|s| {
        s.spawn(|| seeds_check = get_quest_seeds(&env));
    });
    assert_eq!(seeds, seeds_check);
}
