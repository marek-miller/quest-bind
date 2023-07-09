use proptest::prelude::*;
use quest_bind::*;

const TEST_MAX_QUBITS: i32 = 16;

proptest! {
    #[test]
    fn create_diagonal_op_01(n in 1..=TEST_MAX_QUBITS) {
        let env = &QuestEnv::new();
        let _ = DiagonalOp::try_new(n, env).unwrap();

    }

    #[test]
    fn create_diagonal_op_02(n in any::<u32>().prop_map(|x| -(x as i32))) {
        let env = &QuestEnv::new();
        let _ = DiagonalOp::try_new(n, env).unwrap_err();
    }
}

prop_compose! {
    fn strategy_targets()
        (n in 1..=TEST_MAX_QUBITS)
        (
            n in Just(n),
            real in prop::collection::vec(any::<Qreal>(), n as usize),
            imag in prop::collection::vec(any::<Qreal>(), n as usize)
        )
            -> (i32, Vec<Qreal>, Vec<Qreal>) {
        (n, real, imag)
    }
}

proptest! {
    #[test]
    fn set_diagonal_op_elems_01((n, real, imag) in strategy_targets()) {
        let env = &QuestEnv::new();
        let op = &mut DiagonalOp::try_new(n, env).unwrap();

        set_diagonal_op_elems(op, 0, &real, &imag, n as i64).unwrap();
    }
}
