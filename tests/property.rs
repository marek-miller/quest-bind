use proptest::prelude::*;
use quest_bind::*;

const TEST_MAX_QUBITS: i32 = 16;

proptest! {
    #[test]
    fn create_diagonal_op_01(n in 1..TEST_MAX_QUBITS) {
        let env = &QuestEnv::new();
        let _ = DiagonalOp::try_new(n, env).unwrap();

    }

    #[test]
    fn create_diagonal_op_02(n in any::<u32>().prop_map(|x| x as i32 * -1)) {
        let env = &QuestEnv::new();
        let _ = DiagonalOp::try_new(n, env).unwrap_err();
    }
}
