use qst::{
    init_zero_state,
    QuESTEnv,
    Qureg,
};

pub fn main() {
    let env = QuESTEnv::new();
    env.report();

    let mut qureg = Qureg::new(23, &env);

    {
        let qureg = &mut qureg;

        init_zero_state(qureg);
    }

    // qureg.report();
    qureg.destroy();
}
