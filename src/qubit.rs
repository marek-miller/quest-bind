use super::Qureg;

/// Represents a mutable qubit in a register `Qureg`
#[derive(Debug)]
pub struct Qubit<'a, 'b: 'a> {
    qureg: &'a mut Qureg<'b>,
    index: i32,
}

impl<'a, 'b> Qubit<'a, 'b> {
    pub fn new(
        qureg: &'a mut Qureg<'b>,
        index: i32,
    ) -> Option<Self> {
        let num_qubits_represented = qureg.num_qubits_represented();
        if index < 0 || index >= num_qubits_represented {
            None
        } else {
            Some(Self {
                qureg,
                index,
            })
        }
    }

    pub fn index(&self) -> i32 {
        self.index
    }

    pub fn measure(&mut self) -> i32 {
        super::measure(self.qureg, self.index).expect("qubit is a valid index")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        init_zero_state,
        pauli_x,
        QuestEnv,
    };

    #[test]
    fn qubit_init_01() {
        let env = &QuestEnv::new();
        let qureg = &mut Qureg::try_new(2, env).unwrap();
        init_zero_state(qureg);

        assert!(Qubit::new(qureg, 0).is_some());
        assert!(Qubit::new(qureg, 1).is_some());
    }

    #[test]
    fn qubit_init_02() {
        let env = &QuestEnv::new();
        let qureg = &mut Qureg::try_new(2, env).unwrap();
        init_zero_state(qureg);

        assert!(Qubit::new(qureg, 2).is_none());
        assert!(Qubit::new(qureg, 3).is_none());
        assert!(Qubit::new(qureg, -1).is_none());
    }

    #[test]
    fn qubit_into_01() {
        let env = &QuestEnv::new();
        let qureg = &mut Qureg::try_new(2, env).unwrap();
        init_zero_state(qureg);

        assert_eq!(Qubit::new(qureg, 0).unwrap().index(), 0);
        assert_eq!(Qubit::new(qureg, 1).unwrap().index(), 1);
    }

    #[test]
    fn qubit_lifetimes_01() {
        let env = &QuestEnv::new();
        let mut qureg = Qureg::try_new(2, env).unwrap();
        init_zero_state(&mut qureg);

        let qubit = Qubit::new(&mut qureg, 0).unwrap();

        // init_zero_state(&mut qureg);
        // drop(qureg);

        let _ = qubit.index();
    }

    #[test]
    fn qubit_measure_01() {
        let env = &QuestEnv::new();
        let qureg = &mut Qureg::try_new(2, env).unwrap();
        init_zero_state(qureg);
        pauli_x(qureg, 1).unwrap();

        let qubit = &mut Qubit::new(qureg, 0).unwrap();
        let outcome = qubit.measure();
        assert_eq!(outcome, 0);

        let qubit = &mut Qubit::new(qureg, 1).unwrap();
        let outcome = qubit.measure();
        assert_eq!(outcome, 1);
    }
}
