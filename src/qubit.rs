use super::Qureg;
use crate::{
    exceptions::catch_quest_exception,
    ffi,
    QuestError,
};

/// Represents a mutable qubit in a register `Qureg`
#[derive(Debug)]
pub struct Qubit<'a, 'env: 'a> {
    qureg: &'a mut Qureg<'env>,
    index: i32,
}

impl<'a, 'env> Qubit<'a, 'env> {
    pub fn new(
        qureg: &'a mut Qureg<'env>,
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

    /// Measures a single qubit, collapsing it randomly to 0 or 1.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use quest_bind::*;
    /// let env = &QuestEnv::new();
    /// let qureg = &mut Qureg::try_new(2, env).unwrap();
    ///
    /// // Prepare an entangled state `|00> + |11>`
    /// init_zero_state(qureg);
    /// hadamard(qureg, 0).and(controlled_not(qureg, 0, 1)).unwrap();
    ///
    /// // Qubits are entangled now
    /// let outcome0 = qureg.qubit(0).unwrap().measure().unwrap();
    /// let outcome1 = qureg.qubit(1).unwrap().measure().unwrap();
    ///
    /// assert_eq!(outcome0, outcome1);
    /// ```
    ///
    /// See [QuEST API][1] for more information.
    ///
    /// [1]: https://quest-kit.github.io/QuEST/modules.html
    pub fn measure(&mut self) -> Result<i32, QuestError> {
        catch_quest_exception(|| unsafe {
            ffi::measure(self.qureg.reg, self.index)
        })
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
        let outcome = qubit.measure().unwrap();
        assert_eq!(outcome, 0);

        let qubit = &mut Qubit::new(qureg, 1).unwrap();
        let outcome = qubit.measure().unwrap();
        assert_eq!(outcome, 1);
    }
}
