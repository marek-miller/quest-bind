use super::Qureg;
use crate::{
    exceptions::catch_quest_exception,
    ffi,
    QuestError,
};

/// Represents a mutable qubit in a register [`Qureg`][1].
///
/// The type `Qubit` represent a single two-level system in a quantum
/// register `Qureg` at a specified index. Because this
/// qubit might be entangled with other qubits in the register, any
/// operation on the qubit may potentially change the state of the whole quantum
/// register. Hence, in order to even create a qubit, the user has to provide a
/// mutable reference to the underlying register.  A qubit type holding only a
/// shared reference to its `Qureg` would be virtually useless!
///
/// A `Qubit` is bound by two lifetimes: `'a` and `'env`, where `'a` must
/// outlive `'env`.  The lifetime `'a` refers to the mutable reference to a
/// `Qureg` the qubit will hold, and the `'env` is the lifetime of the QuEST
/// global environment the Qureg is bound to.  See also [`Qureg::try_new()`][2]
///
/// [1]: crate::Qureg
/// [2]: crate::Qureg::try_new()
#[derive(Debug)]
pub struct Qubit<'a, 'env: 'a> {
    qureg: &'a mut Qureg<'env>,
    index: i32,
}

impl<'a, 'env> Qubit<'a, 'env> {
    /// Creates a new qubit.
    ///
    /// The `index` must be strictly positive and smaller than
    /// [`qureg.num_qubits_represented()`][1].
    ///
    /// See also [`Qureg::qubit()`][2] for another method to create a `Qubit`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use quest_bind::*;
    /// let env = &QuestEnv::new();
    /// let qureg = &mut Qureg::try_new(2, env).unwrap();
    /// init_zero_state(qureg);
    ///
    /// assert!(Qubit::new(qureg, 0).is_some());
    /// assert!(Qubit::new(qureg, 1).is_some());
    ///
    /// assert!(Qubit::new(qureg, 2).is_none());
    /// assert!(Qubit::new(qureg, -1).is_none());
    /// ```
    ///
    /// # Returns
    ///
    /// Returns a new `Qubit` wrapped in `Some`, if the specified index is
    /// strictly positive and strictly smaller than
    /// [`qureg.num_qubits_represented()`][1].  Otherwise returns `None`.
    ///
    /// [1]: crate::Qureg::num_qubits_represented()
    /// [2]: crate::Qureg::qubit()
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

    /// Index of this qubit in the underlying [`Qureg`](crate::Qureg).
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use quest_bind::*;
    /// let env = &QuestEnv::new();
    /// let qureg = &mut Qureg::try_new(2, env).unwrap();
    /// init_zero_state(qureg);
    ///
    /// let qubit = Qubit::new(qureg, 1).unwrap();
    /// assert_eq!(qubit.index(), 1);
    /// ```
    pub fn index(&self) -> i32 {
        self.index
    }

    /// Measures the qubit, collapsing it randomly to 0 or 1.
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
