use super::Qureg;
use crate::{
    exceptions::catch_quest_exception,
    ffi,
    QuestError,
};

// TODO: Document this!
unsafe impl<'a> Sync for Qureg<'a> {}

/// Represents a mutable qubit in a register [`Qureg`][1].
///
/// The type `Qubit` represent a single two-level system in a quantum
/// register `Qureg` at a specified index. Because this
/// qubit might be entangled with other qubits in the register, any
/// operation on the qubit may potentially change the state of the whole quantum
/// register. TODO: explain shared ref. trick (atomic QuEST API calls).
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
    qureg: &'a Qureg<'env>,
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
        qureg: &'a Qureg<'env>,
        index: i32,
    ) -> Option<Self> {
        if index < 0 || index >= qureg.num_qubits_represented() {
            None
        } else {
            Some(Self::new_unchecked(qureg, index))
        }
    }

    pub fn new_unchecked(
        qureg: &'a Qureg<'env>,
        index: i32,
    ) -> Self {
        Self {
            qureg,
            index,
        }
    }

    /// Check if `other qubit` belongs to the same `Qureg`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use quest_bind::*;
    /// let env = &QuestEnv::new();
    /// let qureg1 = &Qureg::try_new(2, env).unwrap();
    /// let qureg2 = &Qureg::try_new(2, env).unwrap();
    ///
    /// let qubit = qureg1.qubit(0).unwrap();
    /// let other_qubit = qureg1.qubit(1).unwrap();
    /// assert!(qubit.is_same_qureg(&other_qubit));
    ///
    /// let other_qubit = qureg2.qubit(1).unwrap();
    /// assert!(!qubit.is_same_qureg(&other_qubit));
    /// ```
    pub fn is_same_qureg(
        &self,
        other_qubit: &Qubit,
    ) -> bool {
        std::ptr::eq(self.qureg, other_qubit.qureg)
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

/// Apply the single-qubit Hadamard gate.
///
/// # Examples
///
/// ```rust
/// # use quest_bind::*;
/// let env = &QuestEnv::new();
/// let qureg = &mut Qureg::try_new(2, env).unwrap();
/// init_zero_state(qureg);
///
/// let qubit = &mut qureg.qubit(0).unwrap();
/// hadamard(qubit).unwrap();
///
/// let amp = get_real_amp(qureg, 0).unwrap();
/// assert!((amp - SQRT_2.recip()).abs() < EPSILON);
/// ```
///
/// See [QuEST API][1] for more information.
///
/// [1]: https://quest-kit.github.io/QuEST/modules.html
pub fn hadamard(qubit: &mut Qubit) -> Result<(), QuestError> {
    catch_quest_exception(|| unsafe {
        ffi::hadamard(qubit.qureg.reg, qubit.index);
    })
}

/// Apply the controlled not (single control, single target) gate.
///
/// # Examples
///
/// ```rust
/// # use quest_bind::*;
/// let env = &QuestEnv::new();
/// let mut qureg = Qureg::try_new(2, env).unwrap();
/// init_zero_state(&mut qureg);
/// pauli_x(&mut qureg, 0).unwrap();
///
/// let control_qubit = &mut qureg.qubit(0).unwrap();
/// let target_qubit = &mut qureg.qubit(1).unwrap();
/// controlled_not(control_qubit, target_qubit).unwrap();
///
/// let amp = get_real_amp(&mut qureg, 3).unwrap();
/// assert!((amp - 1.).abs() < EPSILON);
/// ```
///
/// # Returns
///
/// This functions returns `QuestError::DifferentQureg`, if
/// `control qubit` and `target_qubit` belong to a different `Qureg`.
///
/// See [QuEST API][1] for more information.
///
/// [1]: https://quest-kit.github.io/QuEST/modules.html
pub fn controlled_not(
    control_qubit: &mut Qubit,
    target_qubit: &mut Qubit,
) -> Result<(), QuestError> {
    if control_qubit.qureg as *const _ != target_qubit.qureg as *const _ {
        return Err(QuestError::DifferentQureg);
    }
    catch_quest_exception(|| unsafe {
        ffi::controlledNot(
            target_qubit.qureg.reg,
            control_qubit.index,
            target_qubit.index,
        );
    })
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
    fn qubit_index_01() {
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

        let qubit = Qubit::new(&qureg, 0).unwrap();

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

    #[test]
    fn controlled_not_different_qureg() {
        let env = &QuestEnv::new();
        let qureg1 = &Qureg::try_new(2, env).unwrap();
        let qureg2 = &Qureg::try_new(2, env).unwrap();

        let qb1 = &mut qureg1.qubit(0).unwrap();
        let qb2 = &mut qureg2.qubit(1).unwrap();

        assert_eq!(controlled_not(qb1, qb2), Err(QuestError::DifferentQureg));
        assert_eq!(controlled_not(qb2, qb1), Err(QuestError::DifferentQureg));
    }
}
