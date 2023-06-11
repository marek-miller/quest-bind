#![allow(clippy::missing_errors_doc)]

use std::ffi::CString;

use exceptions::catch_quest_exception;
pub use ffi::{
    bitEncoding as BitEncoding,
    pauliOpType as PauliOpType,
    phaseFunc as PhaseFunc,
    phaseGateType as PhaseGateType,
};

mod exceptions;
mod ffi;

// TODO: define number abstractions for numerical types
// (use num_traits)
pub type Qreal = f64;

#[derive(Debug, PartialEq, Clone)]
pub enum Error {
    InvalidQuESTInput { err_msg: String, err_func: String },
    NulError(std::ffi::NulError),
    IntoStringError(std::ffi::IntoStringError),
}

#[derive(Debug, Clone, Copy)]
pub struct Complex(ffi::Complex);

impl Complex {
    #[must_use]
    pub fn new(
        real: Qreal,
        imag: Qreal,
    ) -> Self {
        Self(ffi::Complex {
            real,
            imag,
        })
    }

    #[must_use]
    pub fn real(&self) -> Qreal {
        self.0.real
    }

    #[must_use]
    pub fn imag(&self) -> Qreal {
        self.0.imag
    }
}

#[derive(Debug)]
pub struct ComplexMatrix2(ffi::ComplexMatrix2);

#[derive(Debug)]
pub struct ComplexMatrix4(ffi::ComplexMatrix4);

#[derive(Debug)]
pub struct ComplexMatrixN(ffi::ComplexMatrixN);

#[derive(Debug)]
pub struct Vector(ffi::Vector);

#[derive(Debug)]
pub struct PauliHamil(ffi::PauliHamil);

#[derive(Debug)]
pub struct DiagonalOp(ffi::DiagonalOp);

#[derive(Debug)]
pub struct Qureg(ffi::Qureg);

impl Qureg {
    #[must_use]
    pub fn is_density_matrix(&self) -> bool {
        self.0.isDensityMatrix != 0
    }

    #[must_use]
    pub fn num_qubits_represented(&self) -> i32 {
        self.0.numQubitsRepresented
    }
}

#[derive(Debug)]
pub struct QuESTEnv(ffi::QuESTEnv);

#[must_use]
pub fn create_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Result<Qureg, Error> {
    catch_quest_exception(Qureg(unsafe { ffi::createQureg(num_qubits, env.0) }))
}

#[must_use]
pub fn create_density_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Result<Qureg, Error> {
    catch_quest_exception(Qureg(unsafe {
        ffi::createDensityQureg(num_qubits, env.0)
    }))
}

#[must_use]
pub fn create_clone_qureg(
    qureg: &Qureg,
    env: &QuESTEnv,
) -> Result<Qureg, Error> {
    catch_quest_exception(Qureg(unsafe {
        ffi::createCloneQureg(qureg.0, env.0)
    }))
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_qureg(
    qureg: Qureg,
    env: &QuESTEnv,
) {
    unsafe { ffi::destroyQureg(qureg.0, env.0) }
}

#[must_use]
pub fn create_complex_matrix_n(
    num_qubits: i32
) -> Result<ComplexMatrixN, Error> {
    catch_quest_exception(ComplexMatrixN(unsafe {
        ffi::createComplexMatrixN(num_qubits)
    }))
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_complex_matrix_n(matr: ComplexMatrixN) -> Result<(), Error> {
    catch_quest_exception(unsafe { ffi::destroyComplexMatrixN(matr.0) })
}

#[allow(clippy::cast_sign_loss)]
pub fn init_complex_matrix_n(
    m: &mut ComplexMatrixN,
    real: &[&[Qreal]],
    imag: &[&[Qreal]],
) -> Result<(), Error> {
    let n = m.0.numQubits as usize;

    let mut real_ptrs = Vec::with_capacity(n);
    let mut imag_ptrs = Vec::with_capacity(n);
    catch_quest_exception(unsafe {
        for i in 0..n {
            real_ptrs.push(real[i].as_ptr());
            imag_ptrs.push(imag[i].as_ptr());
        }

        ffi::initComplexMatrixN(m.0, real_ptrs.as_ptr(), imag_ptrs.as_ptr())
    })
}

#[must_use]
pub fn create_pauli_hamil(
    num_qubits: i32,
    num_sum_terms: i32,
) -> PauliHamil {
    PauliHamil(unsafe { ffi::createPauliHamil(num_qubits, num_sum_terms) })
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_pauli_hamil(hamil: PauliHamil) {
    unsafe { ffi::destroyPauliHamil(hamil.0) }
}

pub fn create_pauli_hamil_from_file(fn_: &str) -> Result<PauliHamil, Error> {
    let filename = CString::new(fn_).map_err(Error::NulError)?;
    Ok(PauliHamil(unsafe {
        ffi::createPauliHamilFromFile((*filename).as_ptr())
    }))
}

pub fn init_pauli_hamil(
    hamil: &mut PauliHamil,
    coeffs: &[Qreal],
    codes: &[PauliOpType],
) {
    unsafe {
        ffi::initPauliHamil(hamil.0, coeffs.as_ptr(), codes.as_ptr());
    }
}

#[must_use]
pub fn create_diagonal_op(
    num_qubits: i32,
    env: &QuESTEnv,
) -> DiagonalOp {
    DiagonalOp(unsafe { ffi::createDiagonalOp(num_qubits, env.0) })
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_diagonal_op(
    op: DiagonalOp,
    env: &QuESTEnv,
) {
    unsafe {
        ffi::destroyDiagonalOp(op.0, env.0);
    }
}

pub fn sync_diagonal_op(op: &mut DiagonalOp) {
    unsafe {
        ffi::syncDiagonalOp(op.0);
    }
}

/// # Panics
///
/// This function will panic,
/// if either `real` or `imag` have length smaller than
/// `2.pow(num_qubits)`.
#[allow(clippy::cast_sign_loss)]
pub fn init_diagonal_op(
    op: &mut DiagonalOp,
    real: &[Qreal],
    imag: &[Qreal],
) {
    let len_required = 2usize.pow(op.0.numQubits as u32);
    assert!(real.len() >= len_required);
    assert!(imag.len() >= len_required);

    unsafe {
        ffi::initDiagonalOp(op.0, real.as_ptr(), imag.as_ptr());
    }
}

pub fn init_diagonal_op_from_pauli_hamil(
    op: &mut DiagonalOp,
    hamil: &PauliHamil,
) {
    unsafe { ffi::initDiagonalOpFromPauliHamil(op.0, hamil.0) }
}

pub fn create_diagonal_op_from_pauli_hamil_file(
    fn_: &str,
    env: &QuESTEnv,
) -> Result<DiagonalOp, Error> {
    let filename = CString::new(fn_).map_err(Error::NulError)?;
    Ok(DiagonalOp(unsafe {
        ffi::createDiagonalOpFromPauliHamilFile((*filename).as_ptr(), env.0)
    }))
}

/// # Panics
///
/// This function will panic if either
/// `real.len() >= num_elems as usize`, or
/// `imag.len() >= num_elems as usize`.
#[allow(clippy::cast_sign_loss)]
#[allow(clippy::cast_possible_truncation)]
pub fn set_diagonal_op_elems(
    op: &mut DiagonalOp,
    start_ind: i64,
    real: &[Qreal],
    imag: &[Qreal],
    num_elems: i64,
) {
    assert!(real.len() >= num_elems as usize);
    assert!(imag.len() >= num_elems as usize);

    unsafe {
        ffi::setDiagonalOpElems(
            op.0,
            start_ind,
            real.as_ptr(),
            imag.as_ptr(),
            num_elems,
        );
    }
}

pub fn apply_diagonal_op(
    qureg: &mut Qureg,
    op: &DiagonalOp,
) {
    unsafe {
        ffi::applyDiagonalOp(qureg.0, op.0);
    }
}

#[must_use]
pub fn calc_expec_diagonal_op(
    qureg: &Qureg,
    op: &DiagonalOp,
) -> Complex {
    Complex(unsafe { ffi::calcExpecDiagonalOp(qureg.0, op.0) })
}

pub fn report_state(qureg: &Qureg) {
    unsafe { ffi::reportState(qureg.0) }
}

pub fn report_state_to_screen(
    qureg: &Qureg,
    env: &QuESTEnv,
    report_rank: i32,
) {
    unsafe { ffi::reportStateToScreen(qureg.0, env.0, report_rank) }
}

pub fn report_qureg_params(qureg: &Qureg) {
    unsafe {
        ffi::reportQuregParams(qureg.0);
    }
}

pub fn report_pauli_hamil(hamil: &PauliHamil) {
    unsafe {
        ffi::reportPauliHamil(hamil.0);
    }
}

#[must_use]
pub fn get_num_qubits(qureg: &Qureg) -> i32 {
    unsafe { ffi::getNumQubits(qureg.0) }
}

#[must_use]
pub fn get_num_amps(qureg: &Qureg) -> i64 {
    unsafe { ffi::getNumAmps(qureg.0) }
}

pub fn init_blank_state(qureg: &mut Qureg) {
    unsafe {
        ffi::initBlankState(qureg.0);
    }
}

pub fn init_zero_state(qureg: &mut Qureg) {
    unsafe {
        ffi::initZeroState(qureg.0);
    }
}

pub fn init_plus_state(qureg: &mut Qureg) {
    unsafe {
        ffi::initPlusState(qureg.0);
    }
}

pub fn init_classical_state(
    qureg: &mut Qureg,
    state_ind: i64,
) {
    unsafe {
        ffi::initClassicalState(qureg.0, state_ind);
    }
}

pub fn init_pure_state(
    qureg: &mut Qureg,
    pure_: &Qureg,
) {
    unsafe {
        ffi::initPureState(qureg.0, pure_.0);
    }
}

pub fn init_debug_state(qureg: &mut Qureg) {
    unsafe {
        ffi::initDebugState(qureg.0);
    }
}

pub fn init_state_from_amps(
    qureg: &mut Qureg,
    reals: &[Qreal],
    imags: &[Qreal],
) {
    unsafe {
        ffi::initStateFromAmps(qureg.0, reals.as_ptr(), imags.as_ptr());
    }
}

pub fn set_amps(
    qureg: &mut Qureg,
    start_ind: i64,
    reals: &[Qreal],
    imags: &[Qreal],
    num_amps: i64,
) {
    unsafe {
        ffi::setAmps(
            qureg.0,
            start_ind,
            reals.as_ptr(),
            imags.as_ptr(),
            num_amps,
        );
    }
}

pub fn set_density_amps(
    qureg: &mut Qureg,
    start_row: i64,
    start_col: i64,
    reals: &[Qreal],
    imags: &[Qreal],
    num_amps: i64,
) {
    unsafe {
        ffi::setDensityAmps(
            qureg.0,
            start_row,
            start_col,
            reals.as_ptr(),
            imags.as_ptr(),
            num_amps,
        );
    }
}

pub fn clone_qureg(
    target_qureg: &mut Qureg,
    copy_qureg: &Qureg,
) {
    unsafe {
        ffi::cloneQureg(target_qureg.0, copy_qureg.0);
    }
}

pub fn phase_shift(
    qureg: &mut Qureg,
    target_quibit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::phaseShift(qureg.0, target_quibit, angle);
    }
}

pub fn controlled_phase_shift(
    qureg: &mut Qureg,
    id_qubit1: i32,
    id_qubit2: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::controlledPhaseShift(qureg.0, id_qubit1, id_qubit2, angle);
    }
}

pub fn multi_controlled_phase_shift(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    num_control_qubits: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::multiControlledPhaseShift(
            qureg.0,
            control_qubits.as_ptr(),
            num_control_qubits,
            angle,
        );
    }
}

pub fn controlled_phase_flip(
    qureg: &mut Qureg,
    id_qubit1: i32,
    id_qubit2: i32,
) {
    unsafe {
        ffi::controlledPhaseFlip(qureg.0, id_qubit1, id_qubit2);
    }
}

pub fn multi_controlled_phase_flip(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    num_control_qubits: i32,
) {
    unsafe {
        ffi::multiControlledPhaseFlip(
            qureg.0,
            control_qubits.as_ptr(),
            num_control_qubits,
        );
    }
}

pub fn s_gate(
    qureg: &mut Qureg,
    target_qubit: i32,
) {
    unsafe {
        ffi::sGate(qureg.0, target_qubit);
    }
}

pub fn t_gate(
    qureg: &mut Qureg,
    target_qubit: i32,
) {
    unsafe {
        ffi::tGate(qureg.0, target_qubit);
    }
}

#[must_use]
pub fn create_quest_env() -> QuESTEnv {
    QuESTEnv(unsafe { ffi::createQuESTEnv() })
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_quest_env(env: QuESTEnv) {
    unsafe {
        ffi::destroyQuESTEnv(env.0);
    }
}

pub fn sync_quest_env(env: &QuESTEnv) {
    unsafe {
        ffi::syncQuESTEnv(env.0);
    }
}

#[must_use]
pub fn sync_quuest_success(success_code: i32) -> i32 {
    unsafe { ffi::syncQuESTSuccess(success_code) }
}

pub fn report_quest_env(env: &QuESTEnv) {
    unsafe {
        ffi::reportQuESTEnv(env.0);
    }
}

pub fn get_environment_string(env: &QuESTEnv) -> Result<String, Error> {
    let mut cstr =
        CString::new("CUDA=x OpenMP=x MPI=x threads=xxxxxxx ranks=xxxxxxx")
            .map_err(Error::NulError)?;
    unsafe {
        let cstr_ptr = cstr.into_raw();
        ffi::getEnvironmentString(env.0, cstr_ptr);
        cstr = CString::from_raw(cstr_ptr);
    }
    cstr.into_string().map_err(Error::IntoStringError)
}

pub fn copy_state_to_gpu(qureg: &mut Qureg) {
    unsafe {
        ffi::copyStateToGPU(qureg.0);
    }
}

pub fn copy_state_from_gpu(qureg: &mut Qureg) {
    unsafe {
        ffi::copyStateFromGPU(qureg.0);
    }
}

pub fn copy_substate_to_gpu(
    qureg: &mut Qureg,
    start_ind: i64,
    num_amps: i64,
) {
    unsafe {
        ffi::copySubstateToGPU(qureg.0, start_ind, num_amps);
    }
}

pub fn copy_substate_from_gpu(
    qureg: &mut Qureg,
    start_ind: i64,
    num_amps: i64,
) {
    unsafe {
        ffi::copySubstateToGPU(qureg.0, start_ind, num_amps);
    }
}

#[must_use]
pub fn get_amp(
    qureg: &Qureg,
    index: i64,
) -> Complex {
    Complex(unsafe { ffi::getAmp(qureg.0, index) })
}

#[must_use]
pub fn get_real_amp(
    qureg: &Qureg,
    index: i64,
) -> Qreal {
    unsafe { ffi::getRealAmp(qureg.0, index) }
}

#[must_use]
pub fn get_imag_amp(
    qureg: &Qureg,
    index: i64,
) -> Qreal {
    unsafe { ffi::getImagAmp(qureg.0, index) }
}

#[must_use]
pub fn get_prob_amp(
    qureg: &Qureg,
    index: i64,
) -> Qreal {
    unsafe { ffi::getProbAmp(qureg.0, index) }
}

#[must_use]
pub fn get_density_amp(
    qureg: &Qureg,
    row: i64,
    col: i64,
) -> Complex {
    Complex(unsafe { ffi::getDensityAmp(qureg.0, row, col) })
}

#[must_use]
pub fn calc_total_prob(qureg: &Qureg) -> Qreal {
    unsafe { ffi::calcTotalProb(qureg.0) }
}

pub fn compact_unitary(
    qureg: &mut Qureg,
    target_qubit: i32,
    alpha: Complex,
    beta: Complex,
) {
    unsafe {
        ffi::compactUnitary(qureg.0, target_qubit, alpha.0, beta.0);
    }
}

pub fn unitary(
    qureg: &mut Qureg,
    target_qubit: i32,
    u: &ComplexMatrix2,
) {
    unsafe {
        ffi::unitary(qureg.0, target_qubit, u.0);
    }
}

pub fn rotate_x(
    qureg: &mut Qureg,
    rot_qubit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::rotateX(qureg.0, rot_qubit, angle);
    }
}

pub fn rotate_y(
    qureg: &mut Qureg,
    rot_qubit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::rotateY(qureg.0, rot_qubit, angle);
    }
}

pub fn rotate_z(
    qureg: &mut Qureg,
    rot_qubit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::rotateZ(qureg.0, rot_qubit, angle);
    }
}

pub fn rotate_around_axis(
    qureg: &mut Qureg,
    rot_qubit: i32,
    angle: Qreal,
    axis: &Vector,
) {
    unsafe {
        ffi::rotateAroundAxis(qureg.0, rot_qubit, angle, axis.0);
    }
}

pub fn controlled_rotate_x(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::controlledRotateX(qureg.0, control_qubit, target_qubit, angle);
    }
}

pub fn controlled_rotate_y(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::controlledRotateY(qureg.0, control_qubit, target_qubit, angle);
    }
}

pub fn controlled_rotate_z(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::controlledRotateZ(qureg.0, control_qubit, target_qubit, angle);
    }
}

pub fn controlled_rotate_around_axis(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
    angle: Qreal,
    axis: &Vector,
) {
    unsafe {
        ffi::controlledRotateAroundAxis(
            qureg.0,
            control_qubit,
            target_qubit,
            angle,
            axis.0,
        );
    }
}

pub fn controlled_compact_unitary(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
    alpha: Complex,
    beta: Complex,
) {
    unsafe {
        ffi::controlledCompactUnitary(
            qureg.0,
            control_qubit,
            target_qubit,
            alpha.0,
            beta.0,
        );
    }
}

pub fn controlled_unitary(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
    u: &ComplexMatrix2,
) {
    unsafe {
        ffi::controlledUnitary(qureg.0, control_qubit, target_qubit, u.0);
    }
}

pub fn multi_controlled_unitary(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    num_control_qubits: i32,
    target_qubit: i32,
    u: &ComplexMatrix2,
) {
    unsafe {
        ffi::multiControlledUnitary(
            qureg.0,
            control_qubits.as_ptr(),
            num_control_qubits,
            target_qubit,
            u.0,
        );
    }
}

pub fn pauli_x(
    qureg: &mut Qureg,
    target_qubit: i32,
) {
    unsafe {
        ffi::pauliX(qureg.0, target_qubit);
    }
}

pub fn pauli_y(
    qureg: &mut Qureg,
    target_qubit: i32,
) {
    unsafe {
        ffi::pauliY(qureg.0, target_qubit);
    }
}

pub fn pauli_z(
    qureg: &mut Qureg,
    target_qubit: i32,
) {
    unsafe {
        ffi::pauliZ(qureg.0, target_qubit);
    }
}

pub fn hadamard(
    qureg: &mut Qureg,
    target_qubit: i32,
) {
    unsafe {
        ffi::hadamard(qureg.0, target_qubit);
    }
}

pub fn controlled_not(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
) {
    unsafe {
        ffi::controlledNot(qureg.0, control_qubit, target_qubit);
    }
}

pub fn multi_controlled_multi_qubit_not(
    qureg: &mut Qureg,
    ctrls: &[i32],
    num_ctrls: i32,
    targs: &[i32],
    num_targs: i32,
) {
    unsafe {
        ffi::multiControlledMultiQubitNot(
            qureg.0,
            ctrls.as_ptr(),
            num_ctrls,
            targs.as_ptr(),
            num_targs,
        );
    }
}

pub fn multi_qubit_not(
    qureg: &mut Qureg,
    targs: &[i32],
    num_targs: i32,
) {
    unsafe {
        let targs_ptr = targs.as_ptr();
        ffi::multiQubitNot(qureg.0, targs_ptr, num_targs);
    }
}

pub fn controlled_pauli_y(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit: i32,
) {
    unsafe {
        ffi::controlledPauliY(qureg.0, control_qubit, target_qubit);
    }
}

#[must_use]
pub fn calc_prob_of_outcome(
    qureg: &Qureg,
    measure_qubit: i32,
    outcome: i32,
) -> Qreal {
    unsafe { ffi::calcProbOfOutcome(qureg.0, measure_qubit, outcome) }
}

/// # Panics
///
/// This function will panic if
/// `outcome_probs.len() < num_qubits as usize`
#[allow(clippy::cast_sign_loss)]
pub fn calc_prob_of_all_outcomes(
    outcome_probs: &mut [Qreal],
    qureg: &Qureg,
    qubits: &[i32],
    num_qubits: i32,
) {
    assert!(outcome_probs.len() >= num_qubits as usize);
    unsafe {
        ffi::calcProbOfAllOutcomes(
            outcome_probs.as_mut_ptr(),
            qureg.0,
            qubits.as_ptr(),
            num_qubits,
        );
    }
}

pub fn collapse_to_outcome(
    qureg: &mut Qureg,
    measure_qubit: i32,
    outcome: i32,
) -> Qreal {
    unsafe { ffi::collapseToOutcome(qureg.0, measure_qubit, outcome) }
}

pub fn measure(
    qureg: &mut Qureg,
    measure_qubit: i32,
) -> i32 {
    unsafe { ffi::measure(qureg.0, measure_qubit) }
}

pub fn measure_with_stats(
    qureg: &mut Qureg,
    measure_qubit: i32,
    outcome_prob: &mut Qreal,
) -> i32 {
    unsafe {
        let outcome_prob_ptr = outcome_prob as *mut _;
        ffi::measureWithStats(qureg.0, measure_qubit, outcome_prob_ptr)
    }
}

#[must_use]
pub fn calc_inner_product(
    bra: &Qureg,
    ket: &Qureg,
) -> Complex {
    Complex(unsafe { ffi::calcInnerProduct(bra.0, ket.0) })
}

#[must_use]
pub fn calc_density_inner_product(
    rho1: &Qureg,
    rho2: &Qureg,
) -> Qreal {
    unsafe { ffi::calcDensityInnerProduct(rho1.0, rho2.0) }
}

pub fn seed_quest_default(env: &mut QuESTEnv) {
    unsafe {
        let env_ptr = std::ptr::addr_of_mut!(env.0);
        ffi::seedQuESTDefault(env_ptr);
    }
}

pub fn seed_quest(
    env: &mut QuESTEnv,
    seed_array: &[u64],
    num_seeds: i32,
) {
    // QuEST's function signature is `c_ulong`. Let's use u64 for now...
    unsafe {
        let env_ptr = std::ptr::addr_of_mut!(env.0);
        let seed_array_ptr = seed_array.as_ptr();
        ffi::seedQuEST(env_ptr, seed_array_ptr, num_seeds);
    }
}

#[allow(clippy::cast_sign_loss)]
#[must_use]
pub fn get_quest_seeds<'a: 'b, 'b>(env: &'a QuESTEnv) -> (&'b mut [u64], i32) {
    unsafe {
        let seeds_ptr = &mut std::ptr::null_mut();
        let mut num_seeds = 0_i32;
        ffi::getQuESTSeeds(env.0, seeds_ptr, &mut num_seeds);

        let seeds =
            std::slice::from_raw_parts_mut(*seeds_ptr, num_seeds as usize);
        (seeds, num_seeds)
    }
}

pub fn start_recording_qasm(qureg: &mut Qureg) {
    unsafe {
        ffi::startRecordingQASM(qureg.0);
    }
}

pub fn stop_recording_qasm(qureg: &mut Qureg) {
    unsafe {
        ffi::stopRecordingQASM(qureg.0);
    }
}

pub fn clear_recorded_qasm(qureg: &mut Qureg) {
    unsafe {
        ffi::clearRecordedQASM(qureg.0);
    }
}

pub fn print_recorded_qasm(qureg: &mut Qureg) {
    unsafe {
        ffi::printRecordedQASM(qureg.0);
    }
}

pub fn write_recorded_qasm_to_file(
    qureg: &mut Qureg,
    filename: &str,
) -> Result<(), Error> {
    unsafe {
        let filename_cstr = CString::new(filename).map_err(Error::NulError)?;
        ffi::writeRecordedQASMToFile(qureg.0, (*filename_cstr).as_ptr());
    }
    Ok(())
}

pub fn mix_dephasing(
    qureg: &mut Qureg,
    target_qubit: i32,
    prob: Qreal,
) {
    unsafe {
        ffi::mixDephasing(qureg.0, target_qubit, prob);
    }
}

pub fn mix_two_qubit_dephasing(
    qureg: &mut Qureg,
    qubit1: i32,
    qubit2: i32,
    prob: Qreal,
) {
    unsafe {
        ffi::mixTwoQubitDephasing(qureg.0, qubit1, qubit2, prob);
    }
}

pub fn mix_depolarising(
    qureg: &mut Qureg,
    target_qubit: i32,
    prob: Qreal,
) {
    unsafe {
        ffi::mixDepolarising(qureg.0, target_qubit, prob);
    }
}

pub fn mix_damping(
    qureg: &mut Qureg,
    target_qubit: i32,
    prob: Qreal,
) {
    unsafe {
        ffi::mixDamping(qureg.0, target_qubit, prob);
    }
}

pub fn mix_two_qubit_depolarising(
    qureg: &mut Qureg,
    qubit1: i32,
    qubit2: i32,
    prob: Qreal,
) {
    unsafe {
        ffi::mixTwoQubitDepolarising(qureg.0, qubit1, qubit2, prob);
    }
}

pub fn mix_pauli(
    qureg: &mut Qureg,
    target_qubit: i32,
    prob_x: Qreal,
    prob_y: Qreal,
    prob_z: Qreal,
) {
    unsafe {
        ffi::mixPauli(qureg.0, target_qubit, prob_x, prob_y, prob_z);
    }
}

pub fn mix_density_matrix(
    combine_qureg: &mut Qureg,
    prob: Qreal,
    other_qureg: &Qureg,
) {
    unsafe {
        ffi::mixDensityMatrix(combine_qureg.0, prob, other_qureg.0);
    }
}

#[must_use]
pub fn calc_purity(qureg: &Qureg) -> Qreal {
    unsafe { ffi::calcPurity(qureg.0) }
}

#[must_use]
pub fn calc_fidelity(
    qureg: &Qureg,
    pure_state: &Qureg,
) -> Qreal {
    unsafe { ffi::calcFidelity(qureg.0, pure_state.0) }
}

pub fn swap_gate(
    qureg: &mut Qureg,
    qubit1: i32,
    qubit2: i32,
) {
    unsafe {
        ffi::swapGate(qureg.0, qubit1, qubit2);
    }
}

pub fn sqrt_swap_gate(
    qureg: &mut Qureg,
    qb1: i32,
    qb2: i32,
) {
    unsafe {
        ffi::swapGate(qureg.0, qb1, qb2);
    }
}

pub fn multi_state_controlled_unitary(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    control_state: &[i32],
    num_control_qubits: i32,
    target_qubit: i32,
    u: &ComplexMatrix2,
) {
    unsafe {
        ffi::multiStateControlledUnitary(
            qureg.0,
            control_qubits.as_ptr(),
            control_state.as_ptr(),
            num_control_qubits,
            target_qubit,
            u.0,
        );
    }
}

pub fn multi_rotate_z(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::multiRotateZ(qureg.0, qubits.as_ptr(), num_qubits, angle);
    }
}

pub fn multi_rotate_pauli(
    qureg: &mut Qureg,
    target_qubits: &[i32],
    target_paulis: &[PauliOpType],
    num_targets: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::multiRotatePauli(
            qureg.0,
            target_qubits.as_ptr(),
            target_paulis.as_ptr(),
            num_targets,
            angle,
        );
    }
}

pub fn multi_controlled_multi_rotate_z(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    num_controls: i32,
    target_qubits: &[i32],
    num_targets: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::multiControlledMultiRotateZ(
            qureg.0,
            control_qubits.as_ptr(),
            num_controls,
            target_qubits.as_ptr(),
            num_targets,
            angle,
        );
    }
}

pub fn multi_controlled_multi_rotate_pauli(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    num_controls: i32,
    target_qubits: &[i32],
    target_paulis: &[PauliOpType],
    num_targets: i32,
    angle: Qreal,
) {
    unsafe {
        ffi::multiControlledMultiRotatePauli(
            qureg.0,
            control_qubits.as_ptr(),
            num_controls,
            target_qubits.as_ptr(),
            target_paulis.as_ptr(),
            num_targets,
            angle,
        );
    }
}

pub fn calc_expec_pauli_prod(
    qureg: &Qureg,
    target_qubits: &[i32],
    pauli_codes: &[PauliOpType],
    num_targets: i32,
    workspace: &mut Qureg,
) -> Qreal {
    unsafe {
        ffi::calcExpecPauliProd(
            qureg.0,
            target_qubits.as_ptr(),
            pauli_codes.as_ptr(),
            num_targets,
            workspace.0,
        )
    }
}

pub fn calc_expec_pauli_sum(
    qureg: &Qureg,
    all_pauli_codes: &[PauliOpType],
    term_coeffs: &[Qreal],
    num_sum_terms: i32,
    workspace: &mut Qureg,
) -> Qreal {
    unsafe {
        ffi::calcExpecPauliSum(
            qureg.0,
            all_pauli_codes.as_ptr(),
            term_coeffs.as_ptr(),
            num_sum_terms,
            workspace.0,
        )
    }
}

pub fn calc_expec_pauli_hamil(
    qureg: &Qureg,
    hamil: &PauliHamil,
    workspace: &mut Qureg,
) -> Qreal {
    unsafe { ffi::calcExpecPauliHamil(qureg.0, hamil.0, workspace.0) }
}

pub fn two_qubit_unitary(
    qureg: &mut Qureg,
    target_qubit1: i32,
    target_qubit2: i32,
    u: &ComplexMatrix4,
) {
    unsafe {
        ffi::twoQubitUnitary(qureg.0, target_qubit1, target_qubit2, u.0);
    }
}

pub fn controlled_two_qubit_unitary(
    qureg: &mut Qureg,
    control_qubit: i32,
    target_qubit1: i32,
    target_qubit2: i32,
    u: &ComplexMatrix4,
) {
    unsafe {
        ffi::controlledTwoQubitUnitary(
            qureg.0,
            control_qubit,
            target_qubit1,
            target_qubit2,
            u.0,
        );
    }
}

pub fn multi_controlled_two_qubit_unitary(
    qureg: &mut Qureg,
    control_qubits: &[i32],
    num_control_qubits: i32,
    target_qubit1: i32,
    target_qubit2: i32,
    u: &ComplexMatrix4,
) {
    unsafe {
        ffi::multiControlledTwoQubitUnitary(
            qureg.0,
            control_qubits.as_ptr(),
            num_control_qubits,
            target_qubit1,
            target_qubit2,
            u.0,
        );
    }
}

pub fn multi_qubit_unitary(
    qureg: &mut Qureg,
    targs: &[i32],
    num_targs: i32,
    u: &ComplexMatrixN,
) {
    unsafe {
        ffi::multiQubitUnitary(qureg.0, targs.as_ptr(), num_targs, u.0);
    }
}

pub fn controlled_multi_qubit_unitary(
    qureg: &mut Qureg,
    ctrl: i32,
    targs: &[i32],
    num_targs: i32,
    u: &ComplexMatrixN,
) {
    unsafe {
        ffi::controlledMultiQubitUnitary(
            qureg.0,
            ctrl,
            targs.as_ptr(),
            num_targs,
            u.0,
        );
    }
}

pub fn multi_controlled_multi_qubit_unitary(
    qureg: &mut Qureg,
    ctrls: &[i32],
    num_ctrls: i32,
    targs: &[i32],
    num_targs: i32,
    u: &ComplexMatrixN,
) {
    unsafe {
        ffi::multiControlledMultiQubitUnitary(
            qureg.0,
            ctrls.as_ptr(),
            num_ctrls,
            targs.as_ptr(),
            num_targs,
            u.0,
        );
    }
}

pub fn mix_kraus_map(
    qureg: &mut Qureg,
    target: i32,
    ops: &[ComplexMatrix2],
    num_ops: i32,
) {
    let ops_inner = ops.iter().map(|x| x.0).collect::<Vec<_>>();
    unsafe {
        ffi::mixKrausMap(qureg.0, target, ops_inner.as_ptr(), num_ops);
    }
}

pub fn mix_two_qubit_kraus_map(
    qureg: &mut Qureg,
    target1: i32,
    target2: i32,
    ops: &[ComplexMatrix4],
    num_ops: i32,
) {
    let ops_inner = ops.iter().map(|x| x.0).collect::<Vec<_>>();
    unsafe {
        ffi::mixTwoQubitKrausMap(
            qureg.0,
            target1,
            target2,
            ops_inner.as_ptr(),
            num_ops,
        );
    }
}

pub fn mix_multi_qubit_kraus_map(
    qureg: &mut Qureg,
    targets: &[i32],
    num_targets: i32,
    ops: &[ComplexMatrixN],
    num_ops: i32,
) {
    let ops_inner = ops.iter().map(|x| x.0).collect::<Vec<_>>();
    unsafe {
        ffi::mixMultiQubitKrausMap(
            qureg.0,
            targets.as_ptr(),
            num_targets,
            ops_inner.as_ptr(),
            num_ops,
        );
    }
}

pub fn mix_nontp_kraus_map(
    qureg: &mut Qureg,
    target: i32,
    ops: &[ComplexMatrix2],
    num_ops: i32,
) {
    let ops_inner = ops.iter().map(|x| x.0).collect::<Vec<_>>();
    unsafe {
        ffi::mixNonTPKrausMap(qureg.0, target, ops_inner.as_ptr(), num_ops);
    }
}

pub fn mix_nontp_two_qubit_kraus_map(
    qureg: &mut Qureg,
    target1: i32,
    target2: i32,
    ops: &[ComplexMatrix4],
    num_ops: i32,
) {
    let ops_inner = ops.iter().map(|x| x.0).collect::<Vec<_>>();
    unsafe {
        ffi::mixNonTPTwoQubitKrausMap(
            qureg.0,
            target1,
            target2,
            ops_inner.as_ptr(),
            num_ops,
        );
    }
}

pub fn mix_nontp_multi_qubit_kraus_map(
    qureg: &mut Qureg,
    targets: &[i32],
    num_targets: i32,
    ops: &[ComplexMatrixN],
    num_ops: i32,
) {
    let ops_inner = ops.iter().map(|x| x.0).collect::<Vec<_>>();
    unsafe {
        ffi::mixNonTPMultiQubitKrausMap(
            qureg.0,
            targets.as_ptr(),
            num_targets,
            ops_inner.as_ptr(),
            num_ops,
        );
    }
}

#[must_use]
pub fn calc_hilbert_schmidt_distance(
    a: &Qureg,
    b: &Qureg,
) -> Qreal {
    unsafe { ffi::calcHilbertSchmidtDistance(a.0, b.0) }
}

pub fn set_weighted_qureg(
    fac1: Complex,
    qureg1: &Qureg,
    fac2: Complex,
    qureg2: &Qureg,
    fac_out: Complex,
    out: &mut Qureg,
) {
    unsafe {
        ffi::setWeightedQureg(
            fac1.0, qureg1.0, fac2.0, qureg2.0, fac_out.0, out.0,
        );
    }
}

pub fn apply_pauli_sum(
    in_qureg: &Qureg,
    all_pauli_codes: &[PauliOpType],
    term_coeffs: &[Qreal],
    num_sum_terms: i32,
    out_qureg: &mut Qureg,
) {
    unsafe {
        ffi::applyPauliSum(
            in_qureg.0,
            all_pauli_codes.as_ptr(),
            term_coeffs.as_ptr(),
            num_sum_terms,
            out_qureg.0,
        );
    }
}

pub fn apply_pauli_hamil(
    in_qureg: &Qureg,
    hamil: &PauliHamil,
    out_qureg: &mut Qureg,
) {
    unsafe {
        ffi::applyPauliHamil(in_qureg.0, hamil.0, out_qureg.0);
    }
}

pub fn apply_trotter_circuitit(
    qureg: &mut Qureg,
    hamil: &PauliHamil,
    time: Qreal,
    order: i32,
    reps: i32,
) {
    unsafe {
        ffi::applyTrotterCircuit(qureg.0, hamil.0, time, order, reps);
    }
}

pub fn apply_matrix2(
    qureg: &mut Qureg,
    target_qubit: i32,
    u: &ComplexMatrix2,
) {
    unsafe {
        ffi::applyMatrix2(qureg.0, target_qubit, u.0);
    }
}

pub fn apply_matrix4(
    qureg: &mut Qureg,
    target_qubit1: i32,
    target_qubit2: i32,
    u: &ComplexMatrix4,
) {
    unsafe {
        ffi::applyMatrix4(qureg.0, target_qubit1, target_qubit2, u.0);
    }
}

pub fn apply_matrix_n(
    qureg: &mut Qureg,
    targs: &[i32],
    num_targs: i32,
    u: &ComplexMatrixN,
) {
    unsafe {
        ffi::applyMatrixN(qureg.0, targs.as_ptr(), num_targs, u.0);
    }
}

pub fn apply_multi_controlled_matrix_n(
    qureg: &mut Qureg,
    ctrls: &[i32],
    num_ctrls: i32,
    targs: &[i32],
    num_targs: i32,
    u: &ComplexMatrixN,
) {
    unsafe {
        ffi::applyMultiControlledMatrixN(
            qureg.0,
            ctrls.as_ptr(),
            num_ctrls,
            targs.as_ptr(),
            num_targs,
            u.0,
        );
    }
}

pub fn apply_phase_func(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits: i32,
    encoding: BitEncoding,
    coeffs: &[Qreal],
    exponents: &[Qreal],
    num_terms: i32,
) {
    unsafe {
        ffi::applyPhaseFunc(
            qureg.0,
            qubits.as_ptr(),
            num_qubits,
            encoding,
            coeffs.as_ptr(),
            exponents.as_ptr(),
            num_terms,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn apply_phase_func_overrides(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits: i32,
    encoding: BitEncoding,
    coeffs: &[Qreal],
    exponents: &[Qreal],
    num_terms: i32,
    override_inds: &[i64],
    override_phases: &[Qreal],
    num_overrides: i32,
) {
    unsafe {
        ffi::applyPhaseFuncOverrides(
            qureg.0,
            qubits.as_ptr(),
            num_qubits,
            encoding,
            coeffs.as_ptr(),
            exponents.as_ptr(),
            num_terms,
            override_inds.as_ptr(),
            override_phases.as_ptr(),
            num_overrides,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn apply_multi_var_phase_func(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits_per_reg: &[i32],
    num_regs: i32,
    encoding: BitEncoding,
    coeffs: &[Qreal],
    exponents: &[Qreal],
    num_terms_per_reg: &[i32],
) {
    unsafe {
        ffi::applyMultiVarPhaseFunc(
            qureg.0,
            qubits.as_ptr(),
            num_qubits_per_reg.as_ptr(),
            num_regs,
            encoding,
            coeffs.as_ptr(),
            exponents.as_ptr(),
            num_terms_per_reg.as_ptr(),
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn apply_multi_var_phase_func_overrides(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits_per_reg: &[i32],
    num_regs: i32,
    encoding: BitEncoding,
    coeffs: &[Qreal],
    exponents: &[Qreal],
    num_terms_per_reg: &[i32],
    override_inds: &[i64],
    override_phases: &[Qreal],
    num_overrides: i32,
) {
    unsafe {
        ffi::applyMultiVarPhaseFuncOverrides(
            qureg.0,
            qubits.as_ptr(),
            num_qubits_per_reg.as_ptr(),
            num_regs,
            encoding,
            coeffs.as_ptr(),
            exponents.as_ptr(),
            num_terms_per_reg.as_ptr(),
            override_inds.as_ptr(),
            override_phases.as_ptr(),
            num_overrides,
        );
    }
}

pub fn apply_named_phase_func(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits_per_reg: &[i32],
    num_regs: i32,
    encoding: BitEncoding,
    function_name_code: PhaseFunc,
) {
    unsafe {
        ffi::applyNamedPhaseFunc(
            qureg.0,
            qubits.as_ptr(),
            num_qubits_per_reg.as_ptr(),
            num_regs,
            encoding,
            function_name_code,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn apply_named_phase_func_overrides(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits_per_reg: &[i32],
    num_regs: i32,
    encoding: BitEncoding,
    function_name_code: PhaseFunc,
    override_inds: &[i64],
    override_phases: &[Qreal],
    num_overrides: i32,
) {
    unsafe {
        ffi::applyNamedPhaseFuncOverrides(
            qureg.0,
            qubits.as_ptr(),
            num_qubits_per_reg.as_ptr(),
            num_regs,
            encoding,
            function_name_code,
            override_inds.as_ptr(),
            override_phases.as_ptr(),
            num_overrides,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn apply_param_named_phase_func(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits_per_reg: &[i32],
    num_regs: i32,
    encoding: BitEncoding,
    function_name_code: PhaseFunc,
    params: &[Qreal],
    num_params: i32,
) {
    unsafe {
        ffi::applyParamNamedPhaseFunc(
            qureg.0,
            qubits.as_ptr(),
            num_qubits_per_reg.as_ptr(),
            num_regs,
            encoding,
            function_name_code,
            params.as_ptr(),
            num_params,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn apply_param_named_phase_func_overrides(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits_per_reg: &[i32],
    num_regs: i32,
    encoding: BitEncoding,
    function_name_code: PhaseFunc,
    params: &[Qreal],
    num_params: i32,
    override_inds: &[i64],
    override_phases: &[Qreal],
    num_overrides: i32,
) {
    unsafe {
        ffi::applyParamNamedPhaseFuncOverrides(
            qureg.0,
            qubits.as_ptr(),
            num_qubits_per_reg.as_ptr(),
            num_regs,
            encoding,
            function_name_code,
            params.as_ptr(),
            num_params,
            override_inds.as_ptr(),
            override_phases.as_ptr(),
            num_overrides,
        );
    }
}

pub fn apply_full_qft(qureg: &mut Qureg) {
    unsafe {
        ffi::applyFullQFT(qureg.0);
    }
}

pub fn apply_qft(
    qureg: &mut Qureg,
    qubits: &[i32],
    num_qubits: i32,
) {
    unsafe {
        ffi::applyQFT(qureg.0, qubits.as_ptr(), num_qubits);
    }
}

pub fn apply_projector(
    qureg: &mut Qureg,
    qubit: i32,
    outcome: i32,
) {
    unsafe {
        ffi::applyProjector(qureg.0, qubit, outcome);
    }
}

#[cfg(test)]
mod tests;
