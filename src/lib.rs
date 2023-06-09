mod ffi;
use std::ffi::CString;

pub use ffi::{
    bitEncoding as BitEncoding,
    pauliOpType as PauliOpType,
    phaseFunc as PhaseFunc,
    phaseGateType as PhaseGateType,
};

pub type Qreal = f64;

#[derive(Debug, PartialEq)]
pub enum Error {
    InvalidQuESTInput { err_msg: String, err_func: String },
}

#[derive(Debug)]
pub struct Complex(ffi::Complex);

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
    pub fn is_density_matrix(&self) -> bool {
        self.0.isDensityMatrix != 0
    }

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
) -> Qureg {
    Qureg(unsafe { ffi::createQureg(num_qubits, env.0) })
}

#[must_use]
pub fn create_density_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Qureg {
    Qureg(unsafe { ffi::createDensityQureg(num_qubits, env.0) })
}

#[must_use]
pub fn create_clone_qureg(
    qureg: &Qureg,
    env: &QuESTEnv,
) -> Qureg {
    Qureg(unsafe { ffi::createCloneQureg(qureg.0, env.0) })
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_qureg(
    qureg: Qureg,
    env: &QuESTEnv,
) {
    unsafe { ffi::destroyQureg(qureg.0, env.0) }
}

#[must_use]
pub fn create_complex_matrix_n(num_qubits: i32) -> ComplexMatrixN {
    ComplexMatrixN(unsafe { ffi::createComplexMatrixN(num_qubits) })
}

#[allow(clippy::needless_pass_by_value)]
pub fn destroy_complex_matrix_n(matr: ComplexMatrixN) {
    unsafe { ffi::destroyComplexMatrixN(matr.0) }
}

#[allow(clippy::cast_sign_loss)]
pub fn init_complex_matrix_n(
    m: &mut ComplexMatrixN,
    real: &[&[Qreal]],
    imag: &[&[Qreal]],
) {
    let n = m.0.numQubits as usize;
    assert!(real.len() >= n);
    assert!(imag.len() >= n);

    let mut real_ptrs = Vec::with_capacity(n);
    let mut imag_ptrs = Vec::with_capacity(n);
    unsafe {
        for i in 0..n {
            assert!(real[i].len() >= n);
            real_ptrs.push(real[i].as_ptr());
            assert!(imag[i].len() >= n);
            imag_ptrs.push(imag[i].as_ptr());
        }

        ffi::initComplexMatrixN(m.0, real_ptrs.as_ptr(), imag_ptrs.as_ptr());
    }
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

#[must_use]
pub fn create_pauli_hamil_from_file(fn_: &str) -> PauliHamil {
    let filename = CString::new(fn_).unwrap();
    PauliHamil(unsafe { ffi::createPauliHamilFromFile((*filename).as_ptr()) })
}

#[allow(clippy::cast_sign_loss)]
pub fn init_pauli_hamil(
    hamil: &mut PauliHamil,
    coeffs: &[Qreal],
    codes: &[PauliOpType],
) {
    let hamil_len = hamil.0.numSumTerms as usize;
    assert_eq!(coeffs.len(), hamil_len);
    assert_eq!(codes.len(), hamil_len * hamil.0.numQubits as usize);

    unsafe {
        let coeffs_ptr = coeffs.as_ptr();
        let codes_ptr = codes.as_ptr();
        ffi::initPauliHamil(hamil.0, coeffs_ptr, codes_ptr);
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
        let real_ptr = real.as_ptr();
        let imag_ptr = imag.as_ptr();
        ffi::initDiagonalOp(op.0, real_ptr, imag_ptr);
    }
}

pub fn init_diagonal_op_from_pauli_hamil(
    op: &mut DiagonalOp,
    hamil: &PauliHamil,
) {
    assert_eq!(op.0.numQubits, hamil.0.numQubits);
    unsafe { ffi::initDiagonalOpFromPauliHamil(op.0, hamil.0) }
}

#[must_use]
pub fn create_diagonal_op_from_pauli_hamil_file(
    fn_: &str,
    env: &QuESTEnv,
) -> DiagonalOp {
    let filename = CString::new(fn_).unwrap();
    DiagonalOp(unsafe {
        ffi::createDiagonalOpFromPauliHamilFile((*filename).as_ptr(), env.0)
    })
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

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn init_complex_matrix_n_02() {
        let mut m = create_complex_matrix_n(2);
        init_complex_matrix_n(
            &mut m,
            &[&[1., 2.], &[3., 4.]],
            &[&[11., 12.], &[13., 14.]],
        );

        unsafe {
            let row = &*(*m.0.real).cast::<[&[f64; 2]; 2]>();
            assert_eq!(row, &[&[1., 2.,], &[3., 4.]]);
        }
        unsafe {
            let row = &*(*m.0.imag).cast::<[&[f64; 2]; 2]>();
            assert_eq!(row, &[&[11., 12.], &[13., 14.],]);
        }
        destroy_complex_matrix_n(m);
    }

    #[test]
    fn init_complex_matrix_n_03() {
        let mut m = create_complex_matrix_n(3);
        init_complex_matrix_n(
            &mut m,
            &[&[1., 2., 3.], &[4., 5., 6.], &[7., 8., 9.]],
            &[&[11., 12., 13.], &[14., 15., 16.], &[17., 18., 19.]],
        );

        unsafe {
            let row = &*(*m.0.real).cast::<[&[f64; 3]; 3]>();
            assert_eq!(row, &[&[1., 2., 3.], &[4., 5., 6.], &[7., 8., 9.]]);
        }
        unsafe {
            let row = &*(*m.0.imag).cast::<[&[f64; 3]; 3]>();
            assert_eq!(
                row,
                &[&[11., 12., 13.], &[14., 15., 16.], &[17., 18., 19.]]
            );
        }
        destroy_complex_matrix_n(m);
    }
}
