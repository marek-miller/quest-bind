mod ffi;
use std::ffi::CString;

pub use ffi::{
    bitEncoding as BitEncoding,
    pauliOpType as PauliOpType,
    phaseFunc as PhaseFunc,
    phaseGateType as PhaseGateType,
};

pub type Qreal = f64;

#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub real: Qreal,
    pub imag: Qreal,
}

#[derive(Debug)]
pub struct ComplexMatrix2(ffi::ComplexMatrix2);

#[derive(Debug)]
pub struct ComplexMatrix4(ffi::ComplexMatrix4);

#[derive(Debug)]
pub struct ComplexMatrixN(ffi::ComplexMatrixN);

#[derive(Debug, Copy, Clone)]
pub struct Vector {
    pub x: Qreal,
    pub y: Qreal,
    pub z: Qreal,
}

#[derive(Debug)]
pub struct PauliHamil(ffi::PauliHamil);

#[derive(Debug)]
pub struct DiagonalOp(ffi::DiagonalOp);

#[derive(Debug)]
pub struct Qureg(ffi::Qureg);

#[derive(Debug)]
pub struct QuESTEnv(ffi::QuESTEnv);

pub fn create_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Qureg {
    Qureg(unsafe { ffi::createQureg(num_qubits, env.0) })
}

pub fn create_density_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Qureg {
    Qureg(unsafe { ffi::createDensityQureg(num_qubits, env.0) })
}

pub fn create_clone_qureg(
    qureg: Qureg,
    env: &QuESTEnv,
) -> Qureg {
    Qureg(unsafe { ffi::createCloneQureg(qureg.0, env.0) })
}

pub fn destroy_qureg(
    qureg: Qureg,
    env: &QuESTEnv,
) {
    unsafe { ffi::destroyQureg(qureg.0, env.0) }
}

pub fn create_complex_matrix_n(num_qubits: i32) -> ComplexMatrixN {
    ComplexMatrixN(unsafe { ffi::createComplexMatrixN(num_qubits) })
}

pub fn destroy_complex_matrix_n(matr: ComplexMatrixN) {
    unsafe { ffi::destroyComplexMatrixN(matr.0) }
}

// pub fn init_complex_matrix_2();

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
        // SAFETY:
        // QuEST only copies data pointed to,
        // so the references stay immutable.
        for i in 0..n {
            assert!(real[i].len() >= n);
            real_ptrs.push(real[i].as_ptr() as *mut Qreal);
            assert!(imag[i].len() >= n);
            imag_ptrs.push(imag[i].as_ptr() as *mut Qreal);
        }

        ffi::initComplexMatrixN(
            m.0,
            real_ptrs.as_mut_ptr(),
            imag_ptrs.as_mut_ptr(),
        )
    }
}

pub fn create_pauli_hamil(
    num_qubits: i32,
    num_sum_terms: i32,
) -> PauliHamil {
    PauliHamil(unsafe { ffi::createPauliHamil(num_qubits, num_sum_terms) })
}

pub fn destroy_pauli_hamil(hamil: PauliHamil) {
    unsafe { ffi::destroyPauliHamil(hamil.0) }
}

pub fn create_pauli_hamil_from_file(fn_: &str) -> PauliHamil {
    let filename = CString::new(fn_).unwrap();
    PauliHamil(unsafe { ffi::createPauliHamilFromFile((*filename).as_ptr()) })
}

pub fn init_pauli_hamil(
    hamil: PauliHamil,
    coeffs: &[Qreal],
    codes: &[PauliOpType],
) {
    let hamil_len = hamil.0.numSumTerms as usize;
    assert_eq!(coeffs.len(), hamil_len);
    assert_eq!(codes.len(), hamil_len * hamil.0.numQubits as usize);

    unsafe {
        // SAFETY:
        // QuEST copies values of arrays supplied without modyfings them,
        // so refs. can stay immutable.
        let coeffs_ptr = coeffs.as_ptr() as *mut _;
        let codes_ptr = codes.as_ptr() as *mut _;
        ffi::initPauliHamil(hamil.0, coeffs_ptr, codes_ptr);
    }
}

pub fn create_diagonal_op(
    num_qubits: i32,
    env: &QuESTEnv,
) -> DiagonalOp {
    DiagonalOp(unsafe { ffi::createDiagonalOp(num_qubits, env.0) })
}

pub fn destroy_diagonal_op(
    op: DiagonalOp,
    env: &QuESTEnv,
) {
    unsafe { ffi::destroyDiagonalOp(op.0, env.0) }
}

pub fn sync_diagonal_op(op: &mut DiagonalOp) {
    unsafe { ffi::syncDiagonalOp(op.0) }
}

pub fn init_diagonal_op(
    op: &mut DiagonalOp,
    real: &[Qreal],
    imag: &[Qreal],
) {
    assert!(real.len() >= 2usize.pow(op.0.numQubits as u32));
    assert!(imag.len() >= 2usize.pow(op.0.numQubits as u32));

    unsafe {
        // SAFETY:
        // QuEST copies values of arrays using memcpy,
        // so refs. can stay immutable.
        let real_ptr = real.as_ptr() as *mut _;
        let imag_ptr = imag.as_ptr() as *mut _;
        ffi::initDiagonalOp(op.0, real_ptr, imag_ptr)
    }
}

pub fn init_diagonal_op_from_pauli_hamil(
    op: &mut DiagonalOp,
    hamil: &PauliHamil,
) {
    assert_eq!(op.0.numQubits, hamil.0.numQubits);
    unsafe { ffi::initDiagonalOpFromPauliHamil(op.0, hamil.0) }
}

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
            let row: &[&[f64; 2]; 2] = std::mem::transmute(*m.0.real);
            assert_eq!(row, &[&[1., 2.,], &[3., 4.]]);
        }
        unsafe {
            let row: &[&[f64; 2]; 2] = std::mem::transmute(*m.0.imag);
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
            let row: &[&[f64; 3]; 3] = std::mem::transmute(*m.0.real);
            assert_eq!(row, &[&[1., 2., 3.], &[4., 5., 6.], &[7., 8., 9.]]);
        }
        unsafe {
            let row: &[&[f64; 3]; 3] = std::mem::transmute(*m.0.imag);
            assert_eq!(
                row,
                &[&[11., 12., 13.], &[14., 15., 16.], &[17., 18., 19.]]
            );
        }
        destroy_complex_matrix_n(m);
    }
}
