mod ffi;
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
pub struct Qureg<'a> {
    env: &'a QuESTEnv,
    reg: ffi::Qureg,
}

#[derive(Debug)]
pub struct QuESTEnv(ffi::QuESTEnv);

pub fn create_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Qureg<'_> {
    let reg = unsafe { ffi::createQureg(num_qubits, env.0) };
    Qureg {
        env,
        reg,
    }
}

pub fn create_density_qureg(
    num_qubits: i32,
    env: &QuESTEnv,
) -> Qureg {
    let reg = unsafe { ffi::createDensityQureg(num_qubits, env.0) };
    Qureg {
        env,
        reg,
    }
}

pub fn create_clone_qureg(qureg: Qureg) -> Qureg {
    let reg = unsafe { ffi::createCloneQureg(qureg.reg, qureg.env.0) };
    Qureg {
        env: qureg.env,
        reg,
    }
}

pub fn destroy_qureg(qureg: Qureg) {
    unsafe { ffi::destroyQureg(qureg.reg, qureg.env.0) }
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

pub fn init_zero_state(qureg: &mut Qureg) {
    unsafe {
        ffi::initZeroState(qureg.reg);
    }
}

pub fn init_plus_state(qureg: &mut Qureg) {
    unsafe {
        ffi::initPlusState(qureg.reg);
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
