//! FFI bindings for QuEST v3.5.0

#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]

use std::ffi::{
    c_char,
    c_double,
    c_int,
    c_longlong,
    c_ulong,
};

type qreal = c_double;

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub enum phaseGateType {
    SIGMA_Z = 0,
    S_GATE = 1,
    T_GATE = 2,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct QASMLogger {
    /// generated QASM string
    buffer:     *mut c_char,
    /// maximum number of chars before overflow
    bufferSize: c_int,
    /// number of chars currently in buffer
    bufferFill: c_int,
    /// whether gates are being added to buffer
    isLogging:  c_int,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct ComplexArray {
    real: *mut qreal,
    imag: *mut qreal,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub enum pauliOpType {
    PAULI_I = 0,
    PAULI_X = 1,
    PAULI_Y = 2,
    PAULI_Z = 3,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Complex {
    real: qreal,
    imag: qreal,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct ComplexMatrix2 {
    pub real: [[qreal; 2usize]; 2usize],
    pub imag: [[qreal; 2usize]; 2usize],
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct ComplexMatrix4 {
    pub real: [[qreal; 4usize]; 4usize],
    pub imag: [[qreal; 4usize]; 4usize],
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct ComplexMatrixN {
    pub numQubits: c_int,
    pub real:      *mut *mut qreal,
    pub imag:      *mut *mut qreal,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct Vector {
    pub x: qreal,
    pub y: qreal,
    pub z: qreal,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub enum phaseFunc {
    NORM = 0,
    SCALED_NORM = 1,
    INVERSE_NORM = 2,
    SCALED_INVERSE_NORM = 3,
    SCALED_INVERSE_SHIFTED_NORM = 4,
    PRODUCT = 5,
    SCALED_PRODUCT = 6,
    INVERSE_PRODUCT = 7,
    SCALED_INVERSE_PRODUCT = 8,
    DISTANCE = 9,
    SCALED_DISTANCE = 10,
    INVERSE_DISTANCE = 11,
    SCALED_INVERSE_DISTANCE = 12,
    SCALED_INVERSE_SHIFTED_DISTANCE = 13,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub enum bitEncoding {
    UNSIGNED = 0,
    TWOS_COMPLEMENT = 1,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct PauliHamil {
    pub pauliCodes:  *mut pauliOpType,
    pub termCoeffs:  *mut qreal,
    pub numSumTerms: c_int,
    pub numQubits:   c_int,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct DiagonalOp {
    pub numQubits:        c_int,
    pub numElemsPerChunk: c_longlong,
    pub numChunks:        c_int,
    pub chunkId:          c_int,
    pub real:             *mut qreal,
    pub imag:             *mut qreal,
    pub deviceOperator:   ComplexArray,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Qureg {
    isDensityMatrix:      c_int,
    numQubitsRepresented: c_int,
    numQubitsInStateVec:  c_int,
    numAmpsPerChunk:      c_longlong,
    numAmpsTotal:         c_longlong,
    chunkId:              c_int,

    numChunks: c_int,

    stateVec:     ComplexArray,
    pairStateVec: ComplexArray,

    deviceStateVec:       ComplexArray,
    firstLevelReduction:  *mut qreal,
    secondLevelReduction: *mut qreal,

    qasmLog: *mut QASMLogger,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct QuESTEnv {
    rank:     c_int,
    numRanks: c_int,
    seeds:    *mut c_ulong,
    numSeeds: c_int,
}

#[link(name = "QuEST")]
extern "C" {
    pub fn createQuESTEnv() -> QuESTEnv;
    pub fn reportQuESTEnv(env: QuESTEnv);
    pub fn destroyQuESTEnv(env: QuESTEnv);
    pub fn createQureg(
        numQubits: c_int,
        env: QuESTEnv,
    ) -> Qureg;
    pub fn destroyQureg(
        qureg: Qureg,
        env: QuESTEnv,
    );
    pub fn initZeroState(qureg: Qureg);
    pub fn initPlusState(qureg: Qureg);
    pub fn reportState(qureg: Qureg);
}

#[cfg(test)]
mod tests {
    #![allow(deref_nullptr)]

    use super::*;

    #[test]
    fn bindgen_test_layout_QASMLogger() {
        assert_eq!(
            ::std::mem::size_of::<QASMLogger>(),
            24usize,
            concat!("Size of: ", stringify!(QASMLogger))
        );
        assert_eq!(
            ::std::mem::align_of::<QASMLogger>(),
            8usize,
            concat!("Alignment of ", stringify!(QASMLogger))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QASMLogger>())).buffer as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(QASMLogger),
                "::",
                stringify!(buffer)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QASMLogger>())).bufferSize as *const _
                    as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(QASMLogger),
                "::",
                stringify!(bufferSize)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QASMLogger>())).bufferFill as *const _
                    as usize
            },
            12usize,
            concat!(
                "Offset of field: ",
                stringify!(QASMLogger),
                "::",
                stringify!(bufferFill)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QASMLogger>())).isLogging as *const _
                    as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(QASMLogger),
                "::",
                stringify!(isLogging)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_ComplexArray() {
        assert_eq!(
            ::std::mem::size_of::<ComplexArray>(),
            16usize,
            concat!("Size of: ", stringify!(ComplexArray))
        );
        assert_eq!(
            ::std::mem::align_of::<ComplexArray>(),
            8usize,
            concat!("Alignment of ", stringify!(ComplexArray))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexArray>())).real as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexArray),
                "::",
                stringify!(real)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexArray>())).imag as *const _
                    as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexArray),
                "::",
                stringify!(imag)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_Complex() {
        assert_eq!(
            ::std::mem::size_of::<Complex>(),
            16usize,
            concat!("Size of: ", stringify!(Complex))
        );
        assert_eq!(
            ::std::mem::align_of::<Complex>(),
            8usize,
            concat!("Alignment of ", stringify!(Complex))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Complex>())).real as *const _ as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(Complex),
                "::",
                stringify!(real)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Complex>())).imag as *const _ as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(Complex),
                "::",
                stringify!(imag)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_ComplexMatrix2() {
        assert_eq!(
            ::std::mem::size_of::<ComplexMatrix2>(),
            64usize,
            concat!("Size of: ", stringify!(ComplexMatrix2))
        );
        assert_eq!(
            ::std::mem::align_of::<ComplexMatrix2>(),
            8usize,
            concat!("Alignment of ", stringify!(ComplexMatrix2))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrix2>())).real as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrix2),
                "::",
                stringify!(real)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrix2>())).imag as *const _
                    as usize
            },
            32usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrix2),
                "::",
                stringify!(imag)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_ComplexMatrix4() {
        assert_eq!(
            ::std::mem::size_of::<ComplexMatrix4>(),
            256usize,
            concat!("Size of: ", stringify!(ComplexMatrix4))
        );
        assert_eq!(
            ::std::mem::align_of::<ComplexMatrix4>(),
            8usize,
            concat!("Alignment of ", stringify!(ComplexMatrix4))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrix4>())).real as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrix4),
                "::",
                stringify!(real)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrix4>())).imag as *const _
                    as usize
            },
            128usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrix4),
                "::",
                stringify!(imag)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_ComplexMatrixN() {
        assert_eq!(
            ::std::mem::size_of::<ComplexMatrixN>(),
            24usize,
            concat!("Size of: ", stringify!(ComplexMatrixN))
        );
        assert_eq!(
            ::std::mem::align_of::<ComplexMatrixN>(),
            8usize,
            concat!("Alignment of ", stringify!(ComplexMatrixN))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrixN>())).numQubits as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrixN),
                "::",
                stringify!(numQubits)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrixN>())).real as *const _
                    as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrixN),
                "::",
                stringify!(real)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<ComplexMatrixN>())).imag as *const _
                    as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(ComplexMatrixN),
                "::",
                stringify!(imag)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_Vector() {
        assert_eq!(
            ::std::mem::size_of::<Vector>(),
            24usize,
            concat!("Size of: ", stringify!(Vector))
        );
        assert_eq!(
            ::std::mem::align_of::<Vector>(),
            8usize,
            concat!("Alignment of ", stringify!(Vector))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Vector>())).x as *const _ as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(Vector),
                "::",
                stringify!(x)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Vector>())).y as *const _ as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(Vector),
                "::",
                stringify!(y)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Vector>())).z as *const _ as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(Vector),
                "::",
                stringify!(z)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_PauliHamil() {
        assert_eq!(
            ::std::mem::size_of::<PauliHamil>(),
            24usize,
            concat!("Size of: ", stringify!(PauliHamil))
        );
        assert_eq!(
            ::std::mem::align_of::<PauliHamil>(),
            8usize,
            concat!("Alignment of ", stringify!(PauliHamil))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<PauliHamil>())).pauliCodes as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(PauliHamil),
                "::",
                stringify!(pauliCodes)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<PauliHamil>())).termCoeffs as *const _
                    as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(PauliHamil),
                "::",
                stringify!(termCoeffs)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<PauliHamil>())).numSumTerms as *const _
                    as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(PauliHamil),
                "::",
                stringify!(numSumTerms)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<PauliHamil>())).numQubits as *const _
                    as usize
            },
            20usize,
            concat!(
                "Offset of field: ",
                stringify!(PauliHamil),
                "::",
                stringify!(numQubits)
            )
        );
    }
    #[test]
    fn bindgen_test_layout_DiagonalOp() {
        assert_eq!(
            ::std::mem::size_of::<DiagonalOp>(),
            56usize,
            concat!("Size of: ", stringify!(DiagonalOp))
        );
        assert_eq!(
            ::std::mem::align_of::<DiagonalOp>(),
            8usize,
            concat!("Alignment of ", stringify!(DiagonalOp))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).numQubits as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(numQubits)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).numElemsPerChunk
                    as *const _ as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(numElemsPerChunk)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).numChunks as *const _
                    as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(numChunks)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).chunkId as *const _
                    as usize
            },
            20usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(chunkId)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).real as *const _ as usize
            },
            24usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(real)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).imag as *const _ as usize
            },
            32usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(imag)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<DiagonalOp>())).deviceOperator
                    as *const _ as usize
            },
            40usize,
            concat!(
                "Offset of field: ",
                stringify!(DiagonalOp),
                "::",
                stringify!(deviceOperator)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_Qureg() {
        assert_eq!(
            ::std::mem::size_of::<Qureg>(),
            112usize,
            concat!("Size of: ", stringify!(Qureg))
        );
        assert_eq!(
            ::std::mem::align_of::<Qureg>(),
            8usize,
            concat!("Alignment of ", stringify!(Qureg))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).isDensityMatrix as *const _
                    as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(isDensityMatrix)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).numQubitsRepresented
                    as *const _ as usize
            },
            4usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(numQubitsRepresented)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).numQubitsInStateVec
                    as *const _ as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(numQubitsInStateVec)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).numAmpsPerChunk as *const _
                    as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(numAmpsPerChunk)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).numAmpsTotal as *const _
                    as usize
            },
            24usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(numAmpsTotal)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).chunkId as *const _ as usize
            },
            32usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(chunkId)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).numChunks as *const _ as usize
            },
            36usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(numChunks)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).stateVec as *const _ as usize
            },
            40usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(stateVec)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).pairStateVec as *const _
                    as usize
            },
            56usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(pairStateVec)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).deviceStateVec as *const _
                    as usize
            },
            72usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(deviceStateVec)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).firstLevelReduction
                    as *const _ as usize
            },
            88usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(firstLevelReduction)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).secondLevelReduction
                    as *const _ as usize
            },
            96usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(secondLevelReduction)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<Qureg>())).qasmLog as *const _ as usize
            },
            104usize,
            concat!(
                "Offset of field: ",
                stringify!(Qureg),
                "::",
                stringify!(qasmLog)
            )
        );
    }

    #[test]
    fn bindgen_test_layout_QuESTEnv() {
        assert_eq!(
            ::std::mem::size_of::<QuESTEnv>(),
            24usize,
            concat!("Size of: ", stringify!(QuESTEnv))
        );
        assert_eq!(
            ::std::mem::align_of::<QuESTEnv>(),
            8usize,
            concat!("Alignment of ", stringify!(QuESTEnv))
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QuESTEnv>())).rank as *const _ as usize
            },
            0usize,
            concat!(
                "Offset of field: ",
                stringify!(QuESTEnv),
                "::",
                stringify!(rank)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QuESTEnv>())).numRanks as *const _
                    as usize
            },
            4usize,
            concat!(
                "Offset of field: ",
                stringify!(QuESTEnv),
                "::",
                stringify!(numRanks)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QuESTEnv>())).seeds as *const _ as usize
            },
            8usize,
            concat!(
                "Offset of field: ",
                stringify!(QuESTEnv),
                "::",
                stringify!(seeds)
            )
        );
        assert_eq!(
            unsafe {
                &(*(::std::ptr::null::<QuESTEnv>())).numSeeds as *const _
                    as usize
            },
            16usize,
            concat!(
                "Offset of field: ",
                stringify!(QuESTEnv),
                "::",
                stringify!(numSeeds)
            )
        );
    }
}
