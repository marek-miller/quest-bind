#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

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
pub struct QuESTEnv {
    rank:     c_int,
    numRanks: c_int,
    seeds:    *mut c_ulong,
    numSeeds: c_int,
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
