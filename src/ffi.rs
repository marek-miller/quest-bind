//! FFI bindings to libQuEST v3.5.0

#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::ffi::{
    c_char,
    c_double,
    c_int,
    c_longlong,
    c_ulong,
    CStr,
};

pub type qreal = c_double;

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
    pub real: *mut qreal,
    pub imag: *mut qreal,
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

    pub fn createQureg(
        numQubits: c_int,
        env: QuESTEnv,
    ) -> Qureg;

    pub fn createDensityQureg(
        numQubits: c_int,
        env: QuESTEnv,
    ) -> Qureg;

    pub fn createCloneQureg(
        qureg: Qureg,
        env: QuESTEnv,
    ) -> Qureg;

    pub fn destroyQureg(
        qureg: Qureg,
        env: QuESTEnv,
    );

    pub fn createComplexMatrixN(numQubits: c_int) -> ComplexMatrixN;

    pub fn destroyComplexMatrixN(matr: ComplexMatrixN);

    pub fn initComplexMatrixN(
        m: ComplexMatrixN,
        real: *const *const qreal,
        imag: *const *const qreal,
    );

    pub fn createPauliHamil(
        numQubits: c_int,
        numSumTerms: c_int,
    ) -> PauliHamil;

    pub fn destroyPauliHamil(hamil: PauliHamil);

    pub fn createPauliHamilFromFile(fn_: *const c_char) -> PauliHamil;

    pub fn initPauliHamil(
        hamil: PauliHamil,
        coeffs: *const qreal,
        codes: *const pauliOpType,
    );

    pub fn createDiagonalOp(
        numQubits: c_int,
        env: QuESTEnv,
    ) -> DiagonalOp;

    pub fn destroyDiagonalOp(
        op: DiagonalOp,
        env: QuESTEnv,
    );

    pub fn syncDiagonalOp(op: DiagonalOp);

    pub fn initDiagonalOp(
        op: DiagonalOp,
        real: *const qreal,
        imag: *const qreal,
    );

    pub fn initDiagonalOpFromPauliHamil(
        op: DiagonalOp,
        hamil: PauliHamil,
    );

    pub fn createDiagonalOpFromPauliHamilFile(
        fn_: *const c_char,
        env: QuESTEnv,
    ) -> DiagonalOp;

    pub fn setDiagonalOpElems(
        op: DiagonalOp,
        startInd: c_longlong,
        real: *mut qreal,
        imag: *mut qreal,
        numElems: c_longlong,
    );

    pub fn applyDiagonalOp(
        qureg: Qureg,
        op: DiagonalOp,
    );

    pub fn calcExpecDiagonalOp(
        qureg: Qureg,
        op: DiagonalOp,
    ) -> Complex;

    pub fn reportState(qureg: Qureg);

    pub fn reportStateToScreen(
        qureg: Qureg,
        env: QuESTEnv,
        reportRank: c_int,
    );

    pub fn reportQuregParams(qureg: Qureg);

    pub fn reportPauliHamil(hamil: PauliHamil);

    pub fn getNumQubits(qureg: Qureg) -> c_int;

    pub fn getNumAmps(qureg: Qureg) -> c_longlong;

    pub fn initBlankState(qureg: Qureg);

    pub fn initZeroState(qureg: Qureg);

    pub fn initPlusState(qureg: Qureg);

    pub fn initClassicalState(
        qureg: Qureg,
        stateInd: c_longlong,
    );

    pub fn initPureState(
        qureg: Qureg,
        pure_: Qureg,
    );

    pub fn initDebugState(qureg: Qureg);

    pub fn initStateFromAmps(
        qureg: Qureg,
        reals: *mut qreal,
        imags: *mut qreal,
    );

    pub fn setAmps(
        qureg: Qureg,
        startInd: c_longlong,
        reals: *mut qreal,
        imags: *mut qreal,
        numAmps: c_longlong,
    );

    pub fn setDensityAmps(
        qureg: Qureg,
        startRow: c_longlong,
        startCol: c_longlong,
        reals: *mut qreal,
        imags: *mut qreal,
        numAmps: c_longlong,
    );

    pub fn cloneQureg(
        targetQureg: Qureg,
        copyQureg: Qureg,
    );

    pub fn phaseShift(
        qureg: Qureg,
        targetQubit: c_int,
        angle: qreal,
    );

    pub fn controlledPhaseShift(
        qureg: Qureg,
        idQubit1: c_int,
        idQubit2: c_int,
        angle: qreal,
    );

    pub fn multiControlledPhaseShift(
        qureg: Qureg,
        controlQubits: *mut c_int,
        numControlQubits: c_int,
        angle: qreal,
    );

    pub fn controlledPhaseFlip(
        qureg: Qureg,
        idQubit1: c_int,
        idQubit2: c_int,
    );

    pub fn multiControlledPhaseFlip(
        qureg: Qureg,
        controlQubits: *mut c_int,
        numControlQubits: c_int,
    );

    pub fn sGate(
        qureg: Qureg,
        targetQubit: c_int,
    );

    pub fn tGate(
        qureg: Qureg,
        targetQubit: c_int,
    );

    pub fn createQuESTEnv() -> QuESTEnv;

    pub fn destroyQuESTEnv(env: QuESTEnv);

    pub fn syncQuESTEnv(env: QuESTEnv);

    pub fn syncQuESTSuccess(successCode: c_int) -> c_int;

    pub fn reportQuESTEnv(env: QuESTEnv);

    pub fn getEnvironmentString(
        env: QuESTEnv,
        str: *mut c_char,
    );

    pub fn copyStateToGPU(qureg: Qureg);

    pub fn copyStateFromGPU(qureg: Qureg);

    pub fn copySubstateToGPU(
        qureg: Qureg,
        startInd: c_longlong,
        numAmps: c_longlong,
    );

    pub fn copySubstateFromGPU(
        qureg: Qureg,
        startInd: c_longlong,
        numAmps: c_longlong,
    );

    pub fn getAmp(
        qureg: Qureg,
        index: c_longlong,
    ) -> Complex;

    pub fn getRealAmp(
        qureg: Qureg,
        index: c_longlong,
    ) -> qreal;

    pub fn getImagAmp(
        qureg: Qureg,
        index: c_longlong,
    ) -> qreal;

    pub fn getProbAmp(
        qureg: Qureg,
        index: c_longlong,
    ) -> qreal;

    pub fn getDensityAmp(
        qureg: Qureg,
        row: c_longlong,
        col: c_longlong,
    ) -> Complex;

    pub fn calcTotalProb(qureg: Qureg) -> qreal;

    pub fn compactUnitary(
        qureg: Qureg,
        targetQubit: c_int,
        alpha: Complex,
        beta: Complex,
    );

    pub fn unitary(
        qureg: Qureg,
        targetQubit: c_int,
        u: ComplexMatrix2,
    );

    pub fn rotateX(
        qureg: Qureg,
        rotQubit: c_int,
        angle: qreal,
    );

    pub fn rotateY(
        qureg: Qureg,
        rotQubit: c_int,
        angle: qreal,
    );

    pub fn rotateZ(
        qureg: Qureg,
        rotQubit: c_int,
        angle: qreal,
    );

    pub fn rotateAroundAxis(
        qureg: Qureg,
        rotQubit: c_int,
        angle: qreal,
        axis: Vector,
    );

    pub fn controlledRotateX(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
        angle: qreal,
    );

    pub fn controlledRotateY(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
        angle: qreal,
    );

    pub fn controlledRotateZ(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
        angle: qreal,
    );

    pub fn controlledRotateAroundAxis(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
        angle: qreal,
        axis: Vector,
    );

    pub fn controlledCompactUnitary(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
        alpha: Complex,
        beta: Complex,
    );

    pub fn controlledUnitary(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
        u: ComplexMatrix2,
    );

    pub fn multiControlledUnitary(
        qureg: Qureg,
        controlQubits: *mut c_int,
        numControlQubits: c_int,
        targetQubit: c_int,
        u: ComplexMatrix2,
    );

    pub fn pauliX(
        qureg: Qureg,
        targetQubit: c_int,
    );

    pub fn pauliY(
        qureg: Qureg,
        targetQubit: c_int,
    );

    pub fn pauliZ(
        qureg: Qureg,
        targetQubit: c_int,
    );

    pub fn hadamard(
        qureg: Qureg,
        targetQubit: c_int,
    );

    pub fn controlledNot(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
    );

    pub fn multiControlledMultiQubitNot(
        qureg: Qureg,
        ctrls: *mut c_int,
        numCtrls: c_int,
        targs: *mut c_int,
        numTargs: c_int,
    );

    pub fn multiQubitNot(
        qureg: Qureg,
        targs: *mut c_int,
        numTargs: c_int,
    );

    pub fn controlledPauliY(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit: c_int,
    );

    pub fn calcProbOfOutcome(
        qureg: Qureg,
        measureQubit: c_int,
        outcome: c_int,
    ) -> qreal;

    pub fn calcProbOfAllOutcomes(
        outcomeProbs: *mut qreal,
        qureg: Qureg,
        qubits: *mut c_int,
        numQubits: c_int,
    );

    pub fn collapseToOutcome(
        qureg: Qureg,
        measureQubit: c_int,
        outcome: c_int,
    ) -> qreal;

    pub fn measure(
        qureg: Qureg,
        measureQubit: c_int,
    ) -> c_int;

    pub fn measureWithStats(
        qureg: Qureg,
        measureQubit: c_int,
        outcomeProb: *mut qreal,
    ) -> c_int;

    pub fn calcInnerProduct(
        bra: Qureg,
        ket: Qureg,
    ) -> Complex;

    pub fn calcDensityInnerProduct(
        rho1: Qureg,
        rho2: Qureg,
    ) -> qreal;

    pub fn seedQuESTDefault(env: *mut QuESTEnv);

    pub fn seedQuEST(
        env: *mut QuESTEnv,
        seedArray: *mut c_ulong,
        numSeeds: c_int,
    );

    pub fn getQuESTSeeds(
        env: QuESTEnv,
        seeds: *mut *mut c_ulong,
        numSeeds: *mut c_int,
    );

    pub fn startRecordingQASM(qureg: Qureg);

    pub fn stopRecordingQASM(qureg: Qureg);

    pub fn clearRecordedQASM(qureg: Qureg);

    pub fn printRecordedQASM(qureg: Qureg);

    pub fn writeRecordedQASMToFile(
        qureg: Qureg,
        filename: *mut c_char,
    );

    pub fn mixDephasing(
        qureg: Qureg,
        targetQubit: c_int,
        prob: qreal,
    );

    pub fn mixTwoQubitDephasing(
        qureg: Qureg,
        qubit1: c_int,
        qubit2: c_int,
        prob: qreal,
    );

    pub fn mixDepolarising(
        qureg: Qureg,
        targetQubit: c_int,
        prob: qreal,
    );

    pub fn mixDamping(
        qureg: Qureg,
        targetQubit: c_int,
        prob: qreal,
    );

    pub fn mixTwoQubitDepolarising(
        qureg: Qureg,
        qubit1: c_int,
        qubit2: c_int,
        prob: qreal,
    );

    pub fn mixPauli(
        qureg: Qureg,
        targetQubit: c_int,
        probX: qreal,
        probY: qreal,
        probZ: qreal,
    );

    pub fn mixDensityMatrix(
        combineQureg: Qureg,
        prob: qreal,
        otherQureg: Qureg,
    );

    pub fn calcPurity(qureg: Qureg) -> qreal;

    pub fn calcFidelity(
        qureg: Qureg,
        pureState: Qureg,
    ) -> qreal;

    pub fn swapGate(
        qureg: Qureg,
        qubit1: c_int,
        qubit2: c_int,
    );

    pub fn sqrtSwapGate(
        qureg: Qureg,
        qb1: c_int,
        qb2: c_int,
    );

    pub fn multiStateControlledUnitary(
        qureg: Qureg,
        controlQubits: *mut c_int,
        controlState: *mut c_int,
        numControlQubits: c_int,
        targetQubit: c_int,
        u: ComplexMatrix2,
    );

    pub fn multiRotateZ(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubits: c_int,
        angle: qreal,
    );

    pub fn multiRotatePauli(
        qureg: Qureg,
        targetQubits: *mut c_int,
        targetPaulis: *mut pauliOpType,
        numTargets: c_int,
        angle: qreal,
    );

    pub fn multiControlledMultiRotateZ(
        qureg: Qureg,
        controlQubits: *mut c_int,
        numControls: c_int,
        targetQubits: *mut c_int,
        numTargets: c_int,
        angle: qreal,
    );

    pub fn multiControlledMultiRotatePauli(
        qureg: Qureg,
        controlQubits: *mut c_int,
        numControls: c_int,
        targetQubits: *mut c_int,
        targetPaulis: *mut pauliOpType,
        numTargets: c_int,
        angle: qreal,
    );

    pub fn calcExpecPauliProd(
        qureg: Qureg,
        targetQubits: *mut c_int,
        pauliCodes: *mut pauliOpType,
        numTargets: c_int,
        workspace: Qureg,
    ) -> qreal;

    pub fn calcExpecPauliSum(
        qureg: Qureg,
        allPauliCodes: *mut pauliOpType,
        termCoeffs: *mut qreal,
        numSumTerms: c_int,
        workspace: Qureg,
    ) -> qreal;

    pub fn calcExpecPauliHamil(
        qureg: Qureg,
        hamil: PauliHamil,
        workspace: Qureg,
    ) -> qreal;

    pub fn twoQubitUnitary(
        qureg: Qureg,
        targetQubit1: c_int,
        targetQubit2: c_int,
        u: ComplexMatrix4,
    );

    pub fn controlledTwoQubitUnitary(
        qureg: Qureg,
        controlQubit: c_int,
        targetQubit1: c_int,
        targetQubit2: c_int,
        u: ComplexMatrix4,
    );

    pub fn multiControlledTwoQubitUnitary(
        qureg: Qureg,
        controlQubits: *mut c_int,
        numControlQubits: c_int,
        targetQubit1: c_int,
        targetQubit2: c_int,
        u: ComplexMatrix4,
    );

    pub fn multiQubitUnitary(
        qureg: Qureg,
        targs: *mut c_int,
        numTargs: c_int,
        u: ComplexMatrixN,
    );

    pub fn controlledMultiQubitUnitary(
        qureg: Qureg,
        ctrl: c_int,
        targs: *mut c_int,
        numTargs: c_int,
        u: ComplexMatrixN,
    );

    pub fn multiControlledMultiQubitUnitary(
        qureg: Qureg,
        ctrls: *mut c_int,
        numCtrls: c_int,
        targs: *mut c_int,
        numTargs: c_int,
        u: ComplexMatrixN,
    );

    pub fn mixKrausMap(
        qureg: Qureg,
        target: c_int,
        ops: *mut ComplexMatrix2,
        numOps: c_int,
    );

    pub fn mixTwoQubitKrausMap(
        qureg: Qureg,
        target1: c_int,
        target2: c_int,
        ops: *mut ComplexMatrix4,
        numOps: c_int,
    );

    pub fn mixMultiQubitKrausMap(
        qureg: Qureg,
        targets: *mut c_int,
        numTargets: c_int,
        ops: *mut ComplexMatrixN,
        numOps: c_int,
    );

    pub fn mixNonTPKrausMap(
        qureg: Qureg,
        target: c_int,
        ops: *mut ComplexMatrix2,
        numOps: c_int,
    );

    pub fn mixNonTPTwoQubitKrausMap(
        qureg: Qureg,
        target1: c_int,
        target2: c_int,
        ops: *mut ComplexMatrix4,
        numOps: c_int,
    );

    pub fn mixNonTPMultiQubitKrausMap(
        qureg: Qureg,
        targets: *mut c_int,
        numTargets: c_int,
        ops: *mut ComplexMatrixN,
        numOps: c_int,
    );

    pub fn calcHilbertSchmidtDistance(
        a: Qureg,
        b: Qureg,
    ) -> qreal;

    pub fn setWeightedQureg(
        fac1: Complex,
        qureg1: Qureg,
        fac2: Complex,
        qureg2: Qureg,
        facOut: Complex,
        out: Qureg,
    );

    pub fn applyPauliSum(
        inQureg: Qureg,
        allPauliCodes: *mut pauliOpType,
        termCoeffs: *mut qreal,
        numSumTerms: c_int,
        outQureg: Qureg,
    );

    pub fn applyPauliHamil(
        inQureg: Qureg,
        hamil: PauliHamil,
        outQureg: Qureg,
    );

    pub fn applyTrotterCircuit(
        qureg: Qureg,
        hamil: PauliHamil,
        time: qreal,
        order: c_int,
        reps: c_int,
    );

    pub fn applyMatrix2(
        qureg: Qureg,
        targetQubit: c_int,
        u: ComplexMatrix2,
    );

    pub fn applyMatrix4(
        qureg: Qureg,
        targetQubit1: c_int,
        targetQubit2: c_int,
        u: ComplexMatrix4,
    );

    pub fn applyMatrixN(
        qureg: Qureg,
        targs: *mut c_int,
        numTargs: c_int,
        u: ComplexMatrixN,
    );

    pub fn applyMultiControlledMatrixN(
        qureg: Qureg,
        ctrls: *mut c_int,
        numCtrls: c_int,
        targs: *mut c_int,
        numTargs: c_int,
        u: ComplexMatrixN,
    );

    pub fn bindArraysToStackComplexMatrixN(
        numQubits: c_int,
        re: *mut *mut qreal,
        im: *mut *mut qreal,
        reStorage: *mut *mut qreal,
        imStorage: *mut *mut qreal,
    ) -> ComplexMatrixN;

    pub fn applyPhaseFunc(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubits: c_int,
        encoding: bitEncoding,
        coeffs: *mut qreal,
        exponents: *mut qreal,
        numTerms: c_int,
    );

    pub fn applyPhaseFuncOverrides(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubits: c_int,
        encoding: bitEncoding,
        coeffs: *mut qreal,
        exponents: *mut qreal,
        numTerms: c_int,
        overrideInds: *mut c_longlong,
        overridePhases: *mut qreal,
        numOverrides: c_int,
    );

    pub fn applyMultiVarPhaseFunc(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubitsPerReg: *mut c_int,
        numRegs: c_int,
        encoding: bitEncoding,
        coeffs: *mut qreal,
        exponents: *mut qreal,
        numTermsPerReg: *mut c_int,
    );

    pub fn applyMultiVarPhaseFuncOverrides(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubitsPerReg: *mut c_int,
        numRegs: c_int,
        encoding: bitEncoding,
        coeffs: *mut qreal,
        exponents: *mut qreal,
        numTermsPerReg: *mut c_int,
        overrideInds: *mut c_longlong,
        overridePhases: *mut qreal,
        numOverrides: c_int,
    );

    pub fn applyNamedPhaseFunc(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubitsPerReg: *mut c_int,
        numRegs: c_int,
        encoding: bitEncoding,
        functionNameCode: phaseFunc,
    );

    pub fn applyNamedPhaseFuncOverrides(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubitsPerReg: *mut c_int,
        numRegs: c_int,
        encoding: bitEncoding,
        functionNameCode: phaseFunc,
        overrideInds: *mut c_longlong,
        overridePhases: *mut qreal,
        numOverrides: c_int,
    );

    pub fn applyParamNamedPhaseFunc(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubitsPerReg: *mut c_int,
        numRegs: c_int,
        encoding: bitEncoding,
        functionNameCode: phaseFunc,
        params: *mut qreal,
        numParams: c_int,
    );

    pub fn applyParamNamedPhaseFuncOverrides(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubitsPerReg: *mut c_int,
        numRegs: c_int,
        encoding: bitEncoding,
        functionNameCode: phaseFunc,
        params: *mut qreal,
        numParams: c_int,
        overrideInds: *mut c_longlong,
        overridePhases: *mut qreal,
        numOverrides: c_int,
    );

    pub fn applyFullQFT(qureg: Qureg);

    pub fn applyQFT(
        qureg: Qureg,
        qubits: *mut c_int,
        numQubits: c_int,
    );

    pub fn applyProjector(
        qureg: Qureg,
        qubit: c_int,
        outcome: c_int,
    );

}

// TODO: Allow user to override this.
// For now, it is just a direct translation from
// QuEST's default callback.
#[no_mangle]
pub extern "C" fn invalidQuESTInputError(
    errMsg: *const c_char,
    errFunc: *const c_char,
) {
    let err_msg = unsafe { CStr::from_ptr(errMsg) }.to_str().unwrap();
    let err_func = unsafe { CStr::from_ptr(errFunc) }.to_str().unwrap();

    println!("!!!");
    println!("QueST Error in function {err_func}: {err_msg}");
    println!("!!!");
    println!("exiting..");
    panic!();
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

    #[test]
    #[should_panic]
    fn invalid_quest_input_error_throw_01() {
        unsafe {
            createComplexMatrixN(0);
        };
    }
}
