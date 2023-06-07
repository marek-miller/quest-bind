mod ffi;

// type Qreal = f64;

#[derive(Debug)]
pub struct QuESTEnv(ffi::QuESTEnv);

impl QuESTEnv {
    #[must_use]
    pub fn new() -> Self {
        unsafe { QuESTEnv(ffi::createQuESTEnv()) }
    }

    pub fn destroy(self) {
        unsafe { ffi::destroyQuESTEnv(self.0) }
    }

    /// Print report to standard output
    pub fn report(&self) {
        unsafe {
            ffi::reportQuESTEnv(self.0);
        }
    }
}

impl Default for QuESTEnv {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug)]
pub struct Qureg<'a> {
    env: &'a QuESTEnv,
    reg: ffi::Qureg,
}

impl<'a> Qureg<'a> {
    #[must_use]
    pub fn new(
        num_qubits: i32,
        env: &'a QuESTEnv,
    ) -> Self {
        unsafe {
            Qureg {
                env,
                reg: ffi::createQureg(num_qubits, env.0),
            }
        }
    }

    pub fn destroy(self) {
        unsafe { ffi::destroyQureg(self.reg, self.env.0) }
    }

    pub fn report(&self) {
        unsafe {
            ffi::reportState(self.reg);
        }
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
