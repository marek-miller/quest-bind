use std::{
    ffi::{
        c_char,
        CStr,
    },
    sync::{
        Arc,
        Mutex,
    },
};

use lazy_static::lazy_static;

use super::Error;

lazy_static! {
    static ref QUEST_EXCEPTION: Arc<Mutex<Result<(), Error>>> =
        Arc::new(Mutex::new(Ok(())));
}

#[allow(non_snake_case)]
#[no_mangle]
unsafe extern "C" fn invalidQuESTInputError(
    errMsg: *const c_char,
    errFunc: *const c_char,
) {
    let err_msg = unsafe { CStr::from_ptr(errMsg) }.to_str().unwrap();
    let err_func = unsafe { CStr::from_ptr(errFunc) }.to_str().unwrap();

    let mut excep = QUEST_EXCEPTION.lock().unwrap();
    *excep = Err(Error::InvalidQuESTInput {
        err_msg:  err_msg.to_owned(),
        err_func: err_func.to_owned(),
    });

    log::error!("QueST Error in function {err_func}: {err_msg}");
}

pub fn catch_quest_exception<T>(value: T) -> Result<T, Error> {
    let mut guard = QUEST_EXCEPTION.lock().unwrap();
    let excep = guard.clone();
    match excep {
        Err(Error::InvalidQuESTInput {
            ..
        }) => {
            *guard = Ok(());
            excep.map(|_| value)
        }
        _ => Ok(value),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_complex_matrix_n;

    #[test]
    fn catch_exception_01() -> Result<(), Error> {
        let err = create_complex_matrix_n(0).unwrap_err();
        match err {
            Error::InvalidQuESTInput {
                ..
            } => Ok(()),
            _ => panic!(),
        }
    }

    #[test]
    fn catch_exception_02() -> Result<(), Error> {
        let err = create_complex_matrix_n(-1).unwrap_err();
        match err {
            Error::InvalidQuESTInput {
                ..
            } => Ok(()),
            _ => panic!(),
        }
    }
}
