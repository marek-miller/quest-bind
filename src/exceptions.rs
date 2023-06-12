use std::{
    ffi::{
        c_char,
        CStr,
    },
    sync::{
        Arc,
        Mutex,
        OnceLock,
    },
};

use super::QuestError;

static QUEST_EXCEPTION_GUARD: OnceLock<Arc<Mutex<u8>>> = OnceLock::new();
static QUEST_EXCEPTION_ERROR: OnceLock<Arc<Mutex<Option<QuestError>>>> =
    OnceLock::new();

#[allow(non_snake_case)]
#[no_mangle]
unsafe extern "C" fn invalidQuESTInputError(
    errMsg: *const c_char,
    errFunc: *const c_char,
) {
    let err_msg = unsafe { CStr::from_ptr(errMsg) }.to_str().unwrap();
    let err_func = unsafe { CStr::from_ptr(errFunc) }.to_str().unwrap();

    let mut err = QUEST_EXCEPTION_ERROR
        .get_or_init(|| Arc::new(Mutex::new(None)))
        .lock()
        .unwrap();
    *err = Some(QuestError::InvalidQuESTInput {
        err_msg:  err_msg.to_owned(),
        err_func: err_func.to_owned(),
    });

    log::error!("QueST Error in function {err_func}: {err_msg}");
}

pub fn catch_quest_exception<T, F>(f: F) -> Result<T, QuestError>
where
    F: FnOnce() -> T,
{
    let _guard = QUEST_EXCEPTION_GUARD
        .get_or_init(|| Arc::new(Mutex::new(0)))
        .lock()
        .unwrap();
    let res = f();
    let err = {
        let mut err_lock = QUEST_EXCEPTION_ERROR
            .get_or_init(|| Arc::new(Mutex::new(None)))
            .lock()
            .unwrap();
        let err = err_lock.clone();
        *err_lock = None;
        err
    };

    if let Some(e) = err {
        Err(e)
    } else {
        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use std::thread;

    use super::*;
    use crate::create_complex_matrix_n;

    #[test]
    fn catch_exception_01() -> Result<(), QuestError> {
        let err = create_complex_matrix_n(0).unwrap_err();
        match err {
            QuestError::InvalidQuESTInput {
                ..
            } => Ok(()),
            _ => panic!(),
        }
    }

    #[test]
    fn catch_exception_02() {
        let _ = create_complex_matrix_n(1).unwrap();
    }

    #[test]
    fn catch_exception_parallel_01() {
        thread::scope(|s| {
            s.spawn(|| {
                catch_exception_01().unwrap();
                catch_exception_02();
            });
            s.spawn(|| {
                catch_exception_02();
                catch_exception_01().unwrap();
            });
            s.spawn(|| {
                catch_exception_01().unwrap();
                catch_exception_01().unwrap();
            });
            s.spawn(|| {
                catch_exception_02();
                catch_exception_02();
            });
        });
    }
}
