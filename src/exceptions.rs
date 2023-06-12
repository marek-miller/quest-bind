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

struct ExceptGuard(Arc<Mutex<u8>>);

impl Default for ExceptGuard {
    fn default() -> Self {
        Self(Arc::new(Mutex::new(0)))
    }
}

struct ExceptError(Arc<Mutex<Option<QuestError>>>);

impl Default for ExceptError {
    fn default() -> Self {
        Self(Arc::new(Mutex::new(None)))
    }
}

static QUEST_EXCEPTION_GUARD: OnceLock<ExceptGuard> = OnceLock::new();
static QUEST_EXCEPTION_ERROR: OnceLock<ExceptError> = OnceLock::new();

#[allow(non_snake_case)]
#[no_mangle]
unsafe extern "C" fn invalidQuESTInputError(
    errMsg: *const c_char,
    errFunc: *const c_char,
) {
    let err_msg = unsafe { CStr::from_ptr(errMsg) }.to_str().unwrap();
    let err_func = unsafe { CStr::from_ptr(errFunc) }.to_str().unwrap();

    let mut err = QUEST_EXCEPTION_ERROR
        .get_or_init(|| Default::default())
        .0
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
        .get_or_init(|| Default::default())
        .0
        .lock()
        .unwrap();
    let res = f();
    let err = {
        let mut err_lock = QUEST_EXCEPTION_ERROR
            .get_or_init(|| Default::default())
            .0
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

    use crate::create_complex_matrix_n;

    #[test]
    fn catch_exception_01() {
        let _ = create_complex_matrix_n(1).unwrap();
        let _ = create_complex_matrix_n(-1).unwrap_err();
    }

    #[test]
    fn catch_exception_parallel() {
        thread::scope(|s| {
            s.spawn(|| {
                catch_exception_01();
                catch_exception_01();
            });
            s.spawn(|| {
                catch_exception_01();
                catch_exception_01();
            });
            s.spawn(|| {
                catch_exception_01();
                catch_exception_01();
            });
            s.spawn(|| {
                catch_exception_01();
                catch_exception_01();
            });
        });
    }
}
