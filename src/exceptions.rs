//! Catch exceptions thrown by QuEST.
//!
//! On failure, QuEST throws exceptions via user-configurable global
//! [`invalidQuESTInputError()`][1]. By default, this function prints an error
//! message and aborts, which is problematic in a large distributed setup. We
//! opt for catching all exceptions early. The exception handler is locked
//! during an API call. This means that calling QuEST functions is synchronous
//! and should be thread-safe, but comes at the expense of being able to run
//! only one QuEST API call at the time.
//!
//! This is an internal module that doesn't contain any useful user interface.
//!
//! [1]: https://quest-kit.github.io/QuEST/group__debug.html#ga51a64b05d31ef9bcf6a63ce26c0092db

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

#[derive(Default)]
struct ExceptGuard(Arc<Mutex<[u8; 0]>>);

#[derive(Default)]
struct ExceptError(Arc<Mutex<Option<QuestError>>>);

static QUEST_EXCEPT_GUARD: OnceLock<ExceptGuard> = OnceLock::new();
static QUEST_EXCEPT_ERROR: OnceLock<ExceptError> = OnceLock::new();

/// Report error in a QuEST API call.
///
/// This function is called by QuEST whenever an error occurs.  We redefine it
/// to put the error message and site reported into a thread-safe global
/// storage, instead of just aborting.
///
/// # Panics
///
/// This function will panic if strings returned by QuEST are not properly
/// formatted (null terminated) C strings, or if our mutex is poisoned.
#[allow(non_snake_case)]
#[no_mangle]
unsafe extern "C" fn invalidQuESTInputError(
    errMsg: *const c_char,
    errFunc: *const c_char,
) {
    let err_msg = unsafe { CStr::from_ptr(errMsg) }.to_str().unwrap();
    let err_func = unsafe { CStr::from_ptr(errFunc) }.to_str().unwrap();

    let mut err = QUEST_EXCEPT_ERROR
        .get_or_init(Default::default)
        .0
        .lock()
        .unwrap();
    // assert!(
    //     err.is_none(),
    //     "All exceptions must be dealt with. This is a bug in quest_bind.  \
    //     Please report it."
    // );
    *err = Some(QuestError::InvalidQuESTInputError {
        err_msg:  err_msg.to_owned(),
        err_func: err_func.to_owned(),
    });
    eprintln!("except {:?}", err);

    log::error!("QueST Error in function {err_func}: {err_msg}");
}

/// Execute a call to QuEST API and catch exceptions.
///
/// This function achieves synchronous execution between threads
/// by locking the global `QUEST_EXCEPTION_GUARD` each time,
/// then executing the closure supplied, and finally checking the global
/// storage `QUEST_EXCEPTION_ERROR` for any error messages reported downstream.
///
/// This way, interacting with QuEST API should stay thread-safe at all times,
/// at the expense of being able to call only one function at the time.
/// This is not an undesired property and shouldn't matter much for the overall
/// performance of the simulation, since each functions retains access to all
/// parallelism available in the system.
pub fn catch_quest_exception<T, F>(f: F) -> Result<T, QuestError>
where
    F: FnOnce() -> T,
{
    // Lock QuEST to our call
    let guard = QUEST_EXCEPT_GUARD
        .get_or_init(Default::default)
        .0
        .lock()
        .unwrap();

    // Call QuEST API
    let res = f();

    // Catch exception
    let err = {
        // std::thread::sleep(std::time::Duration::from_millis(100));
        // Wait for QueEST to finish reporting
        let mut err = QUEST_EXCEPT_ERROR
            .get_or_init(Default::default)
            .0
            .lock()
            .unwrap();

        // eprintln!("catch : {:?}", err);

        // This is important! Wipe out error message before the next API call
        (*err).take()

        // the lock to `err` is dropped here
    };

    // Drop the guard as soon as we don't need it anymore:
    drop(guard);

    // This might be a little confusing:
    // If there is Some error, report it;
    // if there is None, everything's Ok.
    match err {
        Some(e) => Err(e),
        None => Ok(res),
    }
}

#[cfg(test)]
mod tests {
    use std::thread;

    use crate::{
        ComplexMatrixN,
        PauliHamil,
    };

    #[test]
    fn catch_exception_01() {
        let _ = ComplexMatrixN::try_new(1).unwrap();
        // Seems like supplying other invalid params here, like e.g. -3,
        // causes QuEST to hang.  Or is this a bug on our side?
        let _ = ComplexMatrixN::try_new(0).unwrap_err();
    }

    #[test]
    fn catch_exception_02() {
        let _ = PauliHamil::try_new(2, 2).unwrap();
        let _ = PauliHamil::try_new(-11, -3).unwrap_err();
    }

    #[test]
    fn catch_exception_parallel_01() {
        thread::scope(|s| {
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

    #[test]
    fn catch_exception_parallel_02() {
        thread::scope(|s| {
            s.spawn(|| {
                catch_exception_02();
                catch_exception_02();
            });
            s.spawn(|| {
                catch_exception_02();
                catch_exception_02();
            });
        });
    }
}
