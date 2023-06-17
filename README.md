# quest-bind

A safe wrapper around [QuEST](https://github.com/QuEST-Kit/QuEST/) v3.5.0. As
thin as possible: the API stays nearly identical to the original.

QuEST (Quantum Exact Simulation Toolkit) is a no-fluff, bent-on-speed quantum
circuit simulator. It is distributed under MIT License.

## How to use it

Initialize a new binary crate:

```sh
cargo new testme
cd testme/
```

Add `quest-bind` to your project's dependencies

```sh
cargo add quest_bind --git https://github.com/marek-miller/quest-bind.git
```

Now write some code and put it in `./src/main.rs`:

```rust
use quest_bind::*;


fn main() -> Result<(), QuestError> {
    let env = &QuestEnv::new();
    report_quest_env(env);

    let qureg = &mut Qureg::try_new(0x10, env)?;
    report_qureg_params(qureg);
    init_plus_state(qureg);

    let outcome_prob = &mut -1.;
    let outcome = measure_with_stats(qureg, 0, outcome_prob)?;
    println!("Measure first qubit.");
    println!("Outcome: {outcome} with prob.: {outcome_prob:.2}");

    Ok(())
}
```

You can read the available documentation locally (refer to
[QuEST headers](https://github.com/QuEST-Kit/QuEST/blob/v3.5.0/QuEST/include/QuEST.h)
for the full description of the C API):

```sh
cargo doc --open
```

Lastly, compile the source code and run:

```sh
cargo run
```

You should be able to see something like this:

```text
EXECUTION ENVIRONMENT:
Running locally on one node
Number of ranks is 1
OpenMP enabled
Number of threads available is 8
Precision: size of qreal is 8 bytes
QUBITS:
Number of qubits is 16.
Number of amps is 65536.
Number of amps per rank is 65536.
Measure first qubit.
Outcome: 0 with prob: 0.50
```

## Distributed and GPU-accelerated mode

QuEST support for MPI and GPU accelerated computation is enabled by setting
appropriate feature flags. To use QuEST's MPI mode, enable `mpi` feature for
`quest_bind`. Simply edit `Cargo.toml` of your binary crate:

```toml
[package]
name = "testme"
version = "0.1.0"
edition = "2021"

[dependencies]
quest_bind = { git = "https://github.com/marek-miller/quest-bind.git", features = ["mpi",] }
```

Now if you compile and run the above program again, the output should be:

```text
EXECUTION ENVIRONMENT:
Running distributed (MPI) version
Number of ranks is 1
...
```

The feature `"gpu"` enables the GPU accelerated mode. These features are
mutually exclusive and in case both flags are set, the feature `"mpi"` takes
precedence.

## Testing

To run unit tests for this library, first clone the repository together with
QuEST source code as submodule:

```sh
git clone --recurse-submodules https://github.com/marek-miller/quest-bind.git
cd quest-bind
```

Then run:

```sh
cargo test
```

Please note that `quest-bin` will not run `QuEST`'s test suite, not will it
check `QuEST`'s correctness. These tests are intended to check if the C API is
invoked correctly, and Rust's types are passed safely back and forth across FFI
boundaries.

If you want to run the test suite in the single-precision floating point mode,
make sure the build script recompiles `libQuEST.so` with the right type
definitions:

```sh
cargo clean
cargo test --features=f32
```

By defualt, `quest-bind` uses double precision floating-point types: `f64`. See
[Numercal types](#numerical-types) section below.

You can also check available examples by running, e.g.:

```sh
 cargo run --release --example grovers_search
```

To see all available examples, try:

```sh
cargo run --example
```

## Note on performance

In the typical case when it's the numerical computation that dominates the CPU
usage, and not API calls, there should be no discernible difference in
performance between programs calling QuEST routines directly and analogous
applications using `quest_bind`. Remember, however, to enable compiler
optimizations by running your code using the `--release` profile:

```sh
cargo run --release
```

## Handling exceptions

On failure, QuEST throws exceptions via user-configurable global
[`invalidQuESTInputError()`](https://quest-kit.github.io/QuEST/group__debug.html#ga51a64b05d31ef9bcf6a63ce26c0092db).
By default, this function prints an error message and aborts, which is
problematic in a large distributed setup.

We opt for catching all exceptions early. The exception handler is locked during
an API call. This means that calling QuEST functions is synchronous and should
be thread-safe, but comes at the expense of being able to run only one QuEST API
call at the time. Bear in mind, though, that each QuEST function retains access
to all parallel computation resources available in the system.

Current implementation returns inside `Result<_, QuestError>` only the first
exception caught. All subsequent messages reported by QuEST, together with that
first one, are nevertheless logged as errors. To review them, add a logger as a
dependency to your crate, e.g.:

```sh
cargo add env_logger
```

Then enable logging in your application:

```rust
fn main()  {
    env_logger::init();

    // (...)

}
```

and run:

```sh
RUST_LOG=info cargo run
```

See [`log` crate](https://docs.rs/log/latest/log/) for more on logging in Rust.

The type `QuestError` doesn't contain (possibly malformed) data returned by the
API call on failure. Only successful calls can reach the library user. This is
intentional, following guidelines in QuEST documentation. Upon failure...

> Users must ensure that the triggered API call does not continue (e.g. the user
> exits or throws an exception), else QuEST will continue with the valid [sic!]
> input and likely trigger a seg-fault.

See
[Quest API](https://quest-kit.github.io/QuEST/group__debug.html#ga51a64b05d31ef9bcf6a63ce26c0092db)
for more information.

## Numerical types

For now, numerical types used by `quest-bind` are chosen to match exactly the C
types that QuEST uses on `x86_64`. This is a safe, but not very portable
strategy. We pass Rust types directly to QuEST without casting, assuming the
following type definitions:

```rust
pub type c_float = f32;
pub type c_double = f64;
pub type c_int = i32;
pub type c_longlong = i64;
pub type c_ulong = u64;
```

This should work for many different architectures. If your system uses slightly
different numerical types, `quest-bind` simply won't compile and there isn't
much you can do besides manually altering the source code. In the future,
`quest-bind` will use numerical abstractions from
[`num_traits`](https://docs.rs/num-traits/latest/num_traits/) for greater
portability.

To check what C types are defined by your Rust installation, see local
documentation for the module `std::ffi` in Rust's Standard Library:

```sh
rustup doc
```

## Releases

### v0.2.0 (??/07/2023)

New features/improvements:

- Catch exceptions thrown by QuEST
- Improve documentation
- Add build script
- Constructors/destructors for QuEST structs
- Add Github workflows CT
- Add example: grovers_search.rs
- Use `Complex<f64>` type from `num` crate (as `QComplex`)
- Use compile flag `"f32"` to set floating point precision

### v0.1.0 (11/06/2023)

Initial release.
