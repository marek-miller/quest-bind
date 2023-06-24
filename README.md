# quest-bind

[![Test](https://github.com/marek-miller/quest-bind/actions/workflows/test.yml/badge.svg)](https://github.com/marek-miller/quest-bind/actions/workflows/test.yml)
[![Docs](https://github.com/marek-miller/quest-bind/actions/workflows/docs.yml/badge.svg)](https://github.com/marek-miller/quest-bind/actions/workflows/docs.yml)

A wrapper around [QuEST](https://github.com/QuEST-Kit/QuEST/) v3.5.0.

Quantum Exact Simulation Toolkit (QuEST) is a no-fluff, bent-on-speed quantum
circuit simulator [[1]](https://doi.org/10.1038/s41598-019-47174-9). It is
distributed under MIT License.

## How to use it

Initialize a new binary crate:

```sh
cargo new testme
cd testme/
```

Add `quest-bind` to your project's dependencies:

```sh
cargo add quest_bind --git https://github.com/marek-miller/quest-bind.git
```

Now write some code and put it in `./src/main.rs`:

```rust
use quest_bind::*;

fn main() -> Result<(), QuestError> {
    // Initialize QuEST environment
    let env = &QuestEnv::new();
    report_quest_env(env);

    // Initialize 2-qubit register in |00> state
    let qureg = &mut Qureg::try_new(2, env)?;
    report_qureg_params(qureg);
    init_zero_state(qureg);

    println!("---\nPrepare Bell state: |00> + |11>");
    hadamard(qureg, 0).and(controlled_not(qureg, 0, 1))?;

    // Measure each qubit
    let outcome0 = measure(qureg, 0)?;
    let outcome1 = measure(qureg, 1)?;

    println!("Qubit \"0\" measured in state: |{outcome0}>");
    println!("Qubit \"1\" measured in state: |{outcome1}>");
    assert_eq!(outcome0, outcome1);
    println!("They match!");

    // Both `env` and `qureg` are safely discarded here
    Ok(())
}
```

You can read the available documentation locally or refer to
[QuEST headers](https://github.com/QuEST-Kit/QuEST/blob/v3.5.0/QuEST/include/QuEST.h)
for the full description of the C API:

```sh
cargo doc --open
```

Lastly, compile and run the program:

```sh
cargo run
```

You should be able to see something like:

```text
EXECUTION ENVIRONMENT:
Running locally on one node
Number of ranks is 1
OpenMP enabled
Number of threads available is 8
Precision: size of qreal is 8 bytes
QUBITS:
Number of qubits is 2.
Number of amps is 4.
Number of amps per rank is 4.
---
Prepare Bell state: |00> + |11>
Qubit "0" measured in state: |0>
Qubit "1" measured in state: |0>
They match!
```

## Distributed and GPU-accelerated mode

QuEST support for MPI and GPU-accelerated computation ca be enabled in
`quest-bind` by setting appropriate feature flags. To enable QuEST's MPI mode,
set the `mpi` feature for `quest_bind`. Simply edit `Cargo.toml` of your binary
crate:

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

The feature `"gpu"` enables the GPU-accelerated mode. These two features are
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

Note that
`quest-bind will not run `QuEST`'s test suite, nor will it check `QuEST`'s
correctness. The tests here are intended to check if the C API is invoked
correctly, and if Rust's types are passed safely back and forth across the FFI
boundary.

If you want to run the test suite in the single-precision floating point mode,
make sure the build script recompiles `libQuEST.so` with the right type
definitions:

```sh
cargo clean
cargo test --features=f32
```

By defualt, `quest-bind` uses Rust's double precision floating-point type:
`f64`. See [Numercal types](#numerical-types) section below.

You can also try the available examples by running, e.g.:

```sh
 cargo run --release --example grovers_search
```

To see the list of all available examples, try:

```sh
cargo run --example
```

## Note on performance

In the typical case when it's the numerical computation that dominates the CPU
usage, and not API calls, there should be no discernible difference in
performance between programs calling QuEST routines directly and analogous
applications using `quest_bind`. Remember, however, to enable optimizations for
both `quest-bind` and `QuEST` by compiling your code using the "release"
profile:

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
first one, are nevertheless logged as errors. To be able to see them, add a
logger as a dependency to your crate, e.g.:

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
intentional, following guidelines from the QuEST documentation:

> [*Upon failure*] Users must ensure that the triggered API call does not
> continue (e.g. the user exits or throws an exception), else QuEST will
> continue with the valid [*sic!*] input and likely trigger a seg-fault.

See
[Quest API](https://quest-kit.github.io/QuEST/group__debug.html#ga51a64b05d31ef9bcf6a63ce26c0092db)
for more information.

## Numerical types

For now, numerical types used by `quest-bind` match exactly the C types that
QuEST uses on `x86_64`. This is a safe, but not very portable strategy. We pass
Rust types directly to QuEST without casting, assuming the following type
definitions:

```rust
pub type c_float = f32;
pub type c_double = f64;
pub type c_int = i32;
pub type c_longlong = i64;
pub type c_ulong = u64;
```

This should work for many different architectures. If your system uses slightly
different numerical types, `quest-bind` simply won't compile and there is not
much you can do besides manually altering the source code.

To check what C types are defined by your Rust installation, see the local
documentation for the module `std::ffi` in Rust's Standard Library:

```sh
rustup doc
```

## TODO

- Expand and improve documentation.
- Expand test suite
- Design test for MPI and GPU modes
- Add non-blocking API: `*_nonblk()` functions.
- Publish to `crates.io`
- Generic nummerical traits from `num_traits`

## Releases

### v0.2.1 (??/07/2023)

New features/improvements:

### v0.2.0 (24/06/2023)

New features/improvements:

- Improve documentation
- Catch exceptions thrown by QuEST
- Add build script
- Constructors/destructors for QuEST structs
- Add example: `grovers_search.rs`
- Use `Complex<f64>` type from `num` crate (as `QComplex`)
- Use compile flag `"f32"` to set floating point precision
- Add Github workflows CT

### v0.1.0 (11/06/2023)

Initial release.
