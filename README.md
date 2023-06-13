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
    let env = QuESTEnv::new();
    report_quest_env(&env);

    let mut qureg = Qureg::try_new(0x10, &env)?;
    {
        let qureg = &mut qureg;
        init_plus_state(qureg);
        report_qureg_params(qureg);

        let mut outcome_prob = 0.;
        let outcome = measure_with_stats(qureg, 1, &mut outcome_prob)?;

        println!("Measure first qubit.");
        println!("Outcome: {outcome} with prob: {outcome_prob:.2}");
    }

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

## Distributed and GPU accelerated mode

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

Now if we compile and run the above programme again, the output should be:

```text
EXECUTION ENVIRONMENT:
Running distributed (MPI) version
Number of ranks is 1
...
```

The feature `"gpu"` enables the GPU accelerated mode. These features are
mutually exclusive, so if you set both flags, the feature `"mpi"` takes
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
cargo test --tests
```

To be able to run documentation tests, we need to work around a known issue with
`cargo`: [#8531](https://github.com/rust-lang/cargo/issues/8531).

Make sure you compile QuEST _without_ MPI support:

```sh
cd QuEST
mkdir build
cd build
cmake -DDISTRIBUTED=OFF ..
make
```

Put a symlink to `libQuEST.so` where cargo can find it:

```sh
cd ../../
mkdir -p target/debug/deps
cd target/debug/deps/
ln -s ../../../QuEST/build/QuEST/libQuEST.so .
cd ../../../
```

Now, we can run all unit tests:

```sh
cargo test
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

Then enable logging in your applications:

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

```text
Users must ensure that the triggered API call does not continue (e.g. the user exits or throws an exception), else QuEST will continue with the valid [sic!] input and likely trigger a seg-fault.
```

See
[Quest API](https://quest-kit.github.io/QuEST/group__debug.html#ga51a64b05d31ef9bcf6a63ce26c0092db)
for more information.

## Releases

### v0.2.0 (??/07/2023)

New features/improvements:

- Catch exceptions thrown by QuEST
- Improve documentation
- Add build script
- Constructors/destructors for QuEST structs

### v0.1.0 (11/06/2023)

Initial release.
