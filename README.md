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
    let env = create_quest_env();
    report_quest_env(&env);

    let mut qureg = create_qureg(0x10, &env)?;
    {
        let qureg = &mut qureg;
        init_plus_state(qureg);
        report_qureg_params(qureg);

        let mut outcome_prob = 0.;
        let outcome = measure_with_stats(qureg, 1, &mut outcome_prob);

        println!("Measure first qubit.");
        println!("Outcome: {outcome} with prob: {outcome_prob:.2}");
    }

    destroy_qureg(qureg, &env);
    destroy_quest_env(env);
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

## Testing

To run unit tests for this library, first clone the repository together with QuEST
source code as submodule:

```sh
git clone --recurse-submodules https://github.com/marek-miller/quest-bind.git
```

Then run:

```sh
cargo test --tests
```

To be able to run documentation tests, we need to work around a known issue with `rustdoc`: [#8531](https://github.com/rust-lang/cargo/issues/8531).

Make sure you compile QuEST _without_ MPI support:

```sh
cd quest-bind
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

## Releases

### v0.2.0 (??/07/2023)

New features/improvements:

- Catch exceptions thrown by QuEST
- Improve documentation
- Add build script

### v0.1.0 (11/06/2023)

Initial release.
