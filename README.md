# quest-bind

A safe wrapper around [QuEST](https://github.com/QuEST-Kit/QuEST/) v3.5.0. As
thin as possible: the API stays almost identical to the original.

QuEST (Quantum Exact Simulation Toolkit) is a no-fluff, bent-on-speed quantum
circuit simulator. It is distributed under MIT License.

## Testing

In order to test this library, first clone the repository together with QuEST
source code as submodule:

```sh
git clone --recurse-submodules https://github.com/marek-miller/quest-bind.git
```

Then compile QuEST _without_ MPI support:

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

Now, we can run unit tests:

```sh
cargo test
```

## How to use it

No package management for now. Here's a brief tutorial how to compile everything
manually.

Initialize a new binary crate:

```sh
cargo new testme
cd testme/
```

Clone the dependencies as before:

```sh
git clone --recurse-submodules https://github.com/marek-miller/quest-bind.git
```

Compile QuEST (you can enable MPI this time) and put a link to the shared
library in your crate's `./target/debug/deps`. Then add dependencies to your
project's `Cargo.toml`:

```toml
[package]
name = "testme"
version = "0.1.0"
edition = "2021"

[dependencies]
quest_bind = { path = "./quest-bind" }
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
        println!("Outcome: {} with prob: {:.2}", outcome, outcome_prob);
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

## Releases

### v0.2.0 (??/07/2023)

New features/improvements:

- Catch exceptions thrown by QuEST
- Improve documentation

### v0.1.0 (11/06/2023)

Initial release.
