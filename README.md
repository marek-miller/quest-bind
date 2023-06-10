# quest

A safe wrapper around [QuEST](https://github.com/QuEST-Kit/QuEST/) v3.5.0. As
thin as possible: the API stays almost identical to the original.

`QuEST` (Quantum Exact Simulation Toolkit) is a no-fluff, bent-on-speed quantum
circuit simulator. It is distributed under MIT License.

## Testing

In order to test this library, first clone the repository together with `QuEST`
source code as submodule:

```sh
git clone --recurse-submodules https://github.com/marek-miller/quest.git
```

Then compile `QuEST` _without_ MPI support:

```sh
cd quest
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

No packaging for now. Here's a brief tutorial how to compile everything
manually.

Initialize a new binary crate:

```sh
cargo new testme
cd testme/
```

Clone the dependencies:

```sh
git clone --recurse-submodules https://github.com/marek-miller/quest.git
```

Compile `QuEST` as above (you can enable MPI this time) and put a link to the
shared library in your crate's `./target/debug/deps`. Then add dependencies to
your project's `Cargo.toml`:

```toml
[package]
name = "testme"
version = "0.1.0"
edition = "2021"

[dependencies]
quest = { path = "./quest" }
```

Now write some code and put it in `./src/main.rs`:

```rust
use quest::*;

fn main() {
    let env = create_quest_env();
    report_quest_env(&env);

    let mut qureg = create_qureg(3, &env);
    init_plus_state(&mut qureg);
    report_qureg_params(&qureg);
    destroy_qureg(qureg, &env);

    destroy_quest_env(env);
}
```

You can read the available documentation locally:

```sh
cargo doc --open
```

Refer to [`QuEST` headers](./QuEST/QuEST/include/QuEST.h) for the full
description of the C API.

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
Number of qubits is 3.
Number of amps is 8.
Number of amps per rank is 8.
```
