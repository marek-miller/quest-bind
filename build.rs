extern crate cmake;

fn main() {
    let mut config = cmake::Config::new("QuEST");
    let mut config = config.no_build_target(true);

    config = config
        .define("MULTITHREADED", "ON")
        .define("GPUACCELERATED", "OFF")
        .define("DISTRIBUTED", "OFF");

    if cfg!(feature = "gpu") {
        config = config
            .define("MULTITHREADED", "OFF")
            .define("GPUACCELERATED", "ON")
            .define("DISTRIBUTED", "OFF");
    }

    if cfg!(feature = "mpi") {
        config = config
            .define("MULTITHREADED", "ON")
            .define("GPUACCELERATED", "OFF")
            .define("DISTRIBUTED", "ON");
    }

    let dst = config.build();

    println!(
        "cargo:rustc-link-search=native={}/build/QuEST",
        dst.display()
    );
    println!("cargo:rustc-link-lib=dylib=QuEST");
}
