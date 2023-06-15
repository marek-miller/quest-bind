extern crate cmake;

fn main() {
    let mut config = cmake::Config::new("QuEST");
    let config = config
        .no_build_target(true)
        .define("MULTITHREADED", "ON")
        .define("GPUACCELERATED", "OFF")
        .define("DISTRIBUTED", "OFF");

    #[cfg(feature = "gpu")]
    config
        .define("MULTITHREADED", "OFF")
        .define("GPUACCELERATED", "ON")
        .define("DISTRIBUTED", "OFF");

    #[cfg(feature = "mpi")]
    config
        .define("MULTITHREADED", "ON")
        .define("GPUACCELERATED", "OFF")
        .define("DISTRIBUTED", "ON");

    let dst = config.build();

    println!(
        "cargo:rustc-link-search=native={}/build/QuEST",
        dst.display()
    );
    println!("cargo:rustc-link-lib=dylib=QuEST");
}
